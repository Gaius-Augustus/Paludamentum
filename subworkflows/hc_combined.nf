nextflow.enable.dsl=2

// Unified HC-gene step for evidence modes (rnaseq / isoseq / mixed) and for
// the protein-only mode (proteins).
//
// Replaces the three HC flows of the previous pipeline (TransDecoder,
// Drusilla with proteins, Drusilla without proteins). Runs the same chain
// for every clade:
//
//  1. Run both gene finders on the genome: Tiberius (passed in from main.nf,
//     already executed to feed PROTEIN_EVIDENCE), and Vipsania --finetune.
//  2. On the one StringTie assembly (STRINGTIE_MERGE passes it through and
//     writes the GFF3 for TD2), run both ORF callers:
//     - Drusilla (filter + annotate + collapse isoforms),
//     - TD2 with the partial-ORF rescue: SHORTEN_INCOMPLETE_ORFS + two
//       auxiliary 12-column DIAMOND runs + CDS_CLASSIFY_AND_REVISE, then
//       TD_REVISED_TO_GENOME maps the revised peptides to the genome.
//  3. Extract proteins for all four gene sets (gffread via PROTEIN_FROM_GFF)
//     and run DIAMOND_BLASTP_FULL (14-column) on each against the protein DB.
//  4. BAM2COV turns the merged RNA-Seq / Iso-Seq BAM into a bedGraph.
//  5. combine_genesets.py integrates the four gene sets with intron hints,
//     coverage, DIAMOND, miniprot and miniprothint into combined.gff3.
//     Full-length transcripts (the StringTie -L assembly of long reads, or
//     the files of params.chains) go in as --chains when there are any.
//     Selection is learned (params.combine.classifier, the shipped
//     conf/combiner_classifier_v1.json): a clade-agnostic model trained on
//     15 species outside the benchmark; clades without a classifier_clade in
//     the HC table get its conservative locus threshold (--unknown-clade).
//     params.combine.classifier = false switches back to the tier rules.
//     For Insecta / Lepidoptera clades, --drop-lone-single-exon tiberius
//     is added in both modes.
//
// Without transcript evidence (transcripts = false, mode proteins) steps 2
// and 4 are skipped: the candidates are the Tiberius and Vipsania gene sets
// alone, with their two DIAMOND tables and the miniprot / miniprothint
// evidence, and the empty placeholder stands for the ORF-caller gene sets,
// their DIAMOND tables, the intron hints, the coverage track and the chains.
// Selection then uses the protein-only classifier
// (params.combine.classifier_proteins, the shipped
// conf/combiner_classifier_protonly_v1.json, trained on the same 15 species
// on candidate pools without transcript evidence); false: the tier rules.

include { GENEFINDER as GENEFINDER_VIPSANIA } from './genefinder.nf'
include { STRINGTIE_MERGE } from '../modules/assembly.nf'
include { TD_ALL; SHORTEN_INCOMPLETE_ORFS; CDS_CLASSIFY_AND_REVISE;
          TD_REVISED_TO_GENOME } from '../modules/transdecoder.nf'
include { FILTER_STRINGTIE; DOWNLOAD_DRUSILLA_MODEL; DRUSILLA_ANNOTATE;
          COLLAPSE_ORFS } from '../modules/drusilla.nf'
include { PROTEIN_FROM_GFF as PROT_TIB; PROTEIN_FROM_GFF as PROT_VIP;
          PROTEIN_FROM_GFF as PROT_TD2; PROTEIN_FROM_GFF as PROT_DRU } from '../modules/genefinder.nf'
include { DIAMOND_MAKEDB;
          DIAMOND_BLASTP as DIAMOND_NORM; DIAMOND_BLASTP as DIAMOND_SHORT;
          DIAMOND_BLASTP_FULL as DIAMOND_TIB; DIAMOND_BLASTP_FULL as DIAMOND_VIP;
          DIAMOND_BLASTP_FULL as DIAMOND_TD2; DIAMOND_BLASTP_FULL as DIAMOND_DRU } from '../modules/diamond.nf'
include { BAM2COV } from '../modules/rnacov.nf'
include { COMBINE_GENESETS } from '../modules/combine.nf'
include { EMPTY_FILE as NO_CLASSIFIER; EMPTY_FILE as NO_TRANSCRIPTS } from '../modules/util.nf'
include { hcClade; tiberiusModelValue; drusillaSetting; truthy } from '../lib_nf/functions.nf'

workflow HC_COMBINED {

    take:
    CH_GENOME
    asm_gtf                // channel: one or more StringTie assembly GTFs (empty without transcripts)
    merged_bam             // channel: one coordinate-sorted BAM of all RNA / Iso-Seq reads, or an empty placeholder
    have_bam               // value: false without a BAM (pyVARUS directories, a given assembly without reads, mode proteins)
    intron_hints           // channel: the intron hints of BAM2HINTS (two files in mixed mode), or an empty placeholder
    chains                 // channel: full-length transcripts for --chains (with the empty placeholder), or the placeholder alone
    proteindb             // path: preprocessed protein DB (one FASTA)
    miniprot_scored        // path: scored miniprot alignments (miniprot_parsed.gff)
    prot_train_gff         // path: miniprot_trainingGenes.gff
    prot_hc_hints          // path: miniprothint hc.gff
    tib_gff                // channel: Tiberius ab initio GFF (produced in main.nf for PROTEIN_EVIDENCE)
    tib_proteins           // channel: Tiberius protein FASTA (same, reused as the gffread output)
    params_map
    transcripts            // value: true in the evidence modes (StringTie assembly, ORF callers), false in mode proteins

    main:
    // ---- gene finders
    def vip = GENEFINDER_VIPSANIA(CH_GENOME, params_map, 'vipsania', false)

    // ---- the protein database for the DIAMOND runs of every gene set
    def db = DIAMOND_MAKEDB(proteindb)

    // ---- the ORF callers on the one StringTie assembly (evidence modes), or the
    // empty placeholder for their gene sets and DIAMOND tables (mode proteins)
    def d = params_map.drusilla ?: [:]
    def dru_gtf
    def td2_gff3
    def dia_dru
    def dia_td2
    if( transcripts ) {
        // ---- the one StringTie assembly as GTF and GFF3 (merged only if there are several)
        def merged_asm = STRINGTIE_MERGE(asm_gtf.collect())

        // ---- Drusilla ORFs on the filtered StringTie assembly
        def dru
        if( d.result ) {
            if( !file(d.result).exists() )
                error "params.drusilla.result is set, but the file does not exist: ${d.result}"
            dru = [gtf: channel.fromPath(d.result, checkIfExists: true).toList()]
        } else {
            if( truthy(d.unstranded_both_strands) && truthy(d.drop_unstranded) )
                log.warn "drusilla.unstranded_both_strands has no effect while drusilla.drop_unstranded is true: " +
                         "FILTER_STRINGTIE drops every transcript with strand '.' before Drusilla runs."
            def filt = FILTER_STRINGTIE(merged_asm.gtf)
            def drusillaModel = channel.value([])
            def useModel = false
            if( !d.weights && !d.cache_dir ) {
                def manifest = "${projectDir}/drusilla/model_cfg/${drusillaSetting(params_map, 'model')}.yaml"
                def url = tiberiusModelValue(manifest, 'weights_url')
                if( url ) {
                    drusillaModel = DOWNLOAD_DRUSILLA_MODEL(url, tiberiusModelValue(manifest, 'weights_sha256') ?: '').model
                    useModel = true
                }
            }
            def dru_raw = DRUSILLA_ANNOTATE(filt.gtf, CH_GENOME, useModel, drusillaModel)
            dru        = COLLAPSE_ORFS(dru_raw.gtf)
        }

        // ---- TD2 with the partial-ORF rescue
        def td        = TD_ALL(merged_asm.gtf, CH_GENOME)
        def pep_short = SHORTEN_INCOMPLETE_ORFS(td.pep)
        def dia_norm  = DIAMOND_NORM(td.pep, db.db)
        def dia_short = DIAMOND_SHORT(pep_short.pep_short, db.db)
        def rev       = CDS_CLASSIFY_AND_REVISE(dia_norm.tsv, dia_short.tsv, td.pep, pep_short.pep_short)
        def td2_gff   = TD_REVISED_TO_GENOME(rev.revised_pep, td.cdna, merged_asm.gtf, merged_asm.gff3)

        // ---- proteins of the ORF-caller gene sets and their DIAMOND runs (14-column)
        def pep_td2 = PROT_TD2('td2',      td2_gff.gff3, CH_GENOME)
        def pep_dru = PROT_DRU('drusilla', dru.gtf,     CH_GENOME)
        dru_gtf  = dru.gtf
        td2_gff3 = td2_gff.gff3
        dia_td2  = DIAMOND_TD2(pep_td2, db.db).tsv
        dia_dru  = DIAMOND_DRU(pep_dru, db.db).tsv
    } else {
        if( d.result )
            log.warn "drusilla.result is not used without transcript evidence (mode proteins): no ORF caller runs."
        def none = NO_TRANSCRIPTS()
        dru_gtf  = none
        td2_gff3 = none
        dia_td2  = none
        dia_dru  = none
    }

    // ---- proteins of the Vipsania gene set; the DIAMOND runs of both gene finders (14-column)
    def pep_vip = PROT_VIP('vipsania', vip.gff, CH_GENOME)
    def dia_tib = DIAMOND_TIB(tib_proteins, db.db)
    def dia_vip = DIAMOND_VIP(pep_vip,      db.db)

    // ---- RNA-Seq coverage; without a BAM the placeholder goes on, and
    // COMBINE_GENESETS runs combine_genesets.py without --coverage
    def cov_bg = have_bam ? BAM2COV(merged_bam).bg : merged_bam

    // ---- integrate
    def clade  = hcClade(params_map)
    def insect = clade.name?.toLowerCase() in ['insecta', 'lepidoptera']
    def min_locus = (params_map.combine?.min_locus != null) ? params_map.combine.min_locus : 4
    // learned selection: the model of params.combine.classifier (evidence modes)
    // or params.combine.classifier_proteins (mode proteins); false / null:
    // tier rules, the empty placeholder goes in. The clade is known to the
    // model when its HC-table entry names a classifier_clade; otherwise the
    // model's conservative locus threshold applies (--unknown-clade).
    def clfParam = transcripts ? params_map.combine?.classifier : params_map.combine?.classifier_proteins
    def classifier = clfParam ? channel.value(file(clfParam.toString(), checkIfExists: true)) : NO_CLASSIFIER()
    def unknown_clade = !(clade.entry?.classifier_clade)
    if( clfParam && unknown_clade )
        log.info "combine_genesets: clade ${clade.name ?: 'unknown'} has no classifier_clade in the HC table; " +
                 "the classifier uses its conservative locus threshold (--unknown-clade)."
    // the thresholds of params.combine.t_locus / t_iso override the model's;
    // with -resume a change here reruns only the combine step and what follows
    def clf_opts = [unknown_clade ? '--unknown-clade' : '',
                    params_map.combine?.t_locus != null ? "--t-locus ${params_map.combine.t_locus}" : '',
                    params_map.combine?.t_iso   != null ? "--t-iso ${params_map.combine.t_iso}"     : '']
                   .findAll().join(' ')

    def combined = COMBINE_GENESETS(
        tib_gff, vip.gff, dru_gtf, td2_gff3,
        intron_hints, cov_bg, chains,
        dia_tib.tsv, dia_vip.tsv, dia_dru, dia_td2,
        prot_train_gff, miniprot_scored, prot_hc_hints,
        classifier, insect, min_locus, clf_opts)

    emit:
    gff      = combined.gff3
    report   = combined.report
    chunks   = vip.chunks            // for cmscan in POSTPROCESS
    versions = vip.versions
}
