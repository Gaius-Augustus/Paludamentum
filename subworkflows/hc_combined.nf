nextflow.enable.dsl=2

// Unified HC-gene step for evidence modes (rnaseq / isoseq / mixed).
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
//     For Insecta / Lepidoptera clades, --drop-lone-single-exon tiberius
//     is added.

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
include { hcClade; tiberiusModelValue; drusillaSetting; truthy } from '../lib_nf/functions.nf'

workflow HC_COMBINED {

    take:
    CH_GENOME
    asm_gtf                // channel: one or more StringTie assembly GTFs
    merged_bam             // channel: one coordinate-sorted BAM of all RNA / Iso-Seq reads, or an empty placeholder
    have_bam               // value: false without a BAM (pyVARUS directories, a given assembly without reads)
    intron_hints           // channel: the intron hints of BAM2HINTS (two files in mixed mode), or an empty placeholder
    proteindb              // path: preprocessed protein DB (one FASTA)
    miniprot_scored        // path: scored miniprot alignments (miniprot_parsed.gff)
    prot_train_gff         // path: miniprot_trainingGenes.gff
    prot_hc_hints          // path: miniprothint hc.gff
    tib_gff                // channel: Tiberius ab initio GFF (produced in main.nf for PROTEIN_EVIDENCE)
    tib_proteins           // channel: Tiberius protein FASTA (same, reused as the gffread output)
    params_map

    main:
    // ---- gene finders
    def vip = GENEFINDER_VIPSANIA(CH_GENOME, params_map, 'vipsania', false)

    // ---- the one StringTie assembly as GTF and GFF3 (merged only if there are several)
    def merged_asm = STRINGTIE_MERGE(asm_gtf.collect())

    // ---- Drusilla ORFs on the filtered StringTie assembly
    def d = params_map.drusilla ?: [:]
    def dru
    if( d.result ) {
        if( !file(d.result).exists() )
            error "params.drusilla.result is set, but the file does not exist: ${d.result}"
        dru = [gtf: channel.fromPath(d.result, checkIfExists: true).toList()]
    } else {
        if( truthy(d.unstranded_both_strands) && (d.drop_unstranded == null || truthy(d.drop_unstranded)) )
            log.warn "drusilla.unstranded_both_strands has no effect while drusilla.drop_unstranded is true (the default): " +
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
    def db        = DIAMOND_MAKEDB(proteindb)
    def pep_short = SHORTEN_INCOMPLETE_ORFS(td.pep)
    def dia_norm  = DIAMOND_NORM(td.pep, db.db)
    def dia_short = DIAMOND_SHORT(pep_short.pep_short, db.db)
    def rev       = CDS_CLASSIFY_AND_REVISE(dia_norm.tsv, dia_short.tsv, td.pep, pep_short.pep_short)
    def td2_gff   = TD_REVISED_TO_GENOME(rev.revised_pep, td.cdna, merged_asm.gtf, merged_asm.gff3)

    // ---- proteins for each gene set
    def pep_vip = PROT_VIP('vipsania', vip.gff,     CH_GENOME)
    def pep_td2 = PROT_TD2('td2',      td2_gff.gff3, CH_GENOME)
    def pep_dru = PROT_DRU('drusilla', dru.gtf,     CH_GENOME)

    // ---- four combine_genesets DIAMOND runs (14-column)
    def dia_tib = DIAMOND_TIB(tib_proteins, db.db)
    def dia_vip = DIAMOND_VIP(pep_vip,      db.db)
    def dia_td2 = DIAMOND_TD2(pep_td2,      db.db)
    def dia_dru = DIAMOND_DRU(pep_dru,      db.db)

    // ---- RNA-Seq coverage; without a BAM the placeholder goes on, and
    // COMBINE_GENESETS runs combine_genesets.py without --coverage
    def cov_bg = have_bam ? BAM2COV(merged_bam).bg : merged_bam

    // ---- integrate
    def clade  = hcClade(params_map)
    def insect = clade.name?.toLowerCase() in ['insecta', 'lepidoptera']
    def min_locus = (params_map.combine?.min_locus != null) ? params_map.combine.min_locus : 4

    def combined = COMBINE_GENESETS(
        tib_gff, vip.gff, dru.gtf, td2_gff.gff3,
        intron_hints, cov_bg,
        dia_tib.tsv, dia_vip.tsv, dia_dru.tsv, dia_td2.tsv,
        prot_train_gff, miniprot_scored, prot_hc_hints,
        insect, min_locus)

    emit:
    gff      = combined.gff3
    report   = combined.report
    chunks   = vip.chunks            // for cmscan in POSTPROCESS
    versions = vip.versions
}
