nextflow.enable.dsl=2

// Drusilla flow of evidence runs with transcripts and a vertebrate gene finder
// model (see hcMethod in lib_nf/functions.nf):
//   1. the StringTie assembly is filtered by length, coverage and TPM,
//   2. Drusilla predicts the ORFs of the kept transcripts,
//   3. stop and start codons of the ORFs are fixed with protein alignments
//      (fix_stop, fix_start), subsequence isoforms are collapsed,
//   4. a LightGBM model filters the gene finder predictions with protein
//      alignment (miniprot) and hint (miniprothint) features,
//   5. partial gene finder genes are predicted again by Tiberius with the hints
//      of their best protein chain, ab initio where there is no chain (hint
//      rescue; the ORF-agreement filter rescue_orf_filter is off by default).
// The ORFs, the kept and the rescued gene finder genes are merged by
// MERGE_GENEFINDER_TRAIN in main.nf, as in the TransDecoder flow.

include { FILTER_STRINGTIE; DRUSILLA_ANNOTATE; FIX_ORFS; GENEFINDER_LGB_FILTER;
          HINT_RESCUE_LOCI; HINT_RESCUE_TIBERIUS } from '../modules/drusilla.nf'
include { truthy } from '../lib_nf/functions.nf'

workflow DRUSILLA_HC {

    take:
    asm_gtf          // one StringTie assembly (GTF with cov and TPM)
    CH_GENOME
    genefinder       // name of the gene finder, names the published files
    genefinder_gff   // ab initio predictions of the gene finder
    miniprot_gff     // scored miniprot alignments
    hc_hints         // miniprothint hc.gff

    main:
    def d = params.drusilla
    filtered = FILTER_STRINGTIE(asm_gtf)
    raw      = DRUSILLA_ANNOTATE(filtered.gtf, CH_GENOME)
    orfs     = FIX_ORFS(raw.gtf, raw.partial, raw.partial5, miniprot_gff, hc_hints, CH_GENOME)
    kept     = GENEFINDER_LGB_FILTER(
        genefinder, genefinder_gff, miniprot_gff, hc_hints, CH_GENOME,
        file(d.lgb_model, checkIfExists: true)
    )

    def genefinder_out = kept.gtf
    if( d.rescue == null || truthy(d.rescue) ) {
        // Tiberius model of the rescue: the model of the run, else vertebrates
        def cfg = d.rescue_model_cfg ?: (genefinder == 'tiberius' ? params.tiberius.model_cfg : 'vertebrates')
        def cfgName = new File(cfg.toString()).name.replaceFirst(/\.ya?ml$/, '')
        loci    = HINT_RESCUE_LOCI(kept.partial, kept.gtf, orfs.gtf, miniprot_gff, hc_hints, CH_GENOME)
        // The GPU task runs only with a --hints Tiberius and at least one locus
        def tiberius = d.rescue_tiberius ?: 'the tiberius.py of the image'
        todo = loci.loci
            .filter { fasta, _hints, _manifest, hintsOk ->
                if( hintsOk != 'true' ) {
                    log.warn "Hint rescue skipped: ${tiberius} has no --hints option. " +
                             "The rescue needs the Tiberius branch hint_integration: use the image of the label " +
                             "hint_rescue in conf/base.config, or set drusilla.rescue_tiberius. " +
                             "Set drusilla.rescue = false to silence this warning."
                    return false
                }
                if( fasta.size() == 0 ) {
                    log.info "Hint rescue: no loci"
                    return false
                }
                return true
            }
            .map { fasta, hints, manifest, _hintsOk -> tuple(fasta, hints, manifest) }
        rescued = HINT_RESCUE_TIBERIUS(todo, cfgName)
        genefinder_out = kept.gtf.mix(rescued.gtf).collect()
    }

    emit:
    orfs       = orfs.gtf
    genefinder = genefinder_out
}
