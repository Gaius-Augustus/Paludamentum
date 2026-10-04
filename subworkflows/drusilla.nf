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

include { FILTER_STRINGTIE; DOWNLOAD_DRUSILLA_MODEL; DRUSILLA_ANNOTATE; FIX_ORFS;
          GENEFINDER_LGB_FILTER; HINT_RESCUE_LOCI; HINT_RESCUE_TIBERIUS } from '../modules/drusilla.nf'
include { DOWNLOAD_TIBERIUS_WEIGHTS as DOWNLOAD_RESCUE_WEIGHTS } from '../modules/genefinder.nf'
include { rescueEnabled; rescueModel; tiberiusModelValue; drusillaSetting } from '../lib_nf/functions.nf'

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
    // The released Drusilla model is downloaded once, with the weights_url and
    // weights_sha256 of its manifest in model_cfg/ of the Drusilla submodule
    // (the Drusilla of the image). Not for local weights, nor with a cache_dir,
    // nor for a name that is not there: the Drusilla registry resolves those.
    def drusillaModel = channel.value([])
    def useModel = false
    if( !d.weights && !d.cache_dir ) {
        def manifest = "${projectDir}/drusilla/model_cfg/${drusillaSetting(params, 'model')}.yaml"
        def url = tiberiusModelValue(manifest, 'weights_url')
        if( url ) {
            drusillaModel = DOWNLOAD_DRUSILLA_MODEL(url, tiberiusModelValue(manifest, 'weights_sha256') ?: '').model
            useModel = true
        }
    }
    raw      = DRUSILLA_ANNOTATE(filtered.gtf, CH_GENOME, useModel, drusillaModel)
    orfs     = FIX_ORFS(raw.gtf, raw.partial, raw.partial5, miniprot_gff, hc_hints, CH_GENOME)
    kept     = GENEFINDER_LGB_FILTER(
        genefinder, genefinder_gff, miniprot_gff, hc_hints, CH_GENOME,
        file(drusillaSetting(params, 'lgb_model'), checkIfExists: true)
    )

    def genefinder_out = kept.gtf
    if( rescueEnabled(params) ) {
        // Tiberius model of the rescue (rescueModel in lib_nf/functions.nf): a
        // model configuration file, or the name of a model in model_cfg/ of the
        // rescue image. Its weights are staged: those of params.tiberius.model_dir
        // if it is the model of the run, else the archive of its weights_url,
        // downloaded once. The weights_url of a name is read from model_cfg/ of
        // the Tiberius submodule, which names the same weights as the image; a
        // name that is not there leaves the download to Tiberius in the task.
        def rm = rescueModel(params)
        def cfgFile = rm.file ? file(rm.file, checkIfExists: true) : []
        def weights = channel.value([])
        def useWeights = true
        if( rm.weights ) {
            if( !file(params.tiberius.model_dir).isDirectory() ) error "params.tiberius.model_dir is not a directory: ${params.tiberius.model_dir}"
            weights = channel.fromPath("${params.tiberius.model_dir}/*", type: 'any').collect()
        } else {
            def url = tiberiusModelValue(rm.file ?: "${projectDir}/tiberius/model_cfg/${rm.name}.yaml", 'weights_url')
            if( url ) weights = DOWNLOAD_RESCUE_WEIGHTS(url).collect()
            else      useWeights = false
        }
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
        rescued = HINT_RESCUE_TIBERIUS(todo, rm.name ?: '', cfgFile, useWeights, weights)
        genefinder_out = kept.gtf.mix(rescued.gtf).collect()
    }

    emit:
    orfs       = orfs.gtf
    genefinder = genefinder_out
}
