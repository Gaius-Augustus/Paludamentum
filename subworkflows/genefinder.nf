nextflow.enable.dsl=2

// Runs the selected gene finder on the genome and emits the merged ab initio
// annotation. The gene finder is chosen by resolveGenefinder(): 'tiberius' or
// 'vipsania'. Each has its own params block (params.tiberius, params.vipsania).
// The genome chunks are emitted too: cmscan of the ncRNA annotation
// (subworkflows/postprocess.nf) runs on them instead of a second split.

include { truthy; resolveGenefinder; tiberiusModelValue } from '../lib_nf/functions.nf'
include { RUN_TIBERIUS; DOWNLOAD_TIBERIUS_WEIGHTS; SPLIT_GENOME; MERGE_GENEFINDER;
          MERGE_GENEFINDER_EVI } from '../modules/genefinder.nf'
include { DOWNLOAD_VIPSANIA_MODEL; RUN_VIPSANIA } from '../modules/vipsania.nf'

workflow GENEFINDER {

    take:
    CH_GENOME
    params_map
    publish_top      // true: publish the ab initio annotation in outdir, false: in outdir/intermediate

    main:
    def tool = resolveGenefinder(params_map)
    def cfg  = params_map[tool] ?: [:]
    def predictions
    if( cfg.result && !file(cfg.result).exists() )
        error "params.${tool}.result is set, but the file does not exist: ${cfg.result}"
    def useResult = cfg.result as boolean
    def finetune  = tool == 'vipsania' && truthy(cfg.finetune)

    // One split of the genome, for the gene finder tasks and for cmscan
    // (ncrna.run); a reused result and Vipsania finetuning read no chunks.
    def genome_chunks = channel.empty()
    if( (!useResult && !finetune) || truthy(params_map.ncrna?.run) ) {
        genome_chunks = SPLIT_GENOME(
            CH_GENOME,
            cfg.min_split_size ?: 20000000,
            cfg.max_files ?: 20
        ).chunks
    }

    // Vipsania model files as a value channel: a local directory, or a download.
    def models = channel.empty()
    if( tool == 'vipsania' && !useResult ) {
        if( !cfg.model ) error "params.vipsania.model is required (clade name or model id)."
        if( cfg.model_dir ) {
            if( !file(cfg.model_dir).isDirectory() ) error "params.vipsania.model_dir is not a directory: ${cfg.model_dir}"
            models = channel.fromPath("${cfg.model_dir}/*", type: 'any').collect()
        } else {
            models = DOWNLOAD_VIPSANIA_MODEL(cfg.model).collect()
        }
    }

    if( useResult ) {
        // reuse an existing prediction instead of running the gene finder
        predictions = channel.fromPath(cfg.result, checkIfExists: true).toList()

    } else if( finetune ) {
        // Vipsania finetunes on the FASTA it annotates, so the genome is not split.
        predictions = RUN_VIPSANIA(CH_GENOME, models, cfg.model).gtf.toList()

    } else {
        if( tool == 'vipsania' ) {
            predictions = RUN_VIPSANIA(genome_chunks.flatten(), models, cfg.model).gtf.toList()
        } else {
            if( !cfg.model_cfg ) error "params.tiberius.model_cfg is required."
            if( !file(cfg.model_cfg).exists() ) error "params.tiberius.model_cfg is not a file: ${cfg.model_cfg}"
            // The extracted weights are staged into every task: those of
            // params.tiberius.model_dir, else the archive of weights_url,
            // downloaded once. Only a model configuration without weights_url
            // leaves the download to Tiberius in each task.
            def weights = channel.value([])
            def useWeights = true
            if( cfg.model_dir ) {
                if( !file(cfg.model_dir).isDirectory() ) error "params.tiberius.model_dir is not a directory: ${cfg.model_dir}"
                weights = channel.fromPath("${cfg.model_dir}/*", type: 'any').collect()
            } else {
                def url = tiberiusModelValue(cfg.model_cfg, 'weights_url')
                if( url ) weights = DOWNLOAD_TIBERIUS_WEIGHTS(url).collect()
                else      useWeights = false
            }
            predictions = RUN_TIBERIUS(genome_chunks.flatten(), cfg.model_cfg, useWeights, weights).toList()
        }
    }

    // Published in outdir (ab initio mode) or in outdir/intermediate (evidence modes)
    if( publish_top ) MERGE_GENEFINDER(tool, predictions)
    else              MERGE_GENEFINDER_EVI(tool, predictions)

    emit:
    gff    = publish_top ? MERGE_GENEFINDER.out : MERGE_GENEFINDER_EVI.out
    chunks = genome_chunks   // the genome chunks (one list), or empty without a split
}
