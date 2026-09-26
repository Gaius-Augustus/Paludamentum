nextflow.enable.dsl=2

// Runs the selected gene finder on the genome and emits the merged ab initio
// annotation. The gene finder is chosen by resolveGenefinder(): 'tiberius' or
// 'vipsania'. Each has its own params block (params.tiberius, params.vipsania).

include { truthy; resolveGenefinder } from '../lib_nf/functions.nf'
include { RUN_TIBERIUS; SPLIT_GENOME; MERGE_GENEFINDER; MERGE_GENEFINDER_EVI } from '../modules/genefinder.nf'
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
    def useResult = cfg.result && file(cfg.result).exists()

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

    } else if( tool == 'vipsania' && truthy(cfg.finetune) ) {
        // Vipsania finetunes on the FASTA it annotates, so the genome is not split.
        predictions = RUN_VIPSANIA(CH_GENOME, models, cfg.model).gtf.toList()

    } else {
        def chunks = SPLIT_GENOME(
            CH_GENOME,
            cfg.min_split_size ?: 20000000,
            cfg.max_files ?: 20
        ).chunks.flatten()

        if( tool == 'vipsania' ) {
            predictions = RUN_VIPSANIA(chunks, models, cfg.model).gtf.toList()
        } else {
            if( !cfg.model_cfg ) error "params.tiberius.model_cfg is required."
            predictions = RUN_TIBERIUS(chunks, cfg.model_cfg).toList()
        }
    }

    // Published in outdir (ab initio mode) or in outdir/intermediate (evidence modes)
    if( publish_top ) MERGE_GENEFINDER(tool, predictions)
    else              MERGE_GENEFINDER_EVI(tool, predictions)

    emit:
    publish_top ? MERGE_GENEFINDER.out : MERGE_GENEFINDER_EVI.out
}
