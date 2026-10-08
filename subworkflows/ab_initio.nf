nextflow.enable.dsl=2

include { genefinderEnabled } from '../lib_nf/functions.nf'
include { GENEFINDER } from './genefinder.nf'

// Ab initio mode: run the gene finder only, parallelized over genome chunks.
workflow AB_INITIO {

    take:
    CH_GENOME
    params_map

    main:
    if( !genefinderEnabled(params_map) ) {
        error "The ab initio mode requires a gene finder: set params.tiberius.run=true or params.vipsania.run=true."
    }
    GENEFINDER(CH_GENOME, params_map, true)

    emit:
    gff    = GENEFINDER.out.gff      // <tool>_ab_initio.gff3, the input of POSTPROCESS
    chunks = GENEFINDER.out.chunks   // genome chunks, for cmscan in POSTPROCESS
    versions = GENEFINDER.out.versions   // version line of the gene finder, for POSTPROCESS
}
