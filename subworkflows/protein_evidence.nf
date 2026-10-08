nextflow.enable.dsl=2

include { MINIPROT_ALIGN; MINIPROTHINT_CONVERT; ALN2HINTS; PREPROCESS_PROTEINDB } from '../modules/proteins.nf'
include { PROTEIN_FROM_GFF } from '../modules/genefinder.nf'
include { genefinderEnabled; resolveGenefinder } from '../lib_nf/functions.nf'
include { GENEFINDER } from './genefinder.nf'

workflow PROTEIN_EVIDENCE {

    take:
    CH_GENOME
    CH_PROTEINS
    CH_SCORE
    params_map

    main:
    def proteindb_ch
    def genefinder_gff_ch = channel.empty()
    def scored_ch
    def prot_gtf_ch
    def prot_hints_ch

    if( !params_map.proteins && !params_map.odb12Partitions ) {
        error "params.proteins or params.odb12Partitions is required for protein evidence modes"
    }
    if( !params_map.scoring_matrix )  error "params.scoring_matrix is required for protein evidence modes"

    if( genefinderEnabled(params_map) ) {
        genefinder_gff_ch = GENEFINDER(CH_GENOME, params_map, false)

        // the ab initio proteins help to select source species in large protein databases
        def ab_initio_prot = PROTEIN_FROM_GFF(resolveGenefinder(params_map), genefinder_gff_ch, CH_GENOME)
        proteindb_ch = PREPROCESS_PROTEINDB(CH_PROTEINS, ab_initio_prot)
    }
    else {
        proteindb_ch = CH_PROTEINS
    }

    scored_ch   = MINIPROT_ALIGN(CH_GENOME, proteindb_ch, CH_SCORE)
    prot_gtf_ch = MINIPROTHINT_CONVERT(scored_ch.gff)
    prot_hints_ch = ALN2HINTS(prot_gtf_ch.gtf)

    emit:
    proteindb     = proteindb_ch
    genefinder_gff = genefinder_gff_ch
    scored_gff    = scored_ch.gff
    prot_gtf      = prot_gtf_ch.gtf
    prot_traingff = prot_gtf_ch.traingff
    prot_hc_hints = prot_gtf_ch.hc_hints
    prot_hints    = prot_hints_ch.hints
}
