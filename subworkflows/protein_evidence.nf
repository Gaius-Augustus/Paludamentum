nextflow.enable.dsl=2

include { MINIPROT_ALIGN; MINIPROTHINT_CONVERT; ALN2HINTS; PREPROCESS_PROTEINDB } from '../modules/proteins.nf'

// Protein evidence pipeline. Caller passes the ab initio Tiberius protein
// FASTA (or an empty file placeholder); PREPROCESS_PROTEINDB uses it to
// select source species when the protein database has more than one million
// OrthoDB-style entries, else it is ignored. The gene finders themselves run
// in HC_COMBINED (not here) so Tiberius is never scheduled twice.

workflow PROTEIN_EVIDENCE {

    take:
    CH_GENOME
    CH_PROTEINS
    CH_SCORE
    tib_proteins_ch        // channel: Tiberius protein FASTA (empty file is fine)
    params_map

    main:
    if( !params_map.proteins && !params_map.odb12Partitions ) {
        error "params.proteins or params.odb12Partitions is required for protein evidence modes"
    }
    if( !params_map.scoring_matrix )  error "params.scoring_matrix is required for protein evidence modes"

    def proteindb_ch = PREPROCESS_PROTEINDB(CH_PROTEINS, tib_proteins_ch)

    def scored_ch     = MINIPROT_ALIGN(CH_GENOME, proteindb_ch, CH_SCORE)
    def prot_gtf_ch   = MINIPROTHINT_CONVERT(scored_ch.gff)
    def prot_hints_ch = ALN2HINTS(prot_gtf_ch.gtf)

    emit:
    proteindb     = proteindb_ch
    scored_gff    = scored_ch.gff
    prot_gtf      = prot_gtf_ch.gtf
    prot_traingff = prot_gtf_ch.traingff
    prot_hc_hints = prot_gtf_ch.hc_hints
    prot_hints    = prot_hints_ch.hints
}
