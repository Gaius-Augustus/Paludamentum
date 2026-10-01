nextflow.enable.dsl=2

include { STRINGTIE_MERGE } from '../modules/assembly.nf'
include { TD_ALL; SHORTEN_INCOMPLETE_ORFS; CDS_CLASSIFY_AND_REVISE } from '../modules/transdecoder.nf'
include { DIAMOND_MAKEDB; DIAMOND_BLASTP as DIAMOND_BLASTP_NORM; DIAMOND_BLASTP as DIAMOND_BLASTP_SHORT; DIAMOND_BLASTP as DIAMOND_BLASTP_REV } from '../modules/diamond.nf'
include { HC_SUPPORTED; HC_FORMAT_FILTER } from '../modules/hc.nf'

workflow HC_GENES {

    take:
    asm_gtf_ch
    CH_GENOME
    proteindb

    main:
    // Without any assembly (every library dropped for its alignment rate) the
    // run would end without HC genes and without the final annotation.
    def assemblies = asm_gtf_ch.collect().ifEmpty {
        error "No transcript assembly: every RNA-Seq/Iso-Seq library was dropped (alignment rate below params.min_alignment_rate = ${params.min_alignment_rate}). Check that the reads belong to this genome, or lower params.min_alignment_rate."
    }
    asm      = STRINGTIE_MERGE(assemblies)
    td_all   = TD_ALL(asm.gtf, CH_GENOME)

    pep_short = SHORTEN_INCOMPLETE_ORFS(td_all.pep)
    db        = DIAMOND_MAKEDB(proteindb)

    dia_norm  = DIAMOND_BLASTP_NORM(td_all.pep, db.db)
    dia_short = DIAMOND_BLASTP_SHORT(pep_short.pep_short, db.db)

    rev       = CDS_CLASSIFY_AND_REVISE(dia_norm.tsv, dia_short.tsv, td_all.pep, pep_short.pep_short)
    dia_rev   = DIAMOND_BLASTP_REV(rev.revised_pep, db.db)

    hc = HC_SUPPORTED(
        dia_rev.tsv,
        rev.revised_pep,
        proteindb,
        asm.gtf,
        asm.gff3,    
        td_all.cdna
    )

    HC_FORMAT_FILTER(hc.training_gff, CH_GENOME)

    emit:
    HC_FORMAT_FILTER.out
}
