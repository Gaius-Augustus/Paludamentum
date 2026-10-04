nextflow.enable.dsl=2

include { MINIMAP2_MAP } from '../modules/isoseq.nf'
include { FILTER_ALIGNMENT as FILTER_ISOSEQ; EMPTY_FILE } from '../modules/util.nf'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_ISO; BAM2HINTS as BAM2HINTS_ISO } from '../modules/rnaseq.nf'
include { DOWNLOAD_SRA_ISOSEQ } from '../modules/download.nf'
include { STRINGTIE_ASSEMBLE_ISO } from '../modules/assembly.nf'
include { VARUS_INPUT as VARUS_INPUT_ISO; MERGE_INTRON_HINTS as MERGE_INTRON_HINTS_ISO } from '../modules/varus.nf'
include { asList; varusInputs } from '../lib_nf/functions.nf'

workflow ISOSEQ_EVIDENCE {

    take:
    CH_GENOME
    params_map
    asm_mode     // 'per_sample': one StringTie assembly per BAM; 'merged': one of the merged BAM; 'none'

    main:
    empty_file = EMPTY_FILE()

    def DO_ISO_LOCAL = params_map.isoseq && params_map.isoseq.size() > 0
    def DO_ISO = DO_ISO_LOCAL || (params_map.isoseq_sra && params_map.isoseq_sra.size() > 0)

    // pyVARUS directories (--longreads) bring their StringTie assembly and intron hints
    def varus = varusInputs(params_map.isoseq_varus, 'isoseq_varus', 'longreads')
    if( asm_mode == 'merged' && varus && (varus.size() > 1 || DO_ISO) )
        error "isoseq_varus: the Drusilla flow needs one StringTie assembly of all Iso-Seq reads, and a pyVARUS " +
              "directory holds the assembly of its own reads only. Pass a single pyVARUS directory and no other " +
              "Iso-Seq reads, or merge the BAMs of the runs (`varus replay` rebuilds a dropped one) and pass the " +
              "output of `varus assemble GENOME --long merged.bam --outdir DIR` as the one directory."

    def hints_out    = empty_file
    def asm_gtf_out  = channel.empty()
    def bam_out      = channel.empty()

    if( DO_ISO ) {
        CH_ISO_LOCAL = DO_ISO_LOCAL ? channel.fromPath(params_map.isoseq, checkIfExists:true) : channel.empty()

        CH_ISO = CH_ISO_LOCAL
        if( params_map.isoseq_sra ) {
            CH_ISO_SRA_IDS = channel.fromList(asList(params_map.isoseq_sra))
            CH_RNASEQ_ISO_SRA = DOWNLOAD_SRA_ISOSEQ(CH_ISO_SRA_IDS)
            CH_ISO = CH_ISO_LOCAL.mix(CH_RNASEQ_ISO_SRA.map { _acc, f -> f })
        }

        iso_bam  = MINIMAP2_MAP(CH_GENOME, CH_ISO)
        FILTER_ISOSEQ(iso_bam.bam)

        iso_bams = channel.empty().mix(FILTER_ISOSEQ.out)

        iso_merged = SAMTOOLS_MERGE_ISO(iso_bams.collect())
        iso_hints  = BAM2HINTS_ISO(iso_merged.bam, CH_GENOME)

        if( asm_mode == 'per_sample' )  asm_gtf_out = STRINGTIE_ASSEMBLE_ISO(iso_bams).gtf
        else if( asm_mode == 'merged' ) asm_gtf_out = STRINGTIE_ASSEMBLE_ISO(iso_merged.bam).gtf

        hints_out    = iso_hints.hints
        bam_out      = iso_merged.bam
    }

    if( varus ) {
        vi = VARUS_INPUT_ISO(channel.fromList(varus), channel.value(file(params_map.genome)))
        if( asm_mode == 'per_sample' )  asm_gtf_out = asm_gtf_out.mix(vi.gtf)
        else if( asm_mode == 'merged' ) asm_gtf_out = vi.gtf
        // One hint source is used as it is; several are summed as bam2hints
        // on all reads would count them.
        def sources = DO_ISO ? hints_out.mix(vi.hints) : vi.hints
        hints_out = (varus.size() == 1 && !DO_ISO) ? vi.hints \
            : MERGE_INTRON_HINTS_ISO(sources.collect(sort: { f -> f.name })).hints
    }

    emit:
    hints    = hints_out
    asm_gtf  = asm_gtf_out
    bam      = bam_out
}
