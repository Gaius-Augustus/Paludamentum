nextflow.enable.dsl=2

include { HISAT2_BUILD; HISAT2_MAP_SINGLE; HISAT2_MAP_PAIRED;
          SAMTOOLS_MERGE as SAMTOOLS_MERGE_RNA;
          BAM2HINTS as BAM2HINTS_RNA } from '../modules/rnaseq.nf'

include { FILTER_ALIGNMENT as FILTER_PE; FILTER_ALIGNMENT as FILTER_SE } from '../modules/util.nf'
include { DOWNLOAD_SRA_PAIRED; DOWNLOAD_SRA_SINGLE } from '../modules/download.nf'
include { STRINGTIE_ASSEMBLE_RNA } from '../modules/assembly.nf'
include { EMPTY_FILE } from '../modules/util.nf'
include { asList } from '../lib_nf/functions.nf'

// Library name of a read file: the file name without the FASTQ suffixes and
// without a trailing read-pair marker (_1, _R1, .1). Two libraries with the same
// file name in different directories still collide; see SAMTOOLS_MERGE.
def libraryId(r1) {
    return r1.name
        .replaceFirst(/\.(fastq|fq|fasta|fa)(\.gz)?$/, '')
        .replaceFirst(/([._-]R?1)$/, '')
}

workflow RNASEQ_EVIDENCE {

    take:
    CH_GENOME
    params_map
    asm_mode     // 'per_sample': one StringTie assembly per BAM; 'merged': one of the merged BAM; 'none'

    main:
    empty_file = EMPTY_FILE()

    def DO_SE_LOCAL = params_map.rnaseq_single && params_map.rnaseq_single.size() > 0
    def DO_PE_LOCAL = params_map.rnaseq_paired && params_map.rnaseq_paired.size() > 0
    def DO_SE = DO_SE_LOCAL || (params_map.rnaseq_sra_single && params_map.rnaseq_sra_single.size() > 0)
    def DO_PE = DO_PE_LOCAL || (params_map.rnaseq_sra_paired && params_map.rnaseq_sra_paired.size() > 0)

    def DO_BAM = params_map.rnaseq_bam && params_map.rnaseq_bam.size() > 0

    def hints_out    = empty_file
    def asm_gtf_out  = channel.empty()
    def bam_out      = channel.empty()

    if( DO_SE || DO_PE || DO_BAM ) {

    // Build channels (local + SRA)
    CH_PAIRED_LOCAL = channel.empty()
    if( DO_PE ) {
        def pe = params_map.rnaseq_paired

        /*
        * Case 1: glob pattern or string → use fromFilePairs (NON-flat)
        *   "/path/*_{1,2}.fastq.gz"
        */
        if( pe instanceof CharSequence ) {

            CH_PAIRED_LOCAL =
                channel.fromFilePairs(pe, size: 2, checkIfExists: true)
                // emits: tuple(id, [r1, r2])
        }

        /*
        * Case 2: list of explicit PAIRS
        *   - [r1, r2]
        *   - [r1, r2]
        */
        else if(
            pe instanceof List &&
            pe.every { p -> p instanceof List && p.size() == 2 }
        ) {

            CH_PAIRED_LOCAL =
                channel.fromList(pe)
                .map { pair ->
                    def r1 = file(pair[0])
                    def r2 = file(pair[1])
                    tuple(libraryId(r1), [r1, r2])
                }
        }
        /*
        * Case 3: flat list of exactly two files → treat as ONE library
        *   - r1
        *   - r2
        */
        else if(
            pe instanceof List &&
            pe.size() == 2 &&
            pe.every { p -> p instanceof CharSequence }
        ) {

            def r1 = file(pe[0])
            def r2 = file(pe[1])
            CH_PAIRED_LOCAL =
                channel.of( tuple(libraryId(r1), [r1, r2]) )
        }
        else if (DO_PE_LOCAL) {
            error """Invalid rnaseq_paired format.
Expected one of:
  (a) a single glob string, e.g.
        rnaseq_paired: "/abs/path/RNA/*_{1,2}.fastq.gz"
  (b) a YAML list of [r1, r2] pairs, e.g.
        rnaseq_paired:
          - ["/abs/path/sample1_1.fastq.gz", "/abs/path/sample1_2.fastq.gz"]
          - ["/abs/path/sample2_1.fastq.gz", "/abs/path/sample2_2.fastq.gz"]
  (c) a flat list of exactly two FASTQ files (one library).
Got: ${pe?.getClass()?.simpleName} -> ${pe}"""
        }

        if( params_map.rnaseq_sra_paired ) {
            CH_RNASEQ_SRA_IDS_PAIRED = channel.fromList(asList(params_map.rnaseq_sra_paired))
            CH_RNASEQ_PAIRED_SRA = DOWNLOAD_SRA_PAIRED(CH_RNASEQ_SRA_IDS_PAIRED)
            CH_PAIRED_SRA = CH_RNASEQ_PAIRED_SRA.map { acc, r1, r2 -> tuple(acc, [r1, r2]) }
            CH_PAIRED = CH_PAIRED_LOCAL.mix(CH_PAIRED_SRA)
        } else CH_PAIRED = CH_PAIRED_LOCAL
    }

    CH_SINGLE = channel.empty()
    if( DO_SE ) {
        CH_SINGLE_LOCAL = DO_SE_LOCAL ? channel.fromPath(params_map.rnaseq_single, checkIfExists:true) : channel.empty()
        if( params_map.rnaseq_sra_single ) {
            CH_RNASEQ_SRA_IDS_SINGLE = channel.fromList(asList(params_map.rnaseq_sra_single))
            CH_RNASEQ_SINGLE_SRA = DOWNLOAD_SRA_SINGLE(CH_RNASEQ_SRA_IDS_SINGLE)
            CH_SINGLE = CH_SINGLE_LOCAL.mix(CH_RNASEQ_SINGLE_SRA.map { _acc, f -> f })
        } else CH_SINGLE = CH_SINGLE_LOCAL
    }

    if ( DO_SE || DO_PE) {
        index = HISAT2_BUILD(CH_GENOME)
    }

    rnaseq_bams = channel.empty()

    if( DO_SE ) {
        map_se = HISAT2_MAP_SINGLE(index.idxdir, CH_SINGLE)
        FILTER_SE(map_se.bam)
        rnaseq_bams = rnaseq_bams.mix(FILTER_SE.out)
    }

    if( DO_PE ) {
        map_pe = HISAT2_MAP_PAIRED(index.idxdir, CH_PAIRED)
        FILTER_PE(map_pe.bam)
        rnaseq_bams = rnaseq_bams.mix(FILTER_PE.out)
    }

    if( DO_BAM ) {
        rnaseq_bams = rnaseq_bams.mix(
            channel.fromList(asList(params_map.rnaseq_bam)).map { f -> file(f) }
        )
    }

    rnaseq_merged = SAMTOOLS_MERGE_RNA(rnaseq_bams.collect())
    rnaseq_hints  = BAM2HINTS_RNA(rnaseq_merged.bam, CH_GENOME)

    if( asm_mode == 'per_sample' )  asm_gtf_out = STRINGTIE_ASSEMBLE_RNA(rnaseq_bams).gtf
    else if( asm_mode == 'merged' ) asm_gtf_out = STRINGTIE_ASSEMBLE_RNA(rnaseq_merged.bam).gtf

    hints_out    = rnaseq_hints.hints
    bam_out      = rnaseq_merged.bam
    }

    emit:
    hints    = hints_out
    asm_gtf  = asm_gtf_out
    bam      = bam_out
}
