nextflow.enable.dsl=2

include { CONCAT_HINTS; EMPTY_FILE } from './modules/util.nf'
include { MERGE_GENEFINDER_TRAIN; MERGE_GENEFINDER_TRAIN_PRIO; PROTEIN_FROM_GFF_FINAL } from './modules/genefinder.nf'
include { inferMode; resolveGenefinder; genefinderEnabled; asList; hcMethod } from './lib_nf/functions.nf'
include { HC_FORMAT_FILTER } from './modules/hc.nf'

include { INPUTS } from './subworkflows/inputs.nf'
include { PROTEIN_EVIDENCE } from './subworkflows/protein_evidence.nf'
include { RNASEQ_EVIDENCE } from './subworkflows/rnaseq_evidence.nf'
include { ISOSEQ_EVIDENCE } from './subworkflows/isoseq_evidence.nf'
include { HC_GENES } from './subworkflows/hc_genes.nf'
include { AB_INITIO } from './subworkflows/ab_initio.nf'
include { DRUSILLA_HC } from './subworkflows/drusilla.nf'
include { STRINGTIE_ASSEMBLE_MIX } from './modules/assembly.nf'

workflow {
  main:
    def OUT_CH = null
    def outdir = params.outdir ?: "results"
    file(outdir).mkdirs()

    // infer from params
    def hasPaired   = params.rnaseq_paired?.size()  > 0 || params.rnaseq_sra_paired?.size() > 0
    def hasSingle   = params.rnaseq_single?.size()  > 0 || params.rnaseq_sra_single?.size() > 0
    def hasIso      = params.isoseq?.size()         > 0 || params.isoseq_sra?.size() > 0
    def hasBAM      = asList(params.rnaseq_bam).size() > 0

    def proteinsList = []
    if( params.proteins ) {
      def rawList = (params.proteins instanceof List) ? params.proteins : [params.proteins]
      proteinsList = rawList.findAll { p -> p }
    }
    def odb12List = []
    if( params.odb12Partitions ) {
      def rawOdb = (params.odb12Partitions instanceof List) ? params.odb12Partitions : [params.odb12Partitions]
      odb12List = rawOdb.findAll { p -> p }
    }
    def hasProteins = proteinsList.size() > 0 || odb12List.size() > 0

    def genefinder    = resolveGenefinder(params)
    def genefinderRun = genefinderEnabled(params)

    def MODE = params.mode ?: inferMode(hasPaired, hasSingle, hasIso, hasBAM, hasProteins)
    // 'tiberius' is the historic name of the ab initio mode
    if( MODE in ['tiberius', 'vipsania'] ) MODE = 'abinitio'
    log.info "Running mode: ${MODE}"
    log.info "Gene finder : ${genefinderRun ? genefinder : 'none'}"
    def hc = hcMethod(params, MODE)
    if( hc.note ) log.warn hc.note
    def useDrusilla = hc.method == 'drusilla'
    if( MODE != 'abinitio' && MODE != 'proteins' ) log.info "HC genes    : ${hc.method}"

    def inp  = INPUTS(params)

    if( MODE == 'abinitio' ) {
      OUT_CH = AB_INITIO(inp.genome, params)

    } else {

      def pe = PROTEIN_EVIDENCE(inp.genome, inp.proteins, inp.score, params)

      def empty_file = EMPTY_FILE()

      // StringTie assemblies: one per BAM for TransDecoder. Drusilla uses one
      // assembly of all reads, whose TPM and coverage its pre-filter reads.
      def asm_mode = !useDrusilla ? 'per_sample' : (MODE == 'mixed' ? 'none' : 'merged')

      // Use fully-qualified Channel to avoid any name shadowing issues
      def re = (MODE in ['mixed','rnaseq']) \
        ? RNASEQ_EVIDENCE(inp.genome, params, asm_mode) \
        : [hints: empty_file, asm_gtf: nextflow.Channel.empty(), bam: nextflow.Channel.empty()]

      def ie = (MODE in ['mixed','isoseq']) \
        ? ISOSEQ_EVIDENCE(inp.genome, params, asm_mode) \
        : [hints: empty_file, asm_gtf: nextflow.Channel.empty(), bam: nextflow.Channel.empty()]

      def asm_gtf  = nextflow.Channel.empty()

      if( MODE == 'mixed' ) {
        asm_gtf  = useDrusilla ? STRINGTIE_ASSEMBLE_MIX(re.bam, ie.bam).gtf : re.asm_gtf.mix(ie.asm_gtf)
      } else if( MODE == 'rnaseq' ) {
        asm_gtf  = re.asm_gtf
      } else if( MODE == 'isoseq' ) {
        asm_gtf  = ie.asm_gtf
      }

      def all_hints = CONCAT_HINTS(pe.prot_hints, re.hints, ie.hints)

      def train_final
      def genefinder_final = pe.genefinder_gff
      if( useDrusilla ) {
        def dr = DRUSILLA_HC(asm_gtf, inp.genome, genefinder, pe.genefinder_gff, pe.scored_gff, pe.prot_hc_hints)
        train_final      = dr.orfs
        genefinder_final = dr.genefinder
      } else if( MODE in ['mixed','rnaseq','isoseq'] ) {
        train_final = HC_GENES(asm_gtf, inp.genome, pe.proteindb, pe.scored_gff)
      } else {
        train_final = HC_FORMAT_FILTER(pe.prot_traingff, inp.genome)
      }

      if( genefinderRun ) {
        MERGE_GENEFINDER_TRAIN(genefinder, genefinder_final, train_final)
        PROTEIN_FROM_GFF_FINAL(genefinder, MERGE_GENEFINDER_TRAIN.out.merged, inp.genome)
        // MERGE_GENEFINDER_TRAIN_PRIO(genefinder, pe.genefinder_gff, train_final)
      }

      OUT_CH = all_hints.hints
    }
}
