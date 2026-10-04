nextflow.enable.dsl=2

include { CONCAT_HINTS; EMPTY_FILE } from './modules/util.nf'
include { MERGE_GENEFINDER_TRAIN; PROTEIN_FROM_GFF_FINAL } from './modules/genefinder.nf'
include { inferMode; normalizeMode; resolveGenefinder; genefinderEnabled; asList; hcMethod; varusInputs; orfFinder; drusillaSetting } from './lib_nf/functions.nf'
include { HC_FORMAT_FILTER } from './modules/hc.nf'

include { INPUTS } from './subworkflows/inputs.nf'
include { PROTEIN_EVIDENCE } from './subworkflows/protein_evidence.nf'
include { RNASEQ_EVIDENCE } from './subworkflows/rnaseq_evidence.nf'
include { ISOSEQ_EVIDENCE } from './subworkflows/isoseq_evidence.nf'
include { HC_GENES } from './subworkflows/hc_genes.nf'
include { AB_INITIO } from './subworkflows/ab_initio.nf'
include { DRUSILLA_HC } from './subworkflows/drusilla.nf'
include { STRINGTIE_ASSEMBLE_MIX } from './modules/assembly.nf'
include { VARUS_INPUT as VARUS_INPUT_MIX } from './modules/varus.nf'

workflow {
  main:
    def outdir = params.outdir ?: "results"
    file(outdir).mkdirs()

    // infer from params
    def hasPaired   = params.rnaseq_paired?.size()  > 0 || params.rnaseq_sra_paired?.size() > 0
    def hasSingle   = params.rnaseq_single?.size()  > 0 || params.rnaseq_sra_single?.size() > 0
    def hasIso      = params.isoseq?.size()         > 0 || params.isoseq_sra?.size() > 0
    def hasBAM      = asList(params.rnaseq_bam).size() > 0
    // pyVARUS output directories: short reads, Iso-Seq, varus assemble --mix
    def nVarus      = asList(params.rnaseq_varus).size()
    def nIsoVarus   = asList(params.isoseq_varus).size()
    def hasMixVarus = asList(params.mixed_varus).size() > 0

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

    def MODE = params.mode ? normalizeMode(params.mode) \
      : inferMode(hasPaired, hasSingle, hasIso || nIsoVarus > 0, hasBAM || nVarus > 0, hasProteins)
    // A forced mode without its inputs would end without the main outputs.
    if( MODE in ['rnaseq', 'mixed'] && !(hasPaired || hasSingle || hasBAM || nVarus > 0) )
      error "params.mode '${MODE}' needs short reads (rnaseq_paired, rnaseq_single, rnaseq_bam, rnaseq_varus or rnaseq_sra_*)."
    if( MODE in ['isoseq', 'mixed'] && !(hasIso || nIsoVarus > 0) )
      error "params.mode '${MODE}' needs Iso-Seq reads (isoseq, isoseq_sra or isoseq_varus)."
    if( MODE != 'abinitio' && !hasProteins )
      error "params.mode '${MODE}' needs protein evidence (proteins or odb12Partitions)."
    log.info "Running mode: ${MODE}"
    log.info "Gene finder : ${genefinderRun ? genefinder : 'none'}"
    def hc = hcMethod(params, MODE)
    if( hc.note ) log.warn hc.note
    def useDrusilla = hc.method == 'drusilla'
    if( MODE != 'abinitio' && MODE != 'proteins' ) {
      // conf/hc_genes.yaml (params.hc_table) chooses per clade of the gene finder model
      def clade = hc.clade ? "clade ${hc.clade}, " : ''
      log.info "HC genes    : ${hc.method} (" + clade +
               (useDrusilla ? "Drusilla model ${drusillaSetting(params, 'model') ?: params.drusilla.weights}" : "ORF finder ${orfFinder(params)}") + ")"
    }

    // Mixed mode with Drusilla: one stringtie --mix assembly of both BAMs.
    // pyVARUS directories have no BAM; varus assemble makes that assembly.
    def mixVarus = []
    if( MODE == 'mixed' && useDrusilla ) {
      if( hasMixVarus ) {
        if( asList(params.mixed_varus).size() > 1 )
          error "mixed_varus takes one directory (the output of varus assemble --short --long), got ${asList(params.mixed_varus).size()}."
        if( nVarus != 1 || nIsoVarus != 1 || hasPaired || hasSingle || hasBAM || hasIso )
          error "mixed_varus is the one StringTie assembly of the Drusilla flow, made from one short-read and one " +
                "Iso-Seq pyVARUS BAM. Pass exactly one rnaseq_varus and one isoseq_varus directory (the runs whose " +
                "BAMs varus assemble used) and no other short-read or Iso-Seq input."
        mixVarus = varusInputs(params.mixed_varus, 'mixed_varus', 'mixed')
      } else if( nVarus > 0 || nIsoVarus > 0 ) {
        error "The Drusilla flow in mixed mode needs one StringTie assembly of the short reads and the Iso-Seq " +
              "reads together, and pyVARUS directories have no BAM to make it. Run " +
              "`varus assemble GENOME --short A/VARUS.bam --long B/VARUS.bam --outdir M` " +
              "(`varus replay` rebuilds a dropped BAM) and pass M as mixed_varus."
      }
    } else if( hasMixVarus ) {
      log.warn "mixed_varus is not used: it serves only the Drusilla flow in mixed mode " +
               "(this run: mode ${MODE}${MODE == 'mixed' ? ', TransDecoder HC genes' : ''})."
    }

    def inp  = INPUTS(params)

    if( MODE == 'abinitio' ) {
      AB_INITIO(inp.genome, params)

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
        if( !useDrusilla )  asm_gtf = re.asm_gtf.mix(ie.asm_gtf)
        else if( mixVarus ) asm_gtf = VARUS_INPUT_MIX(nextflow.Channel.fromList(mixVarus), nextflow.Channel.value(file(params.genome))).gtf
        else                asm_gtf = STRINGTIE_ASSEMBLE_MIX(re.bam, ie.bam).gtf
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
        train_final = HC_GENES(asm_gtf, inp.genome, pe.proteindb)
      } else {
        train_final = HC_FORMAT_FILTER(pe.prot_traingff, inp.genome)
      }

      if( genefinderRun ) {
        MERGE_GENEFINDER_TRAIN(genefinder, genefinder_final, train_final)
        PROTEIN_FROM_GFF_FINAL(genefinder, MERGE_GENEFINDER_TRAIN.out.merged, inp.genome)
      }
    }
}
