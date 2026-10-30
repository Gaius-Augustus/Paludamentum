nextflow.enable.dsl=2

include { CONCAT_HINTS; EMPTY_FILE } from './modules/util.nf'
include { PROTEIN_FROM_GFF } from './modules/genefinder.nf'
include { MERGE_BAMS as MERGE_BAMS_ALL } from './modules/rnaseq.nf'
include { inferMode; normalizeMode; resolveGenefinder; genefinderEnabled; asList; varusInputs; hcClade; cladeAdvice; truthy } from './lib_nf/functions.nf'
include { citationsText } from './lib_nf/citations.nf'
include { methodsText } from './lib_nf/methods.nf'

include { INPUTS } from './subworkflows/inputs.nf'
include { PROTEIN_EVIDENCE } from './subworkflows/protein_evidence.nf'
include { GENEFINDER as GENEFINDER_TIBERIUS } from './subworkflows/genefinder.nf'
include { RNASEQ_EVIDENCE } from './subworkflows/rnaseq_evidence.nf'
include { ISOSEQ_EVIDENCE } from './subworkflows/isoseq_evidence.nf'
include { AB_INITIO } from './subworkflows/ab_initio.nf'
include { HC_COMBINED } from './subworkflows/hc_combined.nf'
include { STRINGTIE_ASSEMBLE_MIX } from './modules/assembly.nf'
include { VARUS_INPUT as VARUS_INPUT_MIX } from './modules/varus.nf'
include { POSTPROCESS; buscoLineage } from './subworkflows/postprocess.nf'

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
    // StringTie assemblies of the user (GTF files or globs), used as they are
    def stringtieFiles = asList(params.stringtie).findAll { p -> p }.collectMany { p ->
      def f = file(p.toString())
      def found = (f instanceof List) ? f : [f]
      if( !found || found.any { g -> !g.exists() } ) error "stringtie: no such file: ${p}"
      return found
    }
    def hasStringtie = stringtieFiles.size() > 0

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

    def MODE = params.mode ? normalizeMode(params.mode) \
      : inferMode(hasPaired, hasSingle, hasIso || nIsoVarus > 0, hasBAM || nVarus > 0 || hasStringtie, hasProteins)
    // A forced mode without its inputs would end without the main outputs.
    // A StringTie assembly (stringtie) stands for the reads of any transcript mode.
    if( MODE in ['rnaseq', 'mixed'] && !(hasPaired || hasSingle || hasBAM || nVarus > 0 || hasStringtie) )
      error "params.mode '${MODE}' needs short reads (rnaseq_paired, rnaseq_single, rnaseq_bam, rnaseq_varus or rnaseq_sra_*) or a StringTie assembly (stringtie)."
    if( MODE in ['isoseq', 'mixed'] && !(hasIso || nIsoVarus > 0 || hasStringtie) )
      error "params.mode '${MODE}' needs Iso-Seq reads (isoseq, isoseq_sra or isoseq_varus) or a StringTie assembly (stringtie)."
    if( hasStringtie && MODE in ['abinitio', 'proteins'] )
      log.warn "stringtie is not used in mode '${MODE}'."
    if( MODE == 'proteins' && !hasProteins )
      error "params.mode '${MODE}' needs protein evidence (proteins or odb12Partitions)."
    // Evidence modes (rnaseq / isoseq / mixed) always require proteins now:
    // combine_genesets reads DIAMOND, miniprot and miniprothint.
    def transcriptMode = MODE in ['rnaseq', 'isoseq', 'mixed']
    if( transcriptMode && !hasProteins )
      error "params.mode '${MODE}' needs protein evidence (proteins or odb12Partitions): the unified HC-gene step runs DIAMOND, miniprot and miniprothint on the four gene sets."
    log.info "Running mode: ${MODE}"

    // Both gene finders always run in evidence modes; abinitio / proteins
    // modes still pick one via params.genefinder or <tool>.run.
    if( transcriptMode ) {
      log.info "Gene finders: tiberius + vipsania"
      if( !params.tiberius?.model_cfg && !params.tiberius?.result )
        error "params.tiberius.model_cfg is required (both gene finders run in evidence modes)."
      if( !params.vipsania?.model && !params.vipsania?.result )
        error "params.vipsania.model is required (both gene finders run in evidence modes)."
    } else {
      def gfRun = genefinderEnabled(params)
      def gfSel = gfRun ? resolveGenefinder(params) : null
      log.info "Gene finder : ${gfRun ? gfSel : 'none'}"
    }

    // Clade of the gene finder model (Tiberius target_species or Vipsania
    // model). Insecta / Lepidoptera triggers --drop-lone-single-exon in
    // combine_genesets. cladeAdvice warns when the run differs from the
    // recommended flow of the clade.
    def clade = hcClade(params)
    def advice = cladeAdvice(params, MODE, hasProteins, [clade: clade.name])
    if( advice ) log.warn advice

    // The combine_genesets flow needs ONE StringTie assembly of all reads
    // (with the cov and TPM attributes of StringTie). A given assembly
    // (stringtie) is that one and nothing is assembled; the reads or BAMs it
    // was made from may come with it (they give the coverage track and the
    // intron hints), pyVARUS directories may not (they bring an assembly of
    // their own).
    def useStringtie = hasStringtie && transcriptMode
    if( useStringtie && stringtieFiles.size() > 1 )
      error "stringtie: the unified HC flow needs one StringTie assembly of all reads, got ${stringtieFiles.size()} files. Pass the one assembly (stringtie --mix for short and Iso-Seq reads together), or let the pipeline assemble."
    if( useStringtie && (nVarus > 0 || nIsoVarus > 0 || hasMixVarus) )
      error "stringtie: a pyVARUS directory brings its own StringTie assembly. Pass either the StringTie GTF (with the reads or BAMs it was made from) or the pyVARUS directories, not both."

    // What this run does, for the references (citations.md) and the methods
    // text (methods.md) of the run: lib_nf/citations.nf, lib_nf/methods.nf
    def tibModel = params.tiberius?.model_cfg ? file(params.tiberius.model_cfg.toString()).baseName : null
    def vipModel = params.vipsania?.model?.toString()?.trim() ?: null
    def gfModel = transcriptMode ? (tibModel ?: vipModel) :
        (resolveGenefinder(params) == 'tiberius' ? tibModel : vipModel)
    def run = [
        // keys of citationKeys (lib_nf/citations.nf)
        mode: MODE,
        genefinder: transcriptMode ? 'both' : (genefinderEnabled(params) ? resolveGenefinder(params) : null),
        tiberiusModel: params.tiberius?.model_cfg,
        hc: transcriptMode ? 'combined' : null,
        orfFinder: transcriptMode ? 'td2' : null,
        odb12: odb12List.size() > 0,
        shortFastq: hasPaired || hasSingle, shortBam: hasBAM, shortVarus: nVarus > 0,
        isoFastq: hasIso, isoVarus: nIsoVarus > 0,
        busco: params.qc?.busco_lineage && truthy(params.qc?.busco),
        compleasm: params.qc?.busco_lineage && truthy(params.qc?.compleasm),
        omark: truthy(params.qc?.omark), gffcompare: params.qc?.reference_annotation as boolean,
        ncrna: truthy(params.ncrna?.run), lncrna: truthy(params.ncrna?.lncrna), fantasia: truthy(params.fantasia?.run),
        // further facts of the methods text (methodsText in lib_nf/methods.nf)
        model: gfModel, clade: clade.name,
        tiberiusResult: params.tiberius?.result ? file(params.tiberius.result.toString()).name : null,
        vipsaniaResult: params.vipsania?.result ? file(params.vipsania.result.toString()).name : null,
        finetune: truthy(params.vipsania?.finetune),
        proteinFiles: proteinsList.size(), odb12Partitions: odb12List,
        scoringMatrix: params.scoring_matrix ? file(params.scoring_matrix.toString()).name : null,
        shortPaired: hasPaired, shortSingle: hasSingle,
        shortSra: asList(params.rnaseq_sra_paired).size() > 0 || asList(params.rnaseq_sra_single).size() > 0,
        isoSra: asList(params.isoseq_sra).size() > 0,
        minAlignmentRate: params.min_alignment_rate ?: 80,
        stringtieFiles: stringtieFiles.size(), mixVarus: MODE == 'mixed' && transcriptMode && hasMixVarus,
        td2PredictArgs: params.td2_predict_args ?: null,
        sanityFilter: truthy(params.postprocess?.sanity_filter), utr: truthy(params.postprocess?.utr),
        maxUtrExtension: params.postprocess?.max_utr_extension,
        buscoLineage: buscoLineage(params.qc?.busco_lineage),
        geneSupport: truthy(params.qc?.gene_support), statistics: truthy(params.qc?.statistics),
        trnascanHighConfidence: truthy(params.ncrna?.trnascan_high_confidence),
        fantasiaMinScore: params.fantasia?.min_score,
        minLocus: (params.combine?.min_locus != null) ? params.combine.min_locus : 4,
    ]
    def version = workflow.manifest.version ?: 'unknown'
    file("${outdir}/citations.md").text = citationsText(run, version)
    file("${outdir}/methods.md").text = methodsText(run, version)

    // Mixed mode: one stringtie --mix assembly of both merged BAMs.
    // pyVARUS directories have no BAM; varus assemble makes that assembly.
    def mixVarus = []
    if( MODE == 'mixed' && transcriptMode ) {
      if( hasMixVarus ) {
        if( asList(params.mixed_varus).size() > 1 )
          error "mixed_varus takes one directory (the output of varus assemble --short --long), got ${asList(params.mixed_varus).size()}."
        if( nVarus != 1 || nIsoVarus != 1 || hasPaired || hasSingle || hasBAM || hasIso )
          error "mixed_varus is the one StringTie assembly of the unified HC flow, made from one short-read and one " +
                "Iso-Seq pyVARUS BAM. Pass exactly one rnaseq_varus and one isoseq_varus directory (the runs whose " +
                "BAMs varus assemble used) and no other short-read or Iso-Seq input."
        mixVarus = varusInputs(params.mixed_varus, 'mixed_varus', 'mixed')
      } else if( nVarus > 0 || nIsoVarus > 0 ) {
        error "The unified HC flow in mixed mode needs one StringTie assembly of the short reads and the Iso-Seq " +
              "reads together, and pyVARUS directories have no BAM to make it. Run " +
              "`varus assemble GENOME --short A/VARUS.bam --long B/VARUS.bam --outdir M` " +
              "(`varus replay` rebuilds a dropped BAM) and pass M as mixed_varus."
      }
    } else if( hasMixVarus ) {
      log.warn "mixed_varus is not used in mode '${MODE}'."
    }

    def inp  = INPUTS(params)

    // The final protein-coding annotation of either branch, and the evidence
    // that the post-processing uses (subworkflows/postprocess.nf)
    def final_gff = nextflow.Channel.empty()
    def asm_short = nextflow.Channel.empty()   // StringTie assemblies of short reads (UTRs)
    def asm_long  = nextflow.Channel.empty()   // StringTie assemblies of Iso-Seq reads (UTRs)
    def hints     = nextflow.Channel.empty()   // hintsfile.gff (gene support)
    def genome_chunks = nextflow.Channel.empty()   // genome chunks of the gene finder (cmscan)
    def genefinder_versions = nextflow.Channel.empty()   // version line of the gene finder
    def post_prefix = null                     // stem of the published files

    if( MODE == 'abinitio' ) {
      def ab = AB_INITIO(inp.genome, params)
      final_gff = ab.gff
      genome_chunks = ab.chunks
      genefinder_versions = ab.versions
      post_prefix = resolveGenefinder(params)

    } else if( MODE == 'proteins' ) {
      // Proteins-only mode: run one gene finder as today, feed its proteins
      // to PREPROCESS_PROTEINDB for species selection, publish the miniprot
      // training GFF as the HC gene set. No combine, no transcript evidence.
      def tool = resolveGenefinder(params)
      def gf = GENEFINDER_TIBERIUS(inp.genome, params, tool, false)
      genome_chunks = gf.chunks
      genefinder_versions = gf.versions
      def tib_prot = PROTEIN_FROM_GFF(tool, gf.gff, inp.genome)
      def pe = PROTEIN_EVIDENCE(inp.genome, inp.proteins, inp.score, tib_prot, params)
      // Nothing to merge: the miniprot training genes become the HC genes
      // (keep today's behaviour; the proteins-only path is unchanged).
      // Publish via POSTPROCESS prefixed with the gene finder name.
      final_gff = pe.prot_traingff
      post_prefix = tool
      // The protein hints as hintsfile.gff: the gene support of POSTPROCESS
      // (GENE_SUPPORT, which GENE_SET_STATISTICS waits for) reads them.
      def empty_file = EMPTY_FILE()
      hints = CONCAT_HINTS(pe.prot_hints, empty_file, empty_file).hints

    } else {

      // Evidence modes: rnaseq / isoseq / mixed. Unified HC flow.
      def empty_file = EMPTY_FILE()

      // Run Tiberius first; its proteins feed PREPROCESS_PROTEINDB (species
      // selection in large OrthoDB databases) and HC_COMBINED reuses the
      // same prediction (same task hash => Nextflow cache hit).
      def tib = GENEFINDER_TIBERIUS(inp.genome, params, 'tiberius', false)
      def tib_prot = PROTEIN_FROM_GFF('tiberius', tib.gff, inp.genome)

      def pe = PROTEIN_EVIDENCE(inp.genome, inp.proteins, inp.score, tib_prot, params)

      // One StringTie assembly of all reads; one merged BAM for coverage
      // (BAM2COV) and intron hints (BAM2HINTS). With a given assembly the
      // reads are mapped for the BAM only ('none': nothing is assembled).
      def asmMode = useStringtie ? 'none' : 'merged'
      def re = (MODE in ['mixed','rnaseq']) \
        ? RNASEQ_EVIDENCE(inp.genome, params, asmMode) \
        : [hints: empty_file, asm_gtf: nextflow.Channel.empty(), bam: nextflow.Channel.empty()]
      def ie = (MODE in ['mixed','isoseq']) \
        ? ISOSEQ_EVIDENCE(inp.genome, params, asmMode) \
        : [hints: empty_file, asm_gtf: nextflow.Channel.empty(), bam: nextflow.Channel.empty()]

      // Assembly channel
      def asm_gtf  = nextflow.Channel.empty()
      def stringtie_gtf = nextflow.Channel.fromList(stringtieFiles)
      if( useStringtie ) {
        asm_gtf  = stringtie_gtf
      } else if( MODE == 'mixed' ) {
        if( mixVarus ) asm_gtf = VARUS_INPUT_MIX(nextflow.Channel.fromList(mixVarus), nextflow.Channel.value(file(params.genome))).gtf
        else           asm_gtf = STRINGTIE_ASSEMBLE_MIX(re.bam, ie.bam).gtf
      } else if( MODE == 'rnaseq' ) {
        asm_gtf  = re.asm_gtf
      } else if( MODE == 'isoseq' ) {
        asm_gtf  = ie.asm_gtf
      }

      // One merged BAM of all reads for BAM2COV. In mixed mode merge both
      // evidence BAMs; otherwise use the one that was made. pyVARUS
      // directories and a given assembly without reads have no BAM: the
      // placeholder stands for the coverage track then, and
      // combine_genesets.py runs without --coverage (COMBINE_GENESETS).
      def haveBam = (MODE in ['mixed', 'rnaseq'] && (hasPaired || hasSingle || hasBAM)) ||
                    (MODE in ['mixed', 'isoseq'] && hasIso)
      def merged_bam
      if( !haveBam ) {
        merged_bam = empty_file
      } else if( MODE == 'mixed' ) {
        merged_bam = MERGE_BAMS_ALL(re.bam.mix(ie.bam)).bam
      } else if( MODE == 'rnaseq' ) {
        merged_bam = re.bam
      } else {
        merged_bam = ie.bam
      }

      // Intron hints for combine_genesets: the subworkflows produce these per
      // read type (BAM2HINTS, or the hints of a pyVARUS directory), an empty
      // placeholder without reads. combine_genesets.py --introns accepts
      // several files; COMBINE_GENESETS leaves empty ones out.
      def introns_combine
      if( MODE == 'mixed' )       introns_combine = re.hints.mix(ie.hints).collect()
      else if( MODE == 'rnaseq' ) introns_combine = re.hints
      else                        introns_combine = ie.hints

      def hc = HC_COMBINED(
          inp.genome, asm_gtf, merged_bam, haveBam, introns_combine,
          pe.proteindb, pe.scored_gff, pe.prot_traingff, pe.prot_hc_hints,
          tib.gff, tib_prot, params
      )

      final_gff = hc.gff
      genome_chunks = tib.chunks
      genefinder_versions = tib.versions.mix(hc.versions)
      post_prefix = 'combined'

      // UTRs: use the one StringTie assembly as asm_short (same as today's
      // Drusilla flow). Iso-Seq contributes to the one merged assembly, so
      // asm_long stays empty.
      asm_short = asm_gtf

      // Hints for POSTPROCESS (gene support): protein + RNA + Iso-Seq hints.
      def all_hints = CONCAT_HINTS(pe.prot_hints, re.hints, ie.hints)
      hints = all_hints.hints
    }

    // Sanity filter, UTRs, GTF and sequences, quality control, ncRNA and GO
    // terms of the final annotation. Evidence runs always have a final
    // annotation now (HC_COMBINED or the proteins-only miniprot training GFF).
    if( post_prefix ) {
      def runInfo = [
          version: version, mode: MODE, genefinder: post_prefix, model: gfModel,
          hc: MODE in ['abinitio', 'proteins'] ? null : 'combined',
      ]
      POSTPROCESS(final_gff, inp.genome, genome_chunks, post_prefix, genefinder_versions, MODE,
                  asm_short, asm_long, hints, runInfo)
    }
}
