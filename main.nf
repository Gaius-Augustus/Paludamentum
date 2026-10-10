nextflow.enable.dsl=2

include { CONCAT_HINTS; EMPTY_FILE } from './modules/util.nf'
include { MERGE_GENEFINDER_TRAIN } from './modules/genefinder.nf'
include { inferMode; normalizeMode; resolveGenefinder; genefinderEnabled; asList; hcMethod; rescueEnabled; rescueModel; varusInputs; orfFinder; drusillaSetting; truthy } from './lib_nf/functions.nf'
include { citationsText } from './lib_nf/citations.nf'
include { methodsText } from './lib_nf/methods.nf'
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

    def genefinder    = resolveGenefinder(params)
    def genefinderRun = genefinderEnabled(params)

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
    // The Drusilla flow filters one assembly of all reads by its coverage and
    // TPM, so a given assembly must be that one assembly.
    def stringtieOnly = hasStringtie && useDrusilla
    if( stringtieOnly && (stringtieFiles.size() > 1 || hasPaired || hasSingle || hasBAM || hasIso
                          || nVarus > 0 || nIsoVarus > 0 || hasMixVarus) )
      error "stringtie: the Drusilla flow (clades with hc: drusilla in hc_table) needs one StringTie assembly of all reads, " +
            "with the coverage and TPM that StringTie writes. Pass a single StringTie GTF and no other RNA-Seq or " +
            "Iso-Seq input, or set drusilla.run = false to merge several assemblies in the TransDecoder flow."

    // What this run does, for the references (citations.md) and the methods
    // text (methods.md) of the run: lib_nf/citations.nf, lib_nf/methods.nf
    def gfParams = params[genefinder] ?: [:]
    def gfModel = genefinder == 'tiberius' ?
        (params.tiberius?.model_cfg ? file(params.tiberius.model_cfg.toString()).baseName : null) :
        params.vipsania?.model?.toString()?.trim() ?: null
    def drp = params.drusilla ?: [:]
    def drusillaSettings = null
    if( useDrusilla ) {
        def fixStop = drp.fix_stop == null || truthy(drp.fix_stop)
        def lgbModel = drusillaSetting(params, 'lgb_model')?.toString()
        def rm = rescueModel(params)
        drusillaSettings = [
            model: drp.weights ? null : drusillaSetting(params, 'model'),
            weights: drp.weights ? file(drp.weights.toString()).name : null,
            lgbModel: lgbModel ? lgbModel.tokenize('/').last().replaceFirst(/\.(tar\.gz|tgz|txt)$/, '') : null,
            minLength: drp.min_length, minCov: drp.min_cov, minTpm: drp.min_tpm,
            longLength: drp.long_length, minTpmLong: drp.min_tpm_long, minCodingLength: drp.min_coding_length,
            fixStop: fixStop, fixStart: fixStop && (drp.fix_start == null || truthy(drp.fix_start)),
            lgbThreshold: drp.lgb_threshold, lgbKeep: (drp.lgb_keep ?: 'correct').toString().tokenize(', \t'),
            rescueModel: rm.name ?: file(rm.file.toString()).baseName,
            rescueFlank: drp.rescue_flank, rescueHintWeight: drp.rescue_hint_weight,
            rescueOrfFilter: truthy(drp.rescue_orf_filter),
        ]
    }
    def run = [
        // keys of citationKeys (lib_nf/citations.nf)
        mode: MODE, genefinder: genefinderRun ? genefinder : null, tiberiusModel: params.tiberius?.model_cfg,
        hc: hc.method, orfFinder: hc.method == 'transdecoder' ? orfFinder(params) : null,
        rescue: rescueEnabled(params), odb12: odb12List.size() > 0,
        shortFastq: hasPaired || hasSingle, shortBam: hasBAM, shortVarus: nVarus > 0,
        isoFastq: hasIso, isoVarus: nIsoVarus > 0,
        busco: params.qc?.busco_lineage && truthy(params.qc?.busco),
        compleasm: params.qc?.busco_lineage && truthy(params.qc?.compleasm),
        omark: truthy(params.qc?.omark), gffcompare: params.qc?.reference_annotation as boolean,
        ncrna: truthy(params.ncrna?.run), lncrna: truthy(params.ncrna?.lncrna), fantasia: truthy(params.fantasia?.run),
        // further facts of the methods text (methodsText in lib_nf/methods.nf)
        model: gfModel, clade: hc.clade,
        result: genefinderRun && gfParams.result ? file(gfParams.result.toString()).name : null,
        finetune: genefinder == 'vipsania' && truthy(params.vipsania?.finetune),
        splitMinSize: gfParams.min_split_size ?: 20000000, splitMaxFiles: gfParams.max_files ?: 20,
        proteinFiles: proteinsList.size(), odb12Partitions: odb12List,
        scoringMatrix: params.scoring_matrix ? file(params.scoring_matrix.toString()).name : null,
        shortPaired: hasPaired, shortSingle: hasSingle,
        shortSra: asList(params.rnaseq_sra_paired).size() > 0 || asList(params.rnaseq_sra_single).size() > 0,
        isoSra: asList(params.isoseq_sra).size() > 0,
        minAlignmentRate: params.min_alignment_rate ?: 80,
        stringtieFiles: stringtieFiles.size(), mixVarus: MODE == 'mixed' && useDrusilla && hasMixVarus,
        td2PredictArgs: params.td2_predict_args ?: null, drusillaSettings: drusillaSettings,
        sanityFilter: truthy(params.postprocess?.sanity_filter), utr: truthy(params.postprocess?.utr),
        maxUtrExtension: params.postprocess?.max_utr_extension,
        buscoLineage: buscoLineage(params.qc?.busco_lineage),
        geneSupport: truthy(params.qc?.gene_support), statistics: truthy(params.qc?.statistics),
        trnascanHighConfidence: truthy(params.ncrna?.trnascan_high_confidence),
        fantasiaMinScore: params.fantasia?.min_score,
    ]
    def version = workflow.manifest.version ?: 'unknown'
    file("${outdir}/citations.md").text = citationsText(run, version)
    file("${outdir}/methods.md").text = methodsText(run, version)

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

    // The final protein-coding annotation of either branch, and the evidence
    // that the post-processing uses (subworkflows/postprocess.nf)
    def final_gff = nextflow.Channel.empty()
    def asm_short = nextflow.Channel.empty()   // StringTie assemblies of short reads (UTRs)
    def asm_long  = nextflow.Channel.empty()   // StringTie assemblies of Iso-Seq reads (UTRs)
    def hints     = nextflow.Channel.empty()   // hintsfile.gff (gene support)
    def genome_chunks = nextflow.Channel.empty()   // genome chunks of the gene finder (cmscan)
    def genefinder_versions = nextflow.Channel.empty()   // version line of the gene finder

    if( MODE == 'abinitio' ) {
      def ab = AB_INITIO(inp.genome, params)
      final_gff = ab.gff
      genome_chunks = ab.chunks
      genefinder_versions = ab.versions

    } else {

      def pe = PROTEIN_EVIDENCE(inp.genome, inp.proteins, inp.score, params)
      genome_chunks = pe.genome_chunks
      genefinder_versions = pe.genefinder_versions

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
      def stringtie_gtf = nextflow.Channel.fromList(stringtieFiles)

      if( stringtieOnly ) {
        asm_gtf  = stringtie_gtf
      } else if( MODE == 'mixed' ) {
        if( !useDrusilla )  asm_gtf = re.asm_gtf.mix(ie.asm_gtf)
        else if( mixVarus ) asm_gtf = VARUS_INPUT_MIX(nextflow.Channel.fromList(mixVarus), nextflow.Channel.value(file(params.genome))).gtf
        else                asm_gtf = STRINGTIE_ASSEMBLE_MIX(re.bam, ie.bam).gtf
      } else if( MODE == 'rnaseq' ) {
        asm_gtf  = re.asm_gtf
      } else if( MODE == 'isoseq' ) {
        asm_gtf  = ie.asm_gtf
      }
      // TransDecoder flow: the given assemblies are merged with those of the reads
      if( hasStringtie && !stringtieOnly ) asm_gtf = asm_gtf.mix(stringtie_gtf)

      // UTRs: the Iso-Seq assemblies of the TransDecoder flow are long-read
      // evidence; one assembly of all reads (Drusilla flow) counts as short reads
      if( useDrusilla || stringtieOnly ) {
        asm_short = asm_gtf
      } else {
        asm_short = (MODE in ['mixed','rnaseq']) ? re.asm_gtf : nextflow.Channel.empty()
        asm_long  = (MODE in ['mixed','isoseq']) ? ie.asm_gtf : nextflow.Channel.empty()
        if( hasStringtie ) asm_short = asm_short.mix(stringtie_gtf)
      }

      def all_hints = CONCAT_HINTS(pe.prot_hints, re.hints, ie.hints)
      hints = all_hints.hints

      // The HC genes are evidence-backed gene models that MERGE_GENEFINDER_TRAIN
      // below merges into the final annotation; they are not a training set.
      // No gene finder is trained in this pipeline: Tiberius takes
      // tiberius.model_cfg and Vipsania downloads a pretrained model. The
      // 'TRAIN' of the process name, and training.gff of the HC flow, are
      // inherited from BRAKER, where the HC genes do train AUGUSTUS.
      def hc_gff
      def genefinder_final = pe.genefinder_gff
      if( useDrusilla ) {
        def dr = DRUSILLA_HC(asm_gtf, inp.genome, genefinder, pe.genefinder_gff, pe.scored_gff, pe.prot_hc_hints)
        hc_gff           = dr.orfs
        genefinder_final = dr.genefinder
      } else if( MODE in ['mixed','rnaseq','isoseq'] ) {
        hc_gff = HC_GENES(asm_gtf, inp.genome, pe.proteindb)
      } else {
        hc_gff = HC_FORMAT_FILTER(pe.prot_traingff, inp.genome)
      }

      if( genefinderRun ) {
        final_gff = MERGE_GENEFINDER_TRAIN(genefinder, genefinder_final, hc_gff).merged
      }
    }

    // Sanity filter, UTRs, GTF and sequences, quality control, ncRNA and GO
    // terms of the final annotation. Evidence runs without a gene finder have
    // no final annotation and end with the hints and the HC genes.
    if( genefinderRun ) {
      def runInfo = [
          version: version, mode: MODE, genefinder: genefinder, model: gfModel,
          hc: MODE in ['abinitio', 'proteins'] ? null : hc.method,
      ]
      POSTPROCESS(final_gff, inp.genome, genome_chunks, genefinder, genefinder_versions, MODE, asm_short, asm_long,
                  hints, runInfo)
    }
}
