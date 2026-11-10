nextflow.enable.dsl=2

include { CONCAT_HINTS; EMPTY_FILE } from './modules/util.nf'
include { PROTEIN_FROM_GFF } from './modules/genefinder.nf'
include { MERGE_BAMS as MERGE_BAMS_ALL } from './modules/rnaseq.nf'
include { inferMode; normalizeMode; resolveGenefinder; genefinderEnabled; asList; varusInputs; hcClade; classifierUnknownClade; truthy } from './lib_nf/functions.nf'
include { citationsText } from './lib_nf/citations.nf'
include { methodsText } from './lib_nf/methods.nf'

include { INPUTS } from './subworkflows/inputs.nf'
include { PROTEIN_EVIDENCE } from './subworkflows/protein_evidence.nf'
include { GENEFINDER as GENEFINDER_TIBERIUS; GENEFINDER as GENEFINDER_VIPSANIA_FIRST } from './subworkflows/genefinder.nf'
include { RNASEQ_EVIDENCE } from './subworkflows/rnaseq_evidence.nf'
include { ISOSEQ_EVIDENCE } from './subworkflows/isoseq_evidence.nf'
include { AB_INITIO } from './subworkflows/ab_initio.nf'
include { HC_COMBINED } from './subworkflows/hc_combined.nf'
include { STRINGTIE_ASSEMBLE_MIX; STRINGTIE_ASSEMBLE_ISO as STRINGTIE_ASSEMBLE_CHAINS; ASSEMBLY_HINTS } from './modules/assembly.nf'
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
    // species kept of a large protein database (PREPROCESS_PROTEINDB); 0 keeps it whole
    def proteinTopSpecies = params.protein_top_species != null ? params.protein_top_species : 13
    if( !(proteinTopSpecies.toString().isInteger() && (proteinTopSpecies as int) >= 0) )
      error "protein_top_species: a whole number >= 0, not '${proteinTopSpecies}'."
    // Evidence modes (rnaseq / isoseq / mixed) always require proteins now:
    // combine_genesets reads DIAMOND, miniprot and miniprothint.
    def transcriptMode = MODE in ['rnaseq', 'isoseq', 'mixed']
    if( transcriptMode && !hasProteins )
      error "params.mode '${MODE}' needs protein evidence (proteins or odb12Partitions): the unified HC-gene step runs DIAMOND, miniprot and miniprothint on the four gene sets."
    // The protein-only mode runs the same unified HC step (both gene finders,
    // protein evidence, combine_genesets) without transcript inputs: no
    // assembly, no ORF callers, no coverage track, no intron hints. Mode
    // abinitio with both gene finders given runs it without any evidence (the
    // ab initio combination, classifier_abinitio); with one gene finder it
    // publishes that prediction as before.
    def tibGiven = (params.tiberius?.model_cfg || params.tiberius?.result) as boolean
    def vipGiven = (params.vipsania?.model || params.vipsania?.result) as boolean
    def abBoth = MODE == 'abinitio' && tibGiven && vipGiven
    def combinedMode = transcriptMode || MODE == 'proteins' || abBoth
    log.info "Running mode: ${MODE}"

    // Both gene finders run in evidence modes and in mode proteins; mode
    // abinitio picks one via params.genefinder or <tool>.run. Without a
    // Tiberius model (no tiberius.model_cfg and no tiberius.result: a species
    // outside the Tiberius clades, e.g. Nematoda, Arachnida) an evidence mode
    // runs Vipsania as the only ab initio gene set, with the classifier trained
    // for that (combine.classifier_vipsania). Measured on 31 benchmark species:
    // 57.4 transcript / 71.2 gene F1, Vipsania alone 44.2 / 54.8, with Tiberius
    // 60.1 / 74.6. Mode proteins needs Tiberius: Vipsania with protein evidence
    // alone ends below Vipsania alone (no agreement signal, no ORF callers).
    def runTiberius = combinedMode && tibGiven
    if( combinedMode ) {
      if( !params.vipsania?.model && !params.vipsania?.result )
        error "params.vipsania.model is required (Vipsania runs in every evidence mode and in mode proteins)."
      if( !runTiberius && MODE == 'proteins' )
        error "params.mode 'proteins' needs a Tiberius model (tiberius.model_cfg or tiberius.result): with Vipsania " +
              "alone and no transcript evidence the combined selection falls below Vipsania alone. Give a Tiberius " +
              "model, or RNA-Seq / Iso-Seq reads, a BAM, pyVARUS directories or a StringTie assembly."
      if( runTiberius )
        log.info "Gene finders: tiberius + vipsania" + (abBoth ? " (ab initio combination: both models given, no evidence)" : "")
      else
        log.info "Gene finders: vipsania (no Tiberius model given: tiberius.model_cfg / tiberius.result; " +
                 "the combined step selects with combine.classifier_vipsania)"
    } else {
      def gfRun = genefinderEnabled(params)
      def gfSel = gfRun ? resolveGenefinder(params) : null
      log.info "Gene finder : ${gfRun ? gfSel : 'none'}"
    }

    // Clade of the gene finder model (Tiberius target_species or Vipsania
    // model). Insecta / Lepidoptera triggers --drop-lone-single-exon in
    // combine_genesets; the flow itself is the same for every clade.
    def clade = hcClade(params)

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
    // Reads or BAMs give the one merged BAM of the coverage track (BAM2COV) and
    // the intron hints (BAM2HINTS); pyVARUS directories and a given assembly
    // without reads bring none. With a given assembly alone, ASSEMBLY_HINTS
    // derives both from the assembly: intron hints with the summed StringTie
    // cov of the transcripts sharing an intron as multiplicity, and a coverage
    // track from the exon cov values. Measured on 31 benchmark species against
    // the BAM: -0.4 transcript / -0.2 gene F1 (no RNA evidence: -4.3 / -3.7).
    def haveReadsBam = (MODE in ['mixed', 'rnaseq'] && (hasPaired || hasSingle || hasBAM)) ||
                       (MODE in ['mixed', 'isoseq'] && hasIso)
    def asmHints = useStringtie && !haveReadsBam

    // Transcripts of full-length reads for combine_genesets.py --chains (GTF
    // or GFF3 with exon lines, e.g. a StringTie -L assembly of long reads).
    // Given with `chains`; in mixed mode the pipeline otherwise makes them
    // from the Iso-Seq reads alone (docs/parameters.md, "Long-read transcripts").
    def chainFiles = asList(params.chains).findAll { p -> p }.collectMany { p ->
      def f = file(p.toString())
      def found = (f instanceof List) ? f : [f]
      if( !found || found.any { g -> !g.exists() } ) error "chains: no such file: ${p}"
      return found
    }
    if( chainFiles && !transcriptMode )
      log.warn "chains is not used in mode '${MODE}'."
    def autoChains = MODE == 'mixed' && (hasIso || nIsoVarus > 0)

    // What this run does, for the references (citations.md) and the methods
    // text (methods.md) of the run: lib_nf/citations.nf, lib_nf/methods.nf
    def tibModel = params.tiberius?.model_cfg ? file(params.tiberius.model_cfg.toString()).baseName : null
    def vipModel = params.vipsania?.model?.toString()?.trim() ?: null
    def gfModel = combinedMode ? (tibModel ?: vipModel) :
        (resolveGenefinder(params) == 'tiberius' ? tibModel : vipModel)
    // thresholds of the selection model of combine_genesets.py, read from its JSON file: the model of
    // the evidence modes (combine.classifier) or of the protein-only mode (combine.classifier_proteins)
    def clfParam = transcriptMode ? (runTiberius ? params.combine?.classifier : params.combine?.classifier_vipsania) :
                   MODE == 'proteins' ? params.combine?.classifier_proteins :
                   abBoth ? params.combine?.classifier_abinitio : null
    def clfFile = clfParam ? file(clfParam.toString()) : null
    def clfModel = clfFile?.exists() ? new groovy.json.JsonSlurper().parse(clfFile.toFile()) : null
    def run = [
        // keys of citationKeys (lib_nf/citations.nf)
        mode: MODE,
        genefinder: combinedMode ? (runTiberius ? 'both' : 'vipsania') : (genefinderEnabled(params) ? resolveGenefinder(params) : null),
        tiberiusModel: params.tiberius?.model_cfg,
        hc: combinedMode ? 'combined' : null,
        // Tiberius in the combined gene sets (false: Vipsania is the only ab initio gene set)
        hcTiberius: combinedMode ? runTiberius : null,
        orfFinder: transcriptMode ? 'td2' : null,
        odb12: odb12List.size() > 0,
        shortFastq: hasPaired || hasSingle, shortBam: hasBAM, shortVarus: nVarus > 0,
        isoFastq: hasIso, isoVarus: nIsoVarus > 0,
        // full-length transcripts of the integration: 'given' (chains), 'assembled' (mixed mode), or null
        chains: !transcriptMode ? null : (chainFiles ? 'given' : (autoChains ? 'assembled' : null)),
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
        proteinTopSpecies: proteinTopSpecies as int,
        scoringMatrix: params.scoring_matrix ? file(params.scoring_matrix.toString()).name : null,
        shortPaired: hasPaired, shortSingle: hasSingle,
        shortSra: asList(params.rnaseq_sra_paired).size() > 0 || asList(params.rnaseq_sra_single).size() > 0,
        isoSra: asList(params.isoseq_sra).size() > 0,
        minAlignmentRate: params.min_alignment_rate ?: 80,
        // the strand fix runs on the BAMs that StringTie assembles (not with a given assembly)
        strandFix: !hasStringtie && (params.strand_fix == null || truthy(params.strand_fix)),
        // both strands of '.' transcripts only matter when FILTER_STRINGTIE keeps them
        drusillaBothStrands: truthy(params.drusilla?.unstranded_both_strands) &&
            !truthy(params.drusilla?.drop_unstranded),
        drusillaFlank: Math.max((params.drusilla?.flank_bp ?: 0) as Integer, (params.drusilla?.flank_3p_bp ?: 0) as Integer),
        drusillaOrfsPerTx: params.drusilla?.orfs_per_tx,
        stringtieFiles: stringtieFiles.size(), mixVarus: MODE == 'mixed' && transcriptMode && hasMixVarus,
        // intron hints and coverage track derived from a given assembly without reads (ASSEMBLY_HINTS)
        asmHints: asmHints,
        td2PredictArgs: params.td2_predict_args ?: null,
        sanityFilter: truthy(params.postprocess?.sanity_filter), utr: truthy(params.postprocess?.utr),
        maxUtrExtension: params.postprocess?.max_utr_extension,
        buscoLineage: buscoLineage(params.qc?.busco_lineage),
        geneSupport: truthy(params.qc?.gene_support), statistics: truthy(params.qc?.statistics),
        trnascanHighConfidence: truthy(params.ncrna?.trnascan_high_confidence),
        fantasiaMinScore: params.fantasia?.min_score,
        minLocus: (params.combine?.min_locus != null) ? params.combine.min_locus : 4,
        // learned selection of combine_genesets.py: the thresholds the run applies (the model file's,
        // the conservative locus threshold for a clade new to the classifier, or params.combine.*)
        classifier: clfModel != null,
        clfTLocus: params.combine?.t_locus != null ? params.combine.t_locus :
                   (clfModel == null ? null : (!classifierUnknownClade(clade) ? clfModel.t_locus :
                                               (clfModel.t_locus_unknown_clade ?: clfModel.t_locus))),
        clfTIso: params.combine?.t_iso != null ? params.combine.t_iso : clfModel?.t_iso,
        clfUnknownClade: clfModel != null && params.combine?.t_locus == null && classifierUnknownClade(clade),
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
    def post_prefix = null                     // gene finder of the run for the report: tiberius, vipsania or combined

    if( abBoth ) {
      // Ab initio combination: both gene finders given and no evidence. Both
      // run, and HC_COMBINED integrates the two predictions with the
      // placeholder for every evidence input (no protein database, no DIAMOND,
      // no miniprot, no transcripts), selecting with the classifier trained on
      // evidence-free pools. Measured on 31 benchmark species against Tiberius
      // alone: +1.5 transcript / +4.8 gene F1, precision-driven (sensitivity
      // falls by 3 to 6 points); negative for vertebrates at transcript level.
      def empty_file = EMPTY_FILE()
      def tib = GENEFINDER_TIBERIUS(inp.genome, params, 'tiberius', false)
      def hc = HC_COMBINED(
          inp.genome, nextflow.Channel.empty(), empty_file, false, empty_file, empty_file,
          empty_file, empty_file, empty_file, empty_file,
          tib.gff, empty_file, params, false, true, false, nextflow.Channel.empty(), nextflow.Channel.empty(), false
      )
      final_gff = hc.gff
      genome_chunks = tib.chunks
      genefinder_versions = tib.versions.mix(hc.versions)
      post_prefix = 'combined'

    } else if( MODE == 'abinitio' ) {
      def ab = AB_INITIO(inp.genome, params)
      final_gff = ab.gff
      genome_chunks = ab.chunks
      genefinder_versions = ab.versions
      post_prefix = resolveGenefinder(params)

    } else if( MODE == 'proteins' ) {
      // Protein-only mode: both gene finders and the protein evidence as in the
      // evidence modes, then the unified HC step without transcript inputs
      // (HC_COMBINED with transcripts = false): the Tiberius and Vipsania gene
      // sets, their DIAMOND tables and the miniprot / miniprothint evidence go
      // to combine_genesets.py; the empty placeholder stands for the
      // ORF-caller gene sets, the BAM, the intron hints and the chains.
      // Measured on 31 benchmark species against Tiberius alone: +5.4
      // transcript / +8.9 gene F1 (the RNA-Seq modes reach +12.2 / +15.4).
      def empty_file = EMPTY_FILE()
      def tib = GENEFINDER_TIBERIUS(inp.genome, params, 'tiberius', false)
      def tib_prot = PROTEIN_FROM_GFF('tiberius', tib.gff, inp.genome)
      def pe = PROTEIN_EVIDENCE(inp.genome, inp.proteins, inp.score, tib_prot, params)
      def hc = HC_COMBINED(
          inp.genome, nextflow.Channel.empty(), empty_file, false, empty_file, empty_file,
          pe.proteindb, pe.scored_gff, pe.prot_traingff, pe.prot_hc_hints,
          tib.gff, tib_prot, params, false, true, false, nextflow.Channel.empty(), nextflow.Channel.empty(), true
      )
      final_gff = hc.gff
      genome_chunks = tib.chunks
      genefinder_versions = tib.versions.mix(hc.versions)
      post_prefix = 'combined'
      // The protein hints as hintsfile.gff: the gene support of POSTPROCESS
      // (GENE_SUPPORT, which GENE_SET_STATISTICS waits for) reads them.
      hints = CONCAT_HINTS(pe.prot_hints, empty_file, empty_file).hints

    } else {

      // Evidence modes: rnaseq / isoseq / mixed. Unified HC flow.
      def empty_file = EMPTY_FILE()

      // Run Tiberius first; its proteins feed PREPROCESS_PROTEINDB (species
      // selection in large OrthoDB databases) and HC_COMBINED reuses the
      // same prediction (same task hash => Nextflow cache hit). Without a
      // Tiberius model Vipsania takes that place: it runs here, its proteins
      // select the species, and HC_COMBINED gets it passed in (vip_given).
      def ab = runTiberius ? GENEFINDER_TIBERIUS(inp.genome, params, 'tiberius', false) :
                             GENEFINDER_VIPSANIA_FIRST(inp.genome, params, 'vipsania', false)
      def ab_prot = PROTEIN_FROM_GFF(runTiberius ? 'tiberius' : 'vipsania', ab.gff, inp.genome)

      def pe = PROTEIN_EVIDENCE(inp.genome, inp.proteins, inp.score, ab_prot, params)

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
      // evidence BAMs; otherwise use the one that was made. Without a BAM the
      // slot carries the coverage track to use as it is: the bedGraph that
      // ASSEMBLY_HINTS derives from a given assembly, or, for pyVARUS
      // directories, the placeholder (combine_genesets.py then runs without
      // --coverage, COMBINE_GENESETS).
      def asm_ev = asmHints ? ASSEMBLY_HINTS(stringtie_gtf) : null
      def haveBam = haveReadsBam
      def merged_bam
      if( !haveBam ) {
        merged_bam = asmHints ? asm_ev.bg : empty_file
      } else if( MODE == 'mixed' ) {
        merged_bam = MERGE_BAMS_ALL(re.bam.mix(ie.bam)).bam
      } else if( MODE == 'rnaseq' ) {
        merged_bam = re.bam
      } else {
        merged_bam = ie.bam
      }

      // Intron hints for combine_genesets: the subworkflows produce these per
      // read type (BAM2HINTS, or the hints of a pyVARUS directory), an empty
      // placeholder without reads; with a given assembly alone they are the
      // introns of the assembly (ASSEMBLY_HINTS). combine_genesets.py
      // --introns accepts several files; COMBINE_GENESETS leaves empty ones out.
      def introns_combine
      if( asmHints )              introns_combine = asm_ev.hints
      else if( MODE == 'mixed' )  introns_combine = re.hints.mix(ie.hints).collect()
      else if( MODE == 'rnaseq' ) introns_combine = re.hints
      else                        introns_combine = ie.hints

      // Full-length transcripts for combine_genesets.py --chains: the given
      // files (chains), else in mixed mode the StringTie -L assembly of the
      // Iso-Seq reads alone (ISOSEQ_EVIDENCE makes it next to the --mix
      // assembly, an isoseq_varus directory brings it; with a given assembly
      // it is made here from the Iso-Seq BAM). Measured on H. sapiens:
      // transcript F1 +0.7 over the --mix assembly alone. The placeholder
      // stands for no chains; COMBINE_GENESETS leaves empty files out.
      def chains_combine
      if( chainFiles )
        chains_combine = nextflow.Channel.fromList(chainFiles).mix(empty_file).collect()
      else if( MODE == 'mixed' && useStringtie && hasIso )
        chains_combine = STRINGTIE_ASSEMBLE_CHAINS(ie.bam).gtf.mix(empty_file).collect()
      else if( MODE == 'mixed' )
        chains_combine = ie.asm_gtf.mix(empty_file).collect()
      else
        chains_combine = empty_file

      // Tiberius slots: the prediction and its proteins, or the empty placeholder
      // without a Tiberius model (then Vipsania goes in through the vip_* slots)
      def hc = HC_COMBINED(
          inp.genome, asm_gtf, merged_bam, haveBam, introns_combine, chains_combine,
          pe.proteindb, pe.scored_gff, pe.prot_traingff, pe.prot_hc_hints,
          runTiberius ? ab.gff : empty_file, runTiberius ? ab_prot : empty_file, params, true,
          runTiberius, !runTiberius,
          runTiberius ? nextflow.Channel.empty() : ab.gff, runTiberius ? nextflow.Channel.empty() : ab_prot, true
      )

      final_gff = hc.gff
      genome_chunks = ab.chunks
      genefinder_versions = ab.versions.mix(hc.versions)
      post_prefix = 'combined'

      // UTRs: use the one StringTie assembly as asm_short (same as today's
      // Drusilla flow). Iso-Seq contributes to the one merged assembly, so
      // asm_long stays empty.
      asm_short = asm_gtf

      // Hints for POSTPROCESS (gene support): protein + RNA + Iso-Seq hints;
      // the intron hints of a given assembly take the RNA slot.
      def all_hints = CONCAT_HINTS(pe.prot_hints, asmHints ? asm_ev.hints : re.hints, ie.hints)
      hints = all_hints.hints
    }

    // Sanity filter, UTRs, GTF and sequences, quality control, ncRNA and GO
    // terms of the final annotation. Every run has a final annotation: the
    // ab initio prediction, or combined.gff3 of HC_COMBINED (evidence modes
    // and mode proteins).
    if( post_prefix ) {
      // The data the run used, one line of the report's run summary
      def evidence = []
      if( hasProteins ) evidence << 'proteins'
      if( transcriptMode ) {
        if( MODE in ['rnaseq', 'mixed'] ) {
          if( hasPaired || hasSingle ) evidence << 'short reads'
          if( hasBAM ) evidence << 'short-read BAM'
          if( nVarus > 0 ) evidence << 'pyVARUS short reads'
        }
        if( MODE in ['isoseq', 'mixed'] ) {
          if( hasIso ) evidence << 'Iso-Seq reads'
          if( nIsoVarus > 0 ) evidence << 'pyVARUS Iso-Seq reads'
        }
        if( hasMixVarus ) evidence << 'pyVARUS mixed assembly'
        if( useStringtie )
          evidence << (asmHints ? 'StringTie assembly (intron hints and coverage derived from it; no reads)' : 'StringTie assembly')
        if( chainFiles ) evidence << 'full-length transcripts (chains)'
      }
      def runInfo = [
          version: version, mode: MODE, evidence: evidence ? evidence.join(', ') : 'none (ab initio)',
          genefinder: post_prefix, model: gfModel,
          hc: combinedMode ? 'combined' : null,
      ]
      POSTPROCESS(final_gff, inp.genome, genome_chunks, post_prefix, genefinder_versions, MODE,
                  asm_short, asm_long, hints, runInfo)
    }
}
