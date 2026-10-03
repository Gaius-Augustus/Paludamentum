// Shared functions of the Paludamentum pipeline.

// Gene finders that the pipeline can run. Each one has a params block of the same name.
def genefinders() {
    return ['tiberius', 'vipsania']
}

// A params value as a list: YAML/CLI give a single string or a list.
def asList(v) {
    if( v == null ) return []
    return (v instanceof List) ? v : [v]
}

// Interpret YAML/CLI values such as true, 'true', 1, 'yes' as boolean.
def truthy(v) {
    if( v instanceof Boolean ) return v
    return v?.toString()?.trim()?.toLowerCase() in ['true', '1', 'yes', 'y', 'on']
}

// Name of the selected gene finder: params.genefinder, else the block with run=true.
// Falls back to 'tiberius', which keeps evidence-only runs (no gene finder) working.
def resolveGenefinder(p) {
    if( p.genefinder ) {
        def chosen = p.genefinder.toString().trim().toLowerCase()
        if( !(chosen in genefinders()) )
            error "Unknown params.genefinder '${p.genefinder}'. Supported: ${genefinders().join(', ')}"
        return chosen
    }
    def enabled = genefinders().findAll { name -> truthy(p[name]?.run) }
    if( enabled.size() > 1 )
        error "More than one gene finder has run=true (${enabled.join(', ')}). Set params.genefinder or disable one."
    return enabled ? enabled[0] : 'tiberius'
}

// True if the selected gene finder shall run (<genefinder>.run).
def genefinderEnabled(p) {
    return truthy(p[resolveGenefinder(p)]?.run)
}

// Pipeline modes; 'tiberius' and 'vipsania' are accepted as historic names of 'abinitio'.
def modes() {
    return ['abinitio', 'proteins', 'rnaseq', 'isoseq', 'mixed']
}

// Mode of a run from its inputs. Transcript evidence (short reads, BAM, Iso-Seq)
// needs protein evidence, because the high-confidence genes need both; the
// caller reports that as an error.
def inferMode(boolean hasPaired, boolean hasSingle, boolean hasIso, boolean hasBAM, boolean hasProteins) {
  def hasShortReads = hasPaired || hasSingle || hasBAM
  if( !hasProteins ) {
    if( hasShortReads || hasIso )
      error "RNA-Seq and Iso-Seq evidence need protein evidence (params.proteins or params.odb12Partitions): the high-confidence genes are selected by protein homology."
    return 'abinitio'
  }
  if( hasShortReads && hasIso ) return 'mixed'
  if( hasIso )                  return 'isoseq'
  if( hasShortReads )           return 'rnaseq'
  return 'proteins'
}

// params.mode normalised to one of modes(), or an error for anything else.
def normalizeMode(value) {
  def mode = value.toString().trim().toLowerCase()
  if( mode in ['tiberius', 'vipsania'] ) mode = 'abinitio'
  if( !(mode in modes()) )
    error "Unknown params.mode '${value}'. Supported: ${modes().join(', ')}."
  return mode
}

// True if the selected gene finder model is one that the Drusilla flow serves:
// Tiberius vertebrates and mammalia*, Vipsania Vertebrata (etb1go6q). Drusilla's
// released model and the LightGBM filter are trained on vertebrates.
def drusillaModelEligible(p) {
    def gf = resolveGenefinder(p)
    if( gf == 'tiberius' ) {
        def cfg = p.tiberius?.model_cfg
        if( !cfg ) return false
        def name = new File(cfg.toString()).name.replaceFirst(/\.ya?ml$/, '')
        return name == 'vertebrates' || name.startsWith('mammalia')
    }
    if( gf == 'vipsania' ) {
        return p.vipsania?.model?.toString()?.trim()?.toLowerCase() in ['etb1go6q', 'vertebrata']
    }
    return false
}

// High-confidence gene step of an evidence run: [method: 'drusilla' | 'transdecoder', note: text or null].
// params.drusilla.run 'auto' selects Drusilla for runs with transcripts and a
// vertebrate gene finder model; true and false force the choice.
def hcMethod(p, String mode) {
    def run = p.drusilla?.run
    def auto = run == null || run.toString().trim().toLowerCase() == 'auto'
    def hasTranscripts = mode in ['rnaseq', 'isoseq', 'mixed']
    if( !auto && !truthy(run) ) return [method: 'transdecoder', note: null]
    if( auto ) {
        if( !hasTranscripts || !genefinderEnabled(p) )
            return [method: 'transdecoder', note: null]
        if( resolveGenefinder(p) == 'tiberius' && p.tiberius?.result && !p.tiberius?.model_cfg )
            return [method: 'transdecoder',
                    note: "params.tiberius.result is set without params.tiberius.model_cfg, so the pipeline cannot tell whether the prediction comes from a vertebrate model; using the TransDecoder high-confidence genes. Set params.tiberius.model_cfg, or params.drusilla.run = true, for the Drusilla flow."]
        if( !drusillaModelEligible(p) )
            return [method: 'transdecoder', note: null]
        if( !p.drusilla?.lgb_model )
            return [method: 'transdecoder',
                    note: "The gene finder model is a vertebrate model, but params.drusilla.lgb_model is not set; using the TransDecoder high-confidence genes."]
        return [method: 'drusilla', note: null]
    }
    if( !hasTranscripts )     error "params.drusilla.run = true needs transcripts (mode rnaseq, isoseq or mixed), the mode is '${mode}'."
    if( !genefinderEnabled(p) ) error "params.drusilla.run = true needs a gene finder (tiberius.run or vipsania.run)."
    if( !p.drusilla?.lgb_model ) error "params.drusilla.run = true needs params.drusilla.lgb_model (LightGBM model of the gene finder filter)."
    return [method: 'drusilla', note: null]
}
