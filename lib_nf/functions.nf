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

// Mode of a run from its inputs. Transcript evidence (short reads, BAM, Iso-Seq,
// pyVARUS directories) needs protein evidence, because the high-confidence
// genes need both; the caller reports that as an error. The caller counts
// rnaseq_varus in hasBAM and isoseq_varus in hasIso.
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

// Value of key in a Tiberius model configuration or a Drusilla model manifest
// (YAML), null if the file cannot be read or has no such key.
def tiberiusModelValue(cfg, String key) {
    try {
        def data = new org.yaml.snakeyaml.Yaml().load(file(cfg.toString()).text)
        return (data instanceof Map) ? data[key]?.toString()?.trim() ?: null : null
    } catch( Exception e ) {
        return null
    }
}

// target_species of a Tiberius model configuration, null if it has none.
def tiberiusTargetSpecies(cfg) {
    return tiberiusModelValue(cfg, 'target_species')
}

// True if the selected gene finder model is one that the Drusilla flow serves:
// Tiberius models with target_species Vertebrata or Mammalia (vertebrates,
// mammalia*), Vipsania Vertebrata (etb1go6q). Drusilla's released model and
// the LightGBM filter are trained on vertebrates.
def drusillaModelEligible(p) {
    def gf = resolveGenefinder(p)
    if( gf == 'tiberius' ) {
        def cfg = p.tiberius?.model_cfg
        return cfg && tiberiusTargetSpecies(cfg)?.toLowerCase() in ['vertebrata', 'mammalia']
    }
    if( gf == 'vipsania' ) {
        return p.vipsania?.model?.toString()?.trim()?.toLowerCase() in ['etb1go6q', 'vertebrata']
    }
    return false
}

// True if the Drusilla flow runs the hint rescue (params.drusilla.rescue).
def rescueEnabled(p) {
    return p.drusilla?.rescue == null || truthy(p.drusilla.rescue)
}

// Tiberius model of the hint rescue: params.drusilla.rescue_model_cfg, else the
// model of a Tiberius run, else vertebrates.
// [file: model configuration to stage, or null; name: model in model_cfg/ of the
// rescue image, or null; weights: true if params.tiberius.model_dir holds its
// weights; note: text or null]. A value without '/' and without a .yaml/.yml
// suffix is a name.
def rescueModel(p) {
    def own = p.drusilla?.rescue_model_cfg?.toString()?.trim()
    if( own ) {
        def isName = !own.contains('/') && !(own ==~ /.*\.ya?ml$/)
        return [file: isName ? null : own, name: isName ? own : null, weights: false, note: null]
    }
    def tiberiusRun = resolveGenefinder(p) == 'tiberius'
    if( tiberiusRun && p.tiberius?.model_cfg )
        return [file: p.tiberius.model_cfg.toString(), name: null, weights: p.tiberius?.model_dir as boolean, note: null]
    def note = tiberiusRun ?
        "params.tiberius.result is set without params.tiberius.model_cfg, so the hint rescue uses the Tiberius model vertebrates. Set params.drusilla.rescue_model_cfg for another model." : null
    return [file: null, name: 'vertebrates', weights: false, note: note]
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
        if( resolveGenefinder(p) == 'tiberius' && p.tiberius?.model_cfg && !tiberiusTargetSpecies(p.tiberius.model_cfg) )
            return [method: 'transdecoder',
                    note: "params.tiberius.model_cfg (${p.tiberius.model_cfg}) has no target_species, so the pipeline cannot tell whether it is a vertebrate model; using the TransDecoder high-confidence genes. Set params.drusilla.run = true for the Drusilla flow."]
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
    return [method: 'drusilla', note: rescueEnabled(p) ? rescueModel(p).note : null]
}

// Header keys (#key=value lines) of a pyVARUS manifest; '# ' lines are comments.
def varusHeader(f) {
    def header = [:]
    def lines = f.readLines().findAll { line -> line.trim() && !line.startsWith('# ') }
    lines.takeWhile { line -> line.startsWith('#') }.each { line ->
        def i = line.indexOf('=')
        if( i > 1 ) header[line.substring(1, i)] = line.substring(i + 1)
    }
    return header
}

// pyVARUS output directories of params key (rnaseq_varus, isoseq_varus,
// mixed_varus), checked at workflow construction: one
// [id, dir, manifest, stringtie.gtf, hints.gff or []] per directory.
// pyVARUS writes stringtie.gtf and hints.gff as Paludamentum makes them from a
// BAM; VARUS.manifest.tsv (varus run) or VARUS.assembly.tsv (varus assemble)
// holds the mode and the genome MD5, which VARUS_INPUT checks.
// mode: 'shortreads', 'longreads' or 'mixed'.
def varusInputs(value, String key, String mode) {
    def takes = [rnaseq_varus: 'short-read runs', isoseq_varus: '--longreads runs',
                 mixed_varus: 'the output of varus assemble --short --long']
    def option = mode == 'longreads' ? '--long' : '--short'
    return asList(value).findAll { v -> v }.collect { v ->
        def dir = file(v.toString())
        if( !dir.exists() ) error "${key}: ${dir} does not exist."
        if( !dir.isDirectory() ) error "${key}: ${dir} is not a directory; pass the output directory of pyVARUS."
        def gtf   = dir.resolve('stringtie.gtf')
        def hints = dir.resolve('hints.gff')
        if( mode == 'mixed' && !gtf.exists() )
            error "${key}: no stringtie.gtf in ${dir}. Make it with `varus assemble GENOME --short A/VARUS.bam --long B/VARUS.bam --outdir ${dir}`."
        if( mode != 'mixed' && !(gtf.exists() && hints.exists()) )
            error "${key}: no stringtie.gtf or hints.gff in ${dir} (pyVARUS older than its assembly step, or a failed run). " +
                  "Make them with `varus assemble GENOME ${option} ${dir}/VARUS.bam --outdir ${dir}`."
        def names = ['VARUS.assembly.tsv', 'VARUS.manifest.tsv']
        def manifest = names.collect { n -> dir.resolve(n) }.find { f -> f.exists() }
        if( !manifest ) error "${key}: no ${names.join(' or ')} in ${dir}; it holds the genome MD5 and the mode of the pyVARUS run."
        def have = varusHeader(manifest).mode ?: 'shortreads'
        if( have != mode )
            error "${key}: ${dir} is a pyVARUS run in mode '${have}' (${manifest.name}), expected '${mode}'. " +
                  "rnaseq_varus takes ${takes.rnaseq_varus}, isoseq_varus ${takes.isoseq_varus}, mixed_varus ${takes.mixed_varus}."
        return [dir.name, dir.toString(), manifest, gtf, mode == 'mixed' ? [] : hints]
    }
}
