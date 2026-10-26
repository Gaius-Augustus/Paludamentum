// Shared functions of the Paludamentum pipeline.

// Gene finders that the pipeline can run. Each one has a params block of the same name.
// Note: in evidence modes both gene finders always run (unified HC flow in HC_COMBINED);
// resolveGenefinder / genefinderEnabled are only consulted in abinitio and proteins modes.
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

// Name of the selected gene finder for ab initio / proteins mode: params.genefinder,
// else the block with run=true. Falls back to 'tiberius'.
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

// True if the selected gene finder shall run (<genefinder>.run). In ab initio
// and proteins modes only; evidence modes always run both.
def genefinderEnabled(p) {
    return truthy(p[resolveGenefinder(p)]?.run)
}

// Pipeline modes; 'tiberius' and 'vipsania' are accepted as historic names of 'abinitio'.
def modes() {
    return ['abinitio', 'proteins', 'rnaseq', 'isoseq', 'mixed']
}

// Mode of a run from its inputs.
def inferMode(boolean hasPaired, boolean hasSingle, boolean hasIso, boolean hasBAM, boolean hasProteins) {
  def hasShortReads = hasPaired || hasSingle || hasBAM
  if( hasShortReads && hasIso ) return 'mixed'
  if( hasIso )                  return 'isoseq'
  if( hasShortReads )           return 'rnaseq'
  return hasProteins ? 'proteins' : 'abinitio'
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
    } catch( Exception _e ) {
        return null
    }
}

// target_species of a Tiberius model configuration, null if it has none.
def tiberiusTargetSpecies(cfg) {
    return tiberiusModelValue(cfg, 'target_species')
}

// Clade table (params.hc_table, conf/hc_genes.yaml). Only used to look up
// the clade of the gene finder model (insect flag + recommended warning).
def hcTable(p) {
    def path = p.hc_table?.toString()?.trim() ?: "${projectDir}/conf/hc_genes.yaml"
    def f = file(path)
    if( !f.exists() ) error "params.hc_table: ${path} not found."
    def data = new org.yaml.snakeyaml.Yaml().load(f.text)
    if( !(data instanceof Map) ) error "params.hc_table: ${path} is not a YAML mapping."
    return data
}

// Clade of the gene finder model in the HC table: [name: the clade (the table's
// spelling when listed), entry: its table entry, or null when not listed].
// Tiberius: target_species of the model configuration. Vipsania: the model,
// a clade name or one of the clade's vipsania_models. Case is ignored.
// In evidence modes both gene finders run; the clade is read from whichever
// model is set (Tiberius first, then Vipsania).
def hcClade(p) {
    def tibCfg = p.tiberius?.model_cfg
    def vipMod = p.vipsania?.model?.toString()?.trim()
    def candidates = []
    if( tibCfg ) {
        def ts = tiberiusTargetSpecies(tibCfg)
        if( ts ) candidates << ts
    }
    if( vipMod ) candidates << vipMod
    for( name in candidates ) {
        def key = name.toLowerCase()
        def clades = hcTable(p).clades
        def hit = (clades instanceof Map) ? clades.find { k, v ->
            k.toString().toLowerCase() == key ||
                (v instanceof Map && (v.vipsania_models instanceof List) &&
                 v.vipsania_models.any { m -> m.toString().toLowerCase() == key })
        } : null
        if( hit ) return [name: hit.key.toString(), entry: (hit.value instanceof Map) ? hit.value : [:]]
    }
    return [name: candidates ? candidates[0] : null, entry: null]
}

// Drusilla model (key 'model'): params.drusilla.model, else the clade's
// drusilla_model in the HC table. Needed in the unified HC flow to annotate
// ORFs with the right clade's Drusilla weights.
def drusillaSetting(p, String key) {
    if( key != 'model' )
        error "drusillaSetting: only 'model' is supported, got '${key}'."
    def d = p.drusilla ?: [:]
    if( d.model ) return d.model
    def entry = hcClade(p).entry ?: [:]
    return entry.drusilla_model ?: null
}

// The most accurate flow measured for the clade of the gene finder model, as
// the key 'recommended' of its entry in the HC table (conf/hc_genes.yaml):
// genefinder, finetune (Vipsania), transcripts, proteins (true: given, false:
// not given), and note (the measured accuracy). Returns null if the clade
// has no recommendation or the run follows it, else the warning.
def cladeAdvice(p, String mode, boolean hasProteins, Map hc) {
    def clade = hcClade(p)
    def rec = clade.entry?.recommended
    if( !(rec instanceof Map) ) return null
    def transcripts = mode in ['rnaseq', 'isoseq', 'mixed']
    def diff = []
    if( rec.finetune != null && !p.vipsania?.result &&
        truthy(rec.finetune) != truthy(p.vipsania?.finetune) )
        diff << (truthy(rec.finetune) ? 'Vipsania is not fine-tuned on the genome (vipsania.finetune: true)' :
                                        'Vipsania is fine-tuned on the genome')
    if( rec.transcripts != null && truthy(rec.transcripts) != transcripts )
        diff << (transcripts ? 'there is RNA-Seq or Iso-Seq evidence' : 'there is no RNA-Seq or Iso-Seq evidence')
    if( rec.proteins != null && truthy(rec.proteins) != hasProteins )
        diff << (hasProteins ? 'protein evidence is given (proteins or odb12Partitions)' : 'there is no protein evidence')
    if( !diff ) return null
    def note = rec.note ? " ${rec.note.toString().trim()}" : ''
    return "This run of clade ${clade.name} differs from the flow recommended for it: ${diff.join('; ')}.${note}"
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

// True if the header of a BAM file says SO:coordinate (@HD, its first line).
// No @HD line, another sort order, SAM/CRAM or an unreadable file: false.
// BAM is BGZF, a series of gzip members: magic BAM\1, int32 l_text, header text.
def bamCoordinateSorted(f) {
    try {
        new java.util.zip.GZIPInputStream(f.newInputStream()).withCloseable { s ->
            def head = s.readNBytes(8)
            if( head.length < 8 || new String(head, 0, 4, 'ISO-8859-1') != 'BAM' )
                return false
            def lText = java.nio.ByteBuffer.wrap(head, 4, 4).order(java.nio.ByteOrder.LITTLE_ENDIAN).getInt()
            def text = s.readNBytes(Math.min(Math.max(lText, 0), 65536))
            def first = new String(text, 'ISO-8859-1').tokenize('\n')[0] ?: ''
            return first.startsWith('@HD\t') && first.tokenize('\t').contains('SO:coordinate')
        }
    } catch( Exception _e ) {
        return false
    }
}
