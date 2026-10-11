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
// pyVARUS directories) without protein evidence gives the transcript mode as
// well: main.nf stops such a run unless the Drusilla flow selects the
// high-confidence genes, which needs no proteins (hcMethod). The caller counts
// rnaseq_varus in hasBAM and isoseq_varus in hasIso.
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

// HC gene step per clade (params.hc_table, conf/hc_genes.yaml): orf_finder and
// clades (clade name -> hc, drusilla_model, lgb_model, lgb_model_sha256,
// orf_finder, vipsania_models).
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
def hcClade(p) {
    def gf = resolveGenefinder(p)
    def name = null
    if( gf == 'tiberius' && p.tiberius?.model_cfg ) name = tiberiusTargetSpecies(p.tiberius.model_cfg)
    else if( gf == 'vipsania' ) name = p.vipsania?.model?.toString()?.trim() ?: null
    if( !name ) return [name: null, entry: null]
    def key = name.toLowerCase()
    def clades = hcTable(p).clades
    def hit = (clades instanceof Map) ? clades.find { k, v ->
        k.toString().toLowerCase() == key ||
            (gf == 'vipsania' && v instanceof Map && (v.vipsania_models instanceof List) &&
             v.vipsania_models.any { m -> m.toString().toLowerCase() == key })
    } : null
    if( !hit ) return [name: name, entry: null]
    return [name: hit.key.toString(), entry: (hit.value instanceof Map) ? hit.value : [:]]
}

// ORF finder of the TransDecoder flow: params.transdecoder, else orf_finder of
// the clade in the HC table, else the table's orf_finder, else td2.
def orfFinder(p) {
    def td = p.transdecoder?.toString()?.trim()?.toLowerCase()
    if( !td ) {
        def own = hcClade(p).entry?.orf_finder
        td = (own ?: hcTable(p).orf_finder ?: 'td2').toString().trim().toLowerCase()
    }
    if( !(td in ['td1', 'td2']) )
        error "The ORF finder must be 'td1' or 'td2' (params.transdecoder, or orf_finder in the HC table), got '${td}'."
    return td
}

// Drusilla model (key 'model'), LightGBM model ('lgb_model') and its sha256
// ('lgb_model_sha256'): params.drusilla, else the clade's entry in the HC
// table. With params.drusilla.run = true, a clade without Drusilla models of
// its own takes those of the table's drusilla_forced clade. A
// params.drusilla.lgb_model is checked only against
// params.drusilla.lgb_model_sha256 (null: not checked).
def drusillaSetting(p, String key) {
    def d = p.drusilla ?: [:]
    def entry = hcClade(p).entry ?: [:]
    def run = d.run?.toString()?.trim()?.toLowerCase()
    if( run && run != 'auto' && truthy(d.run) && !entry.drusilla_model && !entry.lgb_model ) {
        def table = hcTable(p)
        def forced = table.drusilla_forced?.toString()?.toLowerCase()
        def fallback = (forced && table.clades instanceof Map) ?
            table.clades.find { k, _v -> k.toString().toLowerCase() == forced }?.value : null
        if( fallback instanceof Map ) entry = fallback
    }
    if( key == 'model' )            return d.model ?: entry.drusilla_model ?: null
    if( key == 'lgb_model' )        return d.lgb_model ?: entry.lgb_model ?: null
    if( key == 'lgb_model_sha256' ) return d.lgb_model_sha256 ?: (d.lgb_model ? null : entry.lgb_model_sha256 ?: null)
    error "drusillaSetting: unknown key '${key}'."
}

// True if the Drusilla flow runs the hint rescue (params.drusilla.rescue).
// The rescue predicts the partial genes of the LightGBM filter again, so it
// runs only where the filter runs (drusillaLgb).
def rescueEnabled(p) {
    return p.drusilla?.rescue == null || truthy(p.drusilla.rescue)
}

// True if the Drusilla flow runs the LightGBM filter of the gene finder genes:
// the clade (or params.drusilla) has a LightGBM model, and there is protein
// evidence, from which the filter computes its features. Without it the
// Drusilla ORFs are merged with all gene finder genes.
def drusillaLgb(p, boolean hasProteins) {
    return hasProteins && drusillaSetting(p, 'lgb_model') as boolean
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

// TransDecoder result of hcMethod.
def hcTransdecoder(clade, note) {
    return [method: 'transdecoder', clade: clade, note: note]
}

// High-confidence gene step of an evidence run: [method: 'drusilla' |
// 'transdecoder', clade: clade of the gene finder model or null, note: text or
// null]. params.drusilla.run 'auto' follows hc of the clade in the HC table
// (conf/hc_genes.yaml) for runs with transcripts and a gene finder; true and
// false force the choice. The Drusilla flow needs a Drusilla model only; a
// clade without a LightGBM model runs it without the filter (drusillaLgb).
def hcMethod(p, String mode) {
    def run = p.drusilla?.run
    def auto = run == null || run.toString().trim().toLowerCase() == 'auto'
    def hasTranscripts = mode in ['rnaseq', 'isoseq', 'mixed']
    def clade = hcClade(p)
    def hc = (clade.entry?.hc ?: 'transdecoder').toString().trim().toLowerCase()
    if( !(hc in ['drusilla', 'transdecoder']) )
        error "hc of clade ${clade.name} in the HC table must be 'drusilla' or 'transdecoder', got '${hc}'."
    def missing = []
    if( !p.drusilla?.weights && !drusillaSetting(p, 'model') ) missing << 'Drusilla model (params.drusilla.model)'
    if( !auto && !truthy(run) ) return hcTransdecoder(clade.name, null)
    if( auto ) {
        if( !hasTranscripts || !genefinderEnabled(p) )
            return hcTransdecoder(clade.name, null)
        if( resolveGenefinder(p) == 'tiberius' && p.tiberius?.result && !p.tiberius?.model_cfg )
            return hcTransdecoder(clade.name, "params.tiberius.result is set without params.tiberius.model_cfg, so the pipeline cannot tell which clade the prediction comes from; using the TransDecoder high-confidence genes. Set params.tiberius.model_cfg, or params.drusilla.run = true, for the Drusilla flow.")
        if( resolveGenefinder(p) == 'tiberius' && p.tiberius?.model_cfg && !clade.name )
            return hcTransdecoder(clade.name, "params.tiberius.model_cfg (${p.tiberius.model_cfg}) has no target_species, so the pipeline cannot tell its clade; using the TransDecoder high-confidence genes. Set params.drusilla.run = true for the Drusilla flow.")
        if( hc != 'drusilla' )
            return hcTransdecoder(clade.name, null)
        if( missing )
            return hcTransdecoder(clade.name, "The HC table selects the Drusilla flow for ${clade.name}, but there is no ${missing.join(' and no ')}, neither for the clade in the HC table nor in the params; using the TransDecoder high-confidence genes.")
        return [method: 'drusilla', clade: clade.name, note: null]
    }
    if( !hasTranscripts )     error "params.drusilla.run = true needs transcripts (mode rnaseq, isoseq or mixed), the mode is '${mode}'."
    if( !genefinderEnabled(p) ) error "params.drusilla.run = true needs a gene finder (tiberius.run or vipsania.run)."
    if( missing ) error "params.drusilla.run = true needs a ${missing.join(' and a ')}; the HC table has none for clade ${clade.name ?: 'unknown'}, nor a drusilla_forced clade with one."
    return [method: 'drusilla', clade: clade.name, note: rescueEnabled(p) ? rescueModel(p).note : null]
}

// The most accurate flow measured for the clade of the gene finder model, as
// the key 'recommended' of its entry in the HC table (conf/hc_genes.yaml):
// genefinder, finetune (Vipsania), transcripts, drusilla (the Drusilla flow
// selects the HC genes), proteins (true: given, false: not given), and note
// (the measured accuracy). Keys that are not set are not checked. Returns
// null if the clade has no recommendation or the run follows it, else the
// warning, which names every difference.
def cladeAdvice(p, String mode, boolean hasProteins, Map hc) {
    def clade = hcClade(p)
    def rec = clade.entry?.recommended
    if( !(rec instanceof Map) ) return null
    def tool = resolveGenefinder(p)
    def running = genefinderEnabled(p)
    def transcripts = mode in ['rnaseq', 'isoseq', 'mixed']
    def diff = []
    def wantTool = rec.genefinder?.toString()?.trim()?.toLowerCase()
    if( wantTool && (!running || tool != wantTool) )
        diff << (running ? "the gene finder is ${tool}, not ${wantTool}" : "no gene finder runs (${wantTool}.run)")
    if( running && tool == 'vipsania' && rec.finetune != null && !p.vipsania?.result &&
        truthy(rec.finetune) != truthy(p.vipsania?.finetune) )
        diff << (truthy(rec.finetune) ? 'Vipsania is not fine-tuned on the genome (vipsania.finetune: true)' :
                                        'Vipsania is fine-tuned on the genome')
    if( rec.transcripts != null && truthy(rec.transcripts) != transcripts )
        diff << (transcripts ? 'there is RNA-Seq or Iso-Seq evidence' : 'there is no RNA-Seq or Iso-Seq evidence')
    if( transcripts && rec.drusilla != null && truthy(rec.drusilla) != (hc.method == 'drusilla') )
        diff << (hc.method == 'drusilla' ? 'the Drusilla flow selects the high-confidence genes' :
                                           "the high-confidence genes are not those of the Drusilla flow (HC genes: ${hc.method})")
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
            if( head.length < 8 || new String(head, 0, 4, 'ISO-8859-1') != 'BAM\u0001' )
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
