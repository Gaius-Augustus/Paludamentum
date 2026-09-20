// Shared functions of the Paludamentum pipeline.

// Gene finders that the pipeline can run. Each one has a params block of the same name.
def genefinders() {
    return ['tiberius', 'vipsania']
}

// Interpret YAML/CLI values such as true, 'true', 1, 'yes' as boolean.
def truthy(v) {
    if( v instanceof Boolean ) return v
    return v?.toString()?.trim()?.toLowerCase() in ['true', '1', 'yes', 'y']
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

def inferMode(boolean hasPaired, boolean hasSingle, boolean hasIso, boolean hasBAM, boolean hasProteins) {
  if ((hasPaired || hasSingle || hasBAM) && hasIso && hasProteins) return 'mixed'
  if (hasIso && hasProteins) return 'isoseq'
  if (hasPaired || hasSingle || hasBAM && hasProteins) return 'rnaseq'
  if (hasProteins && hasProteins) return 'proteins'
  return 'tiberius'
}
