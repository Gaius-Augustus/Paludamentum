nextflow.enable.dsl=2

include { drusillaSetting; truthy } from '../lib_nf/functions.nf'

// Processes of the Drusilla flow: the StringTie pre-filter, model download,
// the ORF annotation and the subsequence collapse. The LightGBM filter and
// the hint rescue of the previous pipeline are gone (replaced by the unified
// combine_genesets integration in subworkflows/hc_combined.nf).

process FILTER_STRINGTIE {
  label 'container'
  input:
    path gtf

  output:
    path "stringtie.filtered.gtf", emit: gtf

  script:
  def dropUnstranded = params.drusilla.drop_unstranded == null || params.drusilla.drop_unstranded.toString().toLowerCase() in ['true', '1', 'yes']
  """
  filter_stringtie_gtf.py \\
      --in-gtf ${gtf} \\
      --out-gtf stringtie.filtered.gtf \\
      --out-tsv stringtie.decisions.tsv \\
      --min-length ${params.drusilla.min_length} \\
      --min-cov ${params.drusilla.min_cov} \\
      --min-tpm ${params.drusilla.min_tpm} \\
      --long-length ${params.drusilla.long_length} \\
      --min-tpm-long ${params.drusilla.min_tpm_long}${dropUnstranded ? ' \\\n      --drop-unstranded' : ''}
  """

  stub:
  """
  touch stringtie.filtered.gtf
  """
}

// Weights and architecture of a released Drusilla model, downloaded once on the
// submitting host (label 'download'), so the GPU nodes need no internet and the
// shards of DRUSILLA_ANNOTATE do not download into one cache at the same time.
// The archive is checked against the weights_sha256 of the model's manifest.
// Output: drusilla_model/ with model.weights.h5 (Keras 3 loads *.weights.h5 with
// the loader the weights were saved with) and arch.yaml.
process DOWNLOAD_DRUSILLA_MODEL {
  tag "${url.tokenize('/')[-1]}"
  label 'download'
  input:
    val url
    val sha256

  output:
    path "drusilla_model", emit: model

  script:
  def check = sha256 ?
      "echo \"${sha256}  model.archive\" | sha256sum -c --quiet - || { echo \"The sha256 of ${url} is not ${sha256}\" >&2; exit 1; }" :
      "echo \"The manifest of ${url} has no weights_sha256, the archive is not checked\" >&2"
  """
  set -euo pipefail
  curl -fsSL -o model.archive "${url}"
  ${check}
  mkdir -p archive drusilla_model
  tar -xf model.archive -C archive
  rm model.archive
  arch=\$(find archive -name arch.yaml | sort | head -n 1)
  if [ -z "\$arch" ]; then
      echo "The model archive ${url} has no arch.yaml" >&2
      exit 1
  fi
  dir=\$(dirname "\$arch")
  weights=\$(find "\$dir" -maxdepth 1 -name '*.weights.h5' | sort | head -n 1)
  [ -n "\$weights" ] || weights="\$dir/weights.h5"
  if [ ! -f "\$weights" ]; then
      echo "The model archive ${url} has no weights next to arch.yaml" >&2
      exit 1
  fi
  mv "\$arch" drusilla_model/arch.yaml
  mv "\$weights" drusilla_model/model.weights.h5
  rm -rf archive
  """

  stub:
  """
  mkdir -p drusilla_model
  touch drusilla_model/arch.yaml drusilla_model/model.weights.h5
  """
}

// All ORF predictions (the subsequence collapse follows in COLLAPSE_ORFS).
process DRUSILLA_ANNOTATE {
  label 'gpu', 'drusilla', 'bigmem'
  input:
    path gtf
    path genome
    // use_model: model holds drusilla_model/ of DOWNLOAD_DRUSILLA_MODEL. false:
    // params.drusilla.weights, or the Drusilla registry resolves the model name
    // (in params.drusilla.cache_dir, else downloaded in the task)
    val use_model
    path model

  output:
    path "drusilla/orfs.gtf", emit: gtf
    path "drusilla/orfs.partial.gtf", emit: partial
    path "drusilla/orfs.partial5.gtf", emit: partial5

  script:
  def name = drusillaSetting(params, 'model')
  def model_args = params.drusilla.weights ?
      "--weights \"${params.drusilla.weights}\"" + (params.drusilla.config ? " --config \"${params.drusilla.config}\"" : '') :
      use_model ? "--weights drusilla_model/model.weights.h5 --config drusilla_model/arch.yaml" :
      "--model ${name}"
  // The registry resolves the model once before the shards start: it deletes
  // the model directory before it extracts, so shards must not extract at once.
  def registry = (params.drusilla.weights || use_model) ? '' :
      "export DRUSILLA_CACHE_DIR=\"${params.drusilla.cache_dir ?: '\$PWD/drusilla_cache'}\"\n" +
      "  drusilla models download ${name} >&2"
  def extra = ''
  if( params.drusilla.batch_size )                extra += " --batch-size ${params.drusilla.batch_size}"
  if( params.drusilla.min_coding_length != null ) extra += " --min-coding-length ${params.drusilla.min_coding_length}"
  // Transcript handling of Drusilla >= 0.2 (docs/drusilla_flow.md): both
  // strands of unstranded transcripts, ORF completion into genomic flanks,
  // several ORFs per transcript. All off unless set.
  if( truthy(params.drusilla.unstranded_both_strands) ) extra += " --unstranded-both-strands"
  def flank  = (params.drusilla.flank_bp ?: 0) as Integer
  def flank3 = (params.drusilla.flank_3p_bp != null ? params.drusilla.flank_3p_bp : flank) as Integer
  if( flank > 0 || flank3 > 0 ) {
      extra += " --flank-bp ${flank} --flank-3p-bp ${flank3} --flank-keep"
      if( params.drusilla.flank_min_tx_frac != null ) extra += " --flank-min-tx-frac ${params.drusilla.flank_min_tx_frac}"
  }
  if( params.drusilla.orfs_per_tx != null )           extra += " --orfs-per-tx ${params.drusilla.orfs_per_tx}"
  // shards > 1: the transcripts are split by gene into this many parts, each
  // annotated by its own process. On CPUs one process uses only 1-2 cores.
  def shards = (params.drusilla.shards ?: 1) as Integer
  def threads = Math.max(1, (task.cpus as Integer).intdiv(shards))
  """
  ${registry}
  mkdir -p drusilla
  annotate() {
      drusilla annotate \\
          --stringtie-gtf "\$1" \\
          --genome ${genome} \\
          ${model_args} \\
          --out-dir "\$2" \\
          --threads ${threads} \\
          --no-subseq-collapse \\
          --lorf-class \\
          --partial-out "\$2/orfs.partial.gtf" \\
          --partial5-out "\$2/orfs.partial5.gtf"${extra}
      touch "\$2/orfs.partial.gtf" "\$2/orfs.partial5.gtf"
  }
  if [ ${shards} -le 1 ]; then
      annotate ${gtf} drusilla
  else
      awk -F'\t' -v n=${shards} '/^#/ {next} {
              match(\$9, /gene_id "[^"]+"/); g = substr(\$9, RSTART, RLENGTH)
              if (!(g in shard)) shard[g] = k++ % n
              print > ("shard_" shard[g] ".gtf") }' ${gtf}
      pids=""
      for f in shard_*.gtf; do
          i=\${f#shard_}; i=\${i%.gtf}
          TF_NUM_INTRAOP_THREADS=${threads} OMP_NUM_THREADS=${threads} \\
              annotate \$f drusilla_\$i > drusilla_\$i.log 2>&1 &
          pids="\$pids \$!"
      done
      for p in \$pids; do wait \$p || { cat drusilla_*.log >&2; exit 1; }; done
      for o in orfs.gtf orfs.partial.gtf orfs.partial5.gtf; do
          cat drusilla_*/\$o > drusilla/\$o
      done
  fi
  """

  stub:
  """
  mkdir -p drusilla
  touch drusilla/orfs.gtf drusilla/orfs.partial.gtf drusilla/orfs.partial5.gtf
  """
}

// Subsequence collapse: isoforms whose CDS is a subsequence of another one are
// dropped. No codon is changed; the unified HC flow (combine_genesets.py)
// handles start/stop refinements through its tiering rules.
process COLLAPSE_ORFS {
  label 'drusilla'
  publishDir "${params.outdir}/intermediate", pattern: "drusilla_orfs.gtf", mode: 'copy', overwrite: true
  input:
    path orfs, stageAs: 'orfs.raw.gtf'

  output:
    path "drusilla_orfs.gtf", emit: gtf

  script:
  """
  drusilla filter-subseq \\
      --orfs-gtf orfs.raw.gtf \\
      --out-gtf drusilla_orfs.gtf \\
      --report-tsv dropped_subsequences.tsv
  """

  stub:
  """
  touch drusilla_orfs.gtf
  """
}
