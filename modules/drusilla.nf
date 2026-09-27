nextflow.enable.dsl=2

// Processes of the Drusilla flow: the StringTie pre-filter, the Drusilla ORF
// annotation, the stop and start codon fix of the ORFs, the LightGBM filter of
// the gene finder predictions, and the hint rescue of partial gene finder genes.

process FILTER_STRINGTIE {
  label 'container'
  input:
    path gtf

  output:
    path "stringtie.filtered.gtf", emit: gtf

  script:
  def d = params.drusilla
  """
  filter_stringtie_gtf.py \\
      --in-gtf ${gtf} \\
      --out-gtf stringtie.filtered.gtf \\
      --out-tsv stringtie.decisions.tsv \\
      --min-length ${d.min_length} \\
      --min-cov ${d.min_cov} \\
      --min-tpm ${d.min_tpm} \\
      --long-length ${d.long_length} \\
      --min-tpm-long ${d.min_tpm_long}
  """

  stub:
  """
  touch stringtie.filtered.gtf
  """
}

// All ORF predictions (the subsequence collapse follows in FIX_ORFS), with LORF
// classes, and the 3' and 5' truncated ORFs that the stop and start fix can recover.
process DRUSILLA_ANNOTATE {
  label 'gpu', 'drusilla'
  input:
    path gtf
    path genome

  output:
    path "drusilla/orfs.gtf", emit: gtf
    path "drusilla/orfs.partial.gtf", emit: partial
    path "drusilla/orfs.partial5.gtf", emit: partial5

  script:
  def d = params.drusilla
  def model = d.weights ? "--weights ${d.weights}" + (d.config ? " --config ${d.config}" : '') \
                        : "--model ${d.model ?: 'vertebrates'}"
  def extra = ''
  if( d.batch_size )                extra += " --batch-size ${d.batch_size}"
  if( d.min_coding_length != null ) extra += " --min-coding-length ${d.min_coding_length}"
  def cache = d.cache_dir ?: '\$PWD/drusilla_cache'
  // shards > 1: the transcripts are split by gene into this many parts, each
  // annotated by its own process. On CPUs one process uses only 1-2 cores.
  def shards = (d.shards ?: 1) as Integer
  def threads = Math.max(1, (task.cpus as Integer).intdiv(shards))
  """
  export DRUSILLA_CACHE_DIR=${cache}
  mkdir -p drusilla
  annotate() {
      drusilla annotate \\
          --stringtie-gtf \$1 \\
          --genome ${genome} \\
          ${model} \\
          --out-dir \$2 \\
          --threads ${threads} \\
          --no-subseq-collapse \\
          --lorf-class \\
          --partial-out \$2/orfs.partial.gtf \\
          --partial5-out \$2/orfs.partial5.gtf${extra}
      touch \$2/orfs.partial.gtf \$2/orfs.partial5.gtf
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

// Stop codon fix (fix_stop): complete ORFs with an early stop and truncated
// ORFs are extended to a stop codon supported by a protein alignment.
// Start codon fix (fix_start): complete ORFs of the LORF classes LORF_NOUPSTOP
// and upLORF are extended to a start codon hint of miniprothint.
// Then isoforms whose CDS is a subsequence of another one are collapsed.
process FIX_ORFS {
  label 'drusilla'
  publishDir "${params.outdir}/intermediate", pattern: "drusilla_orfs.gtf", mode: 'copy'
  input:
    path orfs,     stageAs: 'orfs.raw.gtf'
    path partial,  stageAs: 'orfs.partial.gtf'
    path partial5, stageAs: 'orfs.partial5.gtf'
    path miniprot_gff
    path hints_gff
    path genome

  output:
    path "drusilla_orfs.gtf", emit: gtf

  script:
  def d = params.drusilla
  def fixStop  = d.fix_stop == null || d.fix_stop.toString().toLowerCase() in ['true', '1', 'yes']
  def fixStart = fixStop && (d.fix_start == null || d.fix_start.toString().toLowerCase() in ['true', '1', 'yes'])
  def fix = !fixStop ? "cp orfs.raw.gtf orfs.fixed.gtf" : """\
  fix_stop_by_miniprot.py \\
      --orfs orfs.raw.gtf \\
      --partial orfs.partial.gtf \\
      --partial5 orfs.partial5.gtf \\
      --miniprot ${miniprot_gff} \\
      --hints ${hints_gff} \\
      --genome ${genome} \\
      --out orfs.fixed.gtf""" + (fixStart ? """ \\
      --fix-starts \\
      --fix-starts-classes LORF_NOUPSTOP upLORF \\
      --max-start-scan 30""" : '')
  """
  ${fix}
  drusilla filter-subseq \\
      --orfs-gtf orfs.fixed.gtf \\
      --out-gtf drusilla_orfs.gtf \\
      --report-tsv dropped_subsequences.tsv
  """

  stub:
  """
  touch drusilla_orfs.gtf
  """
}

// Scores each gene finder transcript as wrong, partial or correct from protein
// alignment and hint features. Transcripts with P(partial) + P(correct) >=
// lgb_threshold get their most likely class; those in lgb_keep are kept, those
// of class partial are the candidates of the hint rescue.
process GENEFINDER_LGB_FILTER {
  label 'drusilla'
  publishDir "${params.outdir}/intermediate", pattern: "*_lgb_{filtered.gtf,scores.tsv}", mode: 'copy'
  input:
    val prefix
    path ab_initio
    path miniprot_gff
    path hints_gff
    path genome
    path lgb_model

  output:
    path "${prefix}_lgb_filtered.gtf", emit: gtf
    path "${prefix}_lgb_partial.gtf", emit: partial
    path "${prefix}_lgb_scores.tsv", emit: scores

  script:
  def d = params.drusilla
  def keep = (d.lgb_keep ?: 'correct').toString().split(/[,\s]+/).findAll { k -> k }.join('|')
  """
  gff_to_cds_gtf.py ${ab_initio} > ab_initio.gtf
  compute_orf_features.py \\
      --orfs-gtf ab_initio.gtf \\
      --miniprot-gff ${miniprot_gff} \\
      --hints-gff ${hints_gff} \\
      --genome ${genome} \\
      --out features.tsv
  apply_lgb_model_gtf.py \\
      --model ${lgb_model} \\
      --features features.tsv \\
      --in-gtf ab_initio.gtf \\
      --out-gtf lgb_scored.gtf \\
      --threshold ${d.lgb_threshold}
  grep -E 'lgb_class "(${keep})"' lgb_scored.gtf > ${prefix}_lgb_filtered.gtf || true
  grep -E 'lgb_class "partial"' lgb_scored.gtf > ${prefix}_lgb_partial.gtf || true
  mv lgb_scored.scores.tsv ${prefix}_lgb_scores.tsv
  """

  stub:
  """
  touch ${prefix}_lgb_filtered.gtf ${prefix}_lgb_partial.gtf ${prefix}_lgb_scores.tsv
  """
}

// Hint rescue, part 1: loci of partial gene finder transcripts without a kept
// transcript on the same strand, each with the hints of its best protein chain.
// Loci where a Drusilla ORF already has all introns of the chain are skipped.
process HINT_RESCUE_LOCI {
  label 'container'
  input:
    path partial_gtf
    path correct_gtf
    path orfs_gtf
    path miniprot_gff
    path hints_gff
    path genome

  output:
    path "rescue/combined_loci.fa", emit: fasta
    path "rescue/combined_hints.gff", emit: hints
    path "rescue/loci_manifest.tsv", emit: manifest

  script:
  """
  samtools faidx ${genome}
  chainedHints.py ${hints_gff} ${miniprot_gff} --output chained_hints.gff
  mkdir -p rescue
  prepare_hint_rescue_loci.py \\
      --partial_gtf ${partial_gtf} \\
      --correct_gtf ${correct_gtf} \\
      --chained_hints chained_hints.gff \\
      --orfs_gtf ${orfs_gtf} \\
      --genome ${genome} \\
      --outdir rescue \\
      --flank ${params.drusilla.rescue_flank}
  touch rescue/combined_loci.fa rescue/combined_hints.gff rescue/loci_manifest.tsv
  """

  stub:
  """
  mkdir -p rescue
  touch rescue/combined_loci.fa rescue/combined_hints.gff rescue/loci_manifest.tsv
  """
}

// Hint rescue, part 2: Tiberius predicts the loci with the chain hints; the
// predictions that agree with the hints are mapped back to the genome.
// Needs a Tiberius with --hints (branch hint_integration). Without it, or
// without loci, the rescue is empty.
process HINT_RESCUE_TIBERIUS {
  label 'gpu', 'container'
  publishDir "${params.outdir}/intermediate", pattern: "hint_rescue.gtf", mode: 'copy'
  input:
    path fasta
    path hints
    path manifest
    val model_cfg

  output:
    path "hint_rescue.gtf", emit: gtf

  script:
  def tiberius = params.drusilla.rescue_tiberius ?: '\$(command -v tiberius.py)'
  // Batch size as in RUN_TIBERIUS: params.tiberius.batch_size, else the int32 cap
  // Loci are short: a small seq_len avoids padding every locus to a genome chunk
  def seqLen = params.drusilla.rescue_seq_len
  def batch = params.tiberius?.batch_size ? "BATCH_ARG='--batch_size ${params.tiberius.batch_size}'" :
      "BATCH_ARG=\$(tiberius_batch_size.py --model_cfg \$(dirname ${tiberius})/model_cfg/${model_cfg}.yaml" +
      (seqLen ? " --seq_len ${seqLen}" : '') + ")"
  """
  ${batch}
  if [ ! -s ${fasta} ]; then
      echo "hint rescue: no loci" >&2
      : > hint_rescue.gtf
  elif ! python3 ${tiberius} --help 2>&1 | grep -q -- '--hints'; then
      echo "WARNING: ${tiberius} has no --hints option; hint rescue skipped" >&2
      : > hint_rescue.gtf
  else
      python3 ${tiberius} \\
          --genome ${fasta} \\
          --model_cfg ${model_cfg} \\
          --hints ${hints} \\
          --hint_weight ${params.drusilla.rescue_hint_weight} \\
          --out rescue_raw.gtf${seqLen ? " --seq_len ${seqLen}" : ''} \${BATCH_ARG:-}
      filter_and_merge_rescue_gtf.py rescue_raw.gtf ${hints} ${manifest} hint_rescue.gtf
  fi
  """

  stub:
  """
  touch hint_rescue.gtf
  """
}
