nextflow.enable.dsl=2

// One pyVARUS output directory (rnaseq_varus, isoseq_varus, mixed_varus).
// varusInputs in lib_nf/functions.nf has checked the files and the mode; here
// the genome MD5 of the manifest is compared with the genome of the run. The
// MD5 of a gzipped genome is compared both ways (VARUS may have read either).
// Emits the StringTie assembly and the intron hints (none for mixed).
process VARUS_INPUT {
  tag "${dir}"
  label 'container'
  input:
    tuple val(id), val(dir), path(manifest, stageAs: 'varus/manifest.tsv'), path(gtf, stageAs: 'varus/stringtie.gtf'), path(hints, stageAs: 'varus/hints.gff')
    path genome
  output:
    path "${id}.stringtie.gtf", emit: gtf
    path "${id}.hints.gff", emit: hints, optional: true
  script:
  def copyHints = hints ? "cp varus/hints.gff \"${id}.hints.gff\"" : ""
  """
  want=\$(sed -n '/^#genome_md5=/{s///p;q}' varus/manifest.tsv)
  have=\$(md5sum < ${genome} | cut -d' ' -f1)
  if [ "\$want" != "\$have" ] && [[ "${genome.name}" == *.gz ]]; then
      have=\$(gzip -dc ${genome} | md5sum | cut -d' ' -f1)
  fi
  if [ -z "\$want" ]; then
      echo "${dir}: the VARUS manifest has no genome_md5; cannot check the genome." >&2
      exit 1
  fi
  if [ "\$want" != "\$have" ]; then
      echo "${dir}: the VARUS run used another genome (MD5 \$want in its manifest, \$have of ${genome.name}). Pass the genome FASTA that VARUS aligned the reads to." >&2
      exit 1
  fi
  cp varus/stringtie.gtf "${id}.stringtie.gtf"
  ${copyHints}
  """

  stub:
  def copyHints = hints ? "cp varus/hints.gff \"${id}.hints.gff\"" : ""
  """
  cp varus/stringtie.gtf "${id}.stringtie.gtf"
  ${copyHints}
  """
}

// Intron hints of several sources (pyVARUS directories, BAM2HINTS of the
// merged BAM) as one file, equal to bam2hints on all reads.
process MERGE_INTRON_HINTS {
  label 'container'
  input:
    path hints, stageAs: "?/*"
  output:
    path "intron_hints.gff", emit: hints
  script:
  """
  merge_intron_hints.py ${hints} > intron_hints.gff
  """

  stub:
  """
  touch intron_hints.gff
  """
}

// pyVARUS inside the pipeline (params.varus.run with params.varus.species; off
// by default): `varus runlist` asks NCBI SRA for the short-read RNA-Seq runs of
// the species, `varus index` builds a HISAT2 index, and `varus run` samples
// reads until the genes are covered (downloads batches of the runs and aligns
// them). --keep-bam keeps VARUS.bam, which the pipeline reads like a BAM of
// rnaseq_bam (strand fix, merge, StringTie, intron hints, coverage track).
// The task needs internet access (NCBI Entrez, SRA; Logan's S3 bucket with
// varus.logan) and the CPUs and memory of an alignment task: it runs on the
// default executor, not on the submitting host like the downloads. A site
// config moves it with withName: 'RUN_PYVARUS' (docs/hpc.md).
// Published to <outdir>/varus/: the manifest and the splice-site DB log, from
// which `varus replay` rebuilds the BAM, the run list, the statistics and
// pyVARUS' own CITATIONS.md; the BAM itself only with params.keep_downloads.
// `varus run` exits with 3 when no batch passed its quality gate (no usable
// RNA-Seq of the species in SRA).
process RUN_PYVARUS {
  tag "${species}"
  label 'pyvarus'
  // hisat2-build needs about 8 x the genome FASTA in memory (pyVARUS, measured on wheat)
  memory { "${Math.max(32, Math.ceil(9 * genome.size() / 1e9 + 8) as int) * task.attempt} GB" }
  time '72h'

  publishDir "${params.outdir}/varus", mode: 'copy', overwrite: true,
    pattern: "{VARUS.manifest.tsv,VARUS.splicedb.log.gz,CITATIONS.md,Runlist.tsv,RunStatistics.csv,Coverage.csv}"
  publishDir "${params.outdir}/varus", mode: 'copy', overwrite: true, pattern: "pyvarus_run.bam",
    enabled: params.keep_downloads as boolean

  input:
    path genome
    val species

  output:
    path "pyvarus_run.bam", emit: bam
    path "VARUS.manifest.tsv", emit: manifest
    path "VARUS.splicedb.log.gz", emit: splicedb
    path "Runlist.tsv", emit: runlist
    path "CITATIONS.md", optional: true
    path "RunStatistics.csv", optional: true
    path "Coverage.csv", optional: true

  script:
  def v = params.varus ?: [:]
  def email = v.email ? "--email '${v.email}'" : ''
  def maxRuns = v.max_runs ? "--max-runs ${v.max_runs as int}" : ''
  def maxBatches = (v.max_batches ?: 1000) as int
  def logan = (v.logan == true || v.logan?.toString()?.toLowerCase() == 'true') ? '--logan' : ''
  def extra = v.extra_args ?: ''
  """
  set -euo pipefail
  varus runlist '${species}' --outdir . ${email} ${maxRuns}
  if ! grep -q '^[^@#]' Runlist.tsv 2>/dev/null; then
      echo "pyVARUS: NCBI SRA lists no short-read RNA-Seq run of '${species}' (varus.species). Check the name, or give reads, a BAM or SRA accessions instead." >&2
      exit 1
  fi
  varus index ${genome} --outdir genome_index --threads ${task.cpus}
  set +e
  varus run '${species}' ${genome} \\
      --runlist Runlist.tsv --index genome_index/hisatidx --outdir . \\
      --threads ${task.cpus} --max-batches ${maxBatches} --keep-bam ${logan} ${extra}
  rc=\$?
  set -e
  if [ "\$rc" -eq 3 ]; then
      echo "pyVARUS: no batch of RNA-Seq reads of '${species}' passed the quality gate (varus run exit status 3): SRA holds no usable RNA-Seq for this genome. Give reads, a BAM or SRA accessions instead, or run without varus.run." >&2
      exit 1
  fi
  [ "\$rc" -eq 0 ] || exit "\$rc"
  mv VARUS.bam pyvarus_run.bam
  # the HISAT2 index and the per-batch directories are not needed any more
  rm -r genome_index
  if [ -d batches ]; then rm -r batches; fi
  """

  stub:
  """
  touch pyvarus_run.bam VARUS.manifest.tsv VARUS.splicedb.log.gz Runlist.tsv CITATIONS.md
  """
}
