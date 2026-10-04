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
