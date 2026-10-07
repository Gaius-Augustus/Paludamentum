nextflow.enable.dsl=2

include { truthy } from '../lib_nf/functions.nf'

// Final protein-coding annotation (docs/postprocessing.md): the sanity filter,
// UTRs from the StringTie assemblies, the published GFF3, GTF and sequences,
// and the longest isoform per gene for the completeness assessment.
// `stem` names the published files: <tool>_evidence or <tool>_ab_initio.

process SANITY_FILTER {
  label 'container'
  publishDir "${params.outdir}/qc", mode:'copy', overwrite: true, pattern: 'sanity_filter.tsv'
  publishDir "${params.outdir}/intermediate", mode:'copy', overwrite: true, pattern: '*_sanity_filtered.gff3'

  input:
    val stem
    path gff
    path genome

  output:
    path "${stem}_sanity_filtered.gff3", emit: gff
    path "sanity_filter.tsv", emit: report

  script:
  // postprocess.sanity_filter = false: report only, nothing is removed
  def keep = truthy(params.postprocess.sanity_filter) ? '' : '--keep'
  """
  sanity_filter_gff3.py \\
      --gff3 ${gff} \\
      --genome ${genome} \\
      --out ${stem}_sanity_filtered.gff3 \\
      --report sanity_filter.tsv ${keep}
  """

  stub:
  """
  touch ${stem}_sanity_filtered.gff3 sanity_filter.tsv
  """
}

// UTRs from the StringTie assemblies for transcripts without UTRs; the CDS
// does not change. Short-read assemblies (and the one assembly of the
// Drusilla flow or of params.stringtie) go to --stringtie, Iso-Seq assemblies
// to --longread, whose matches win.
process ADD_UTRS {
  label 'container'
  publishDir "${params.outdir}/qc", mode:'copy', overwrite: true, pattern: 'utr_report.tsv'

  input:
    val stem
    path gff
    path short_gtfs, stageAs: 'stringtie/?/*'
    path long_gtfs, stageAs: 'longread/?/*'

  output:
    path "${stem}_utr.gff3", emit: gff
    path "utr_report.tsv", emit: report

  script:
  def shortArg = short_gtfs ? "--stringtie ${short_gtfs}" : ''
  def longArg  = long_gtfs ? "--longread ${long_gtfs}" : ''
  """
  add_utrs_from_stringtie.py \\
      --gff3 ${gff} ${shortArg} ${longArg} \\
      --max-utr-extension ${params.postprocess.max_utr_extension} \\
      --out utr.gff3 \\
      --report utr_report.tsv
  normalize_gff3.py --gff3 utr.gff3 --out ${stem}_utr.gff3
  """

  stub:
  """
  touch ${stem}_utr.gff3 utr_report.tsv
  """
}

// The published annotation: <stem>.gff3 (validated against the GFF3
// contract), its GTF, the proteins of all isoforms ('*' for stop codons)
// and the coding sequences.
process FINALIZE_ANNOTATION {
  label 'container'
  publishDir "${params.outdir}/", mode:'copy', overwrite: true, pattern: "${stem}*"

  input:
    val stem
    path gff
    path genome

  output:
    path "${stem}.gff3", emit: gff3
    path "${stem}.gtf", emit: gtf
    path "${stem}_proteins.fa", emit: proteins
    path "${stem}_cds.fa", emit: cds
    path "versions.tsv", emit: versions

  script:
  """
  cp ${gff} ${stem}.gff3
  validate_gff3.sh ${stem}.gff3
  gffread ${stem}.gff3 -T -o ${stem}.gtf
  gffread ${stem}.gff3 -g ${genome} -y ${stem}_proteins.fa -S
  gffread ${stem}.gff3 -g ${genome} -x ${stem}_cds.fa
  printf 'GffRead\\t%s\\t%s\\n' "\$(gffread --version 2>&1)" "${task.container ?: 'none'}" > versions.tsv
  printf 'GenomeTools\\t%s\\t%s\\n' "\$(gt --version 2>&1 | head -n 1 | awk '{print \$NF}')" "${task.container ?: 'none'}" >> versions.tsv
  """

  stub:
  """
  touch ${stem}.gff3 ${stem}.gtf ${stem}_proteins.fa ${stem}_cds.fa versions.tsv
  """
}

// Longest coding isoform per gene and its proteins (BUSCO, compleasm)
process LONGEST_ISOFORM {
  label 'container'
  publishDir "${params.outdir}/qc", mode:'copy', overwrite: true

  input:
    val stem
    path gff
    path genome

  output:
    path "${stem}_longest_isoform.gff3", emit: gff3
    path "${stem}_longest_isoform_proteins.fa", emit: proteins

  script:
  """
  longest_isoform.py --gff3 ${gff} --out ${stem}_longest_isoform.gff3
  validate_gff3.sh ${stem}_longest_isoform.gff3
  gffread ${stem}_longest_isoform.gff3 -g ${genome} -y ${stem}_longest_isoform_proteins.fa -S
  """

  stub:
  """
  touch ${stem}_longest_isoform.gff3 ${stem}_longest_isoform_proteins.fa
  """
}
