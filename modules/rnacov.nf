nextflow.enable.dsl=2

// RNA-Seq coverage track for combine_genesets.py: bedGraph with one line per
// region of constant depth (bedtools genomecov -bg -split), gzipped. -split
// counts N-cigar segments as introns, which matches bam2hints introns and the
// expectation of combine_genesets.py --coverage.
process BAM2COV {
  label 'container'
  publishDir "${params.outdir}/intermediate", pattern: "rna_cov.bg.gz", mode: 'copy', overwrite: true

  input:
    path bam

  output:
    path "rna_cov.bg.gz", emit: bg

  script:
  """
  ${params.tools.bedtools} genomecov -ibam ${bam} -bg -split | gzip > rna_cov.bg.gz
  """

  stub:
  """
  : | gzip > rna_cov.bg.gz
  """
}
