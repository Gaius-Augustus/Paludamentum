process DOWNLOAD_SRA_PAIRED {
  label 'container', 'download'

  publishDir "${params.outdir}/sra_downloads/rnaseq_sra_paired/", mode: 'copy'

  input:
    val acc

  output:
    tuple val(acc), path("${acc}_1.fastq.gz"), path("${acc}_2.fastq.gz")

  script:
  """
  fasterq-dump --split-files --threads ${task.cpus} ${acc}
  gzip ${acc}_1.fastq ${acc}_2.fastq
  """

  stub:
  """
  touch ${acc}_1.fastq.gz ${acc}_2.fastq.gz
  """
}

process DOWNLOAD_SRA_SINGLE {
  label 'container', 'download'

  publishDir "${params.outdir}/sra_downloads/rnaseq_sra_single/", mode: 'copy'

  input:
    val acc

  output:
    tuple val(acc), path("${acc}.fastq.gz")

  script:
  """
  fasterq-dump --threads ${task.cpus} ${acc}
  gzip ${acc}.fastq
  """

  stub:
  """
  touch ${acc}.fastq.gz
  """
}

process DOWNLOAD_SRA_ISOSEQ {
  label 'container', 'download'

  publishDir "${params.outdir}/sra_downloads/isoseq_sra/", mode: 'copy'

  input:
    val acc

  output:
    tuple val(acc), path("${acc}.fastq.gz")

  script:
  """
  fasterq-dump --threads ${task.cpus} ${acc}
  gzip ${acc}.fastq
  """

  stub:
  """
  touch ${acc}.fastq.gz
  """
}
