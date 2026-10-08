// SRA reads: prefetch downloads the run into the task directory and checks it,
// fasterq-dump converts the local copy. fasterq-dump on an accession streams it
// instead and, in a Singularity container, prints hundreds of harmless
// 'storage exhausted ... bad file descriptor' lines although the reads are
// complete (sra-tools 3.3.0). --max-size u lifts the 20 GB limit of prefetch.
// The reads are published to sra_downloads/ only with params.keep_downloads;
// the copy in the work directory is what the pipeline reads, and a second copy
// of dozens of libraries is the largest single cost of a run. The published
// copies are overwritten on -resume, so that a copy cut short by an aborted run
// is replaced.

process DOWNLOAD_SRA_PAIRED {
  tag "${acc}"
  label 'container', 'download'

  publishDir "${params.outdir}/sra_downloads/rnaseq_sra_paired/", mode: 'copy', overwrite: true,
    enabled: params.keep_downloads as boolean

  input:
    val acc

  output:
    tuple val(acc), path("${acc}_1.fastq.gz"), path("${acc}_2.fastq.gz")

  script:
  """
  prefetch --max-size u -O . ${acc}
  fasterq-dump --split-files --threads ${task.cpus} ./${acc}
  rm -r ${acc}
  gzip ${acc}_1.fastq ${acc}_2.fastq
  """

  stub:
  """
  touch ${acc}_1.fastq.gz ${acc}_2.fastq.gz
  """
}

process DOWNLOAD_SRA_SINGLE {
  tag "${acc}"
  label 'container', 'download'

  publishDir "${params.outdir}/sra_downloads/rnaseq_sra_single/", mode: 'copy', overwrite: true,
    enabled: params.keep_downloads as boolean

  input:
    val acc

  output:
    tuple val(acc), path("${acc}.fastq.gz")

  script:
  """
  prefetch --max-size u -O . ${acc}
  fasterq-dump --threads ${task.cpus} ./${acc}
  rm -r ${acc}
  gzip ${acc}.fastq
  """

  stub:
  """
  touch ${acc}.fastq.gz
  """
}

process DOWNLOAD_SRA_ISOSEQ {
  tag "${acc}"
  label 'container', 'download'

  publishDir "${params.outdir}/sra_downloads/isoseq_sra/", mode: 'copy', overwrite: true,
    enabled: params.keep_downloads as boolean

  input:
    val acc

  output:
    tuple val(acc), path("${acc}.fastq.gz")

  script:
  """
  prefetch --max-size u -O . ${acc}
  fasterq-dump --threads ${task.cpus} ./${acc}
  rm -r ${acc}
  gzip ${acc}.fastq
  """

  stub:
  """
  touch ${acc}.fastq.gz
  """
}
