nextflow.enable.dsl=2

process MINIMAP2_MAP {
  label 'container', 'bigmem'
  input: path genome; path reads
  output: path "isoseq/${reads.baseName}.bam", emit: bam
  script: """
  mkdir -p isoseq
  ${params.tools.minimap2} -ax splice:hq -uf ${genome} ${reads} -t ${task.cpus} \
    | ${params.tools.samtools} sort -@ ${task.cpus} -o isoseq/${reads.baseName}.bam
  """

  stub:
  """
  mkdir -p isoseq
  touch isoseq/${reads.baseName}.bam
  """
}
