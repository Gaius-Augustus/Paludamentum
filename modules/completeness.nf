nextflow.enable.dsl=2

// Completeness of the genome and of the final proteome (docs/postprocessing.md):
// compleasm and BUSCO with one lineage (params.qc.busco_lineage). Both read
// the lineage from one cache directory (params.qc.busco_download_path,
// default ~/.cache/paludamentum/busco) with the layout
// <dir>/lineages/<lineage>/: compleasm with -L <dir>/lineages, BUSCO with
// --download_path <dir> --offline.

// The lineage, once, on the submitting host (label 'download'), written
// straight into the cache directory (staged as `cache`): compleasm downloads
// the BUSCO lineage tarball, the file_versions.tsv, the eukaryota lineage
// and the placement files it checks for (with their .done markers), and the
// tarball is then extracted completely, so that BUSCO finds dataset.cfg.
process DOWNLOAD_BUSCO_LINEAGE {
  tag "${lineage}"
  label 'download', 'postprocess'

  input:
    val lineage
    path cache

  output:
    path "downloaded.txt", emit: done

  script:
  def name = lineage.tokenize('_')[0]
  def odb  = lineage.tokenize('_')[1]
  """
  set -euo pipefail
  mkdir -p ${cache}/lineages
  compleasm_wrapper.py "\$(command -v ${params.tools.compleasm})" download ${name} --odb ${odb} -L ${cache}/lineages
  for tarball in ${cache}/lineages/${lineage}.*.tar.gz; do
      [ -e "\$tarball" ] || continue
      tar -xzf "\$tarball" -C ${cache}/lineages
      rm -f "\$tarball"
  done
  rm -f ${cache}/lineages/eukaryota_${odb}.*.tar.gz
  if [ ! -f ${cache}/lineages/${lineage}/dataset.cfg ]; then
      echo "No ${cache}/lineages/${lineage}/dataset.cfg after the download of ${lineage}." >&2
      exit 1
  fi
  echo "${lineage}" > downloaded.txt
  """

  stub:
  """
  echo "${lineage}" > downloaded.txt
  """
}

// Proteins for hmmsearch: no record of 100,000 aa or more (hmmsearch refuses
// them), no stop codon symbols
process FILTER_BUSCO_PROTEINS {
  label 'postprocess'

  input:
    path proteins

  output:
    path "busco_proteins.fa", emit: proteins

  script:
  """
  filter_proteins.py --in ${proteins} --out busco_proteins.fa --max-length 100000 --strip-stop
  """

  stub:
  """
  touch busco_proteins.fa
  """
}

process COMPLEASM_GENOME {
  label 'postprocess', 'bigmem'
  publishDir "${params.outdir}/qc", mode:'copy', overwrite: true, pattern: 'compleasm_genome/summary.txt'

  input:
    path genome
    val lineage
    path cache

  output:
    path "compleasm_genome/summary.txt", emit: summary
    path "versions.tsv", emit: versions

  script:
  def name = lineage.tokenize('_')[0]
  def odb  = lineage.tokenize('_')[1]
  """
  compleasm_wrapper.py "\$(command -v ${params.tools.compleasm})" run -a ${genome} -o compleasm_genome -t ${task.cpus} \\
      -l ${name} --odb ${odb} -L ${cache}/lineages
  # miniprot alignments and hmmsearch output are large; only the summary is kept
  rm -rf compleasm_genome/*/hmmer_output compleasm_genome/*/miniprot_output.gff
  printf 'compleasm\\t%s\\t%s\\n' "\$(${params.tools.compleasm} --version 2>&1 | awk '{print \$NF}')" "${task.container ?: 'none'}" > versions.tsv
  """

  stub:
  """
  mkdir -p compleasm_genome
  touch compleasm_genome/summary.txt versions.tsv
  """
}

process COMPLEASM_PROTEINS {
  label 'postprocess'
  publishDir "${params.outdir}/qc", mode:'copy', overwrite: true, pattern: 'compleasm_proteins/summary.txt'

  input:
    path proteins
    val lineage
    path cache

  output:
    path "compleasm_proteins/summary.txt", emit: summary

  script:
  def name = lineage.tokenize('_')[0]
  def odb  = lineage.tokenize('_')[1]
  """
  compleasm_wrapper.py "\$(command -v ${params.tools.compleasm})" protein -p ${proteins} -o compleasm_proteins -t ${task.cpus} \\
      -l ${name} --odb ${odb} -L ${cache}/lineages
  rm -rf compleasm_proteins/*/hmmer_output
  """

  stub:
  """
  mkdir -p compleasm_proteins
  touch compleasm_proteins/summary.txt
  """
}

process BUSCO_GENOME {
  label 'busco', 'bigmem'
  publishDir "${params.outdir}/qc", mode:'copy', overwrite: true, pattern: 'busco_genome_short_summary.txt'

  input:
    path genome
    val lineage
    path cache

  output:
    path "busco_genome_short_summary.txt", emit: summary
    path "versions.tsv", emit: versions

  script:
  """
  busco -i ${genome} -o genome --out_path . -l ${lineage} -m genome -c ${task.cpus} \\
      --download_path ${cache} --offline
  cp genome/short_summary*.txt busco_genome_short_summary.txt
  rm -rf genome
  printf 'BUSCO\\t%s\\t%s\\n' "\$(busco --version 2>&1 | head -n 1 | awk '{print \$NF}')" "${task.container ?: 'none'}" > versions.tsv
  """

  stub:
  """
  touch busco_genome_short_summary.txt versions.tsv
  """
}

process BUSCO_PROTEINS {
  label 'busco'
  publishDir "${params.outdir}/qc", mode:'copy', overwrite: true, pattern: 'busco_proteins_short_summary.txt'

  input:
    path proteins
    val lineage
    path cache

  output:
    path "busco_proteins_short_summary.txt", emit: summary

  script:
  """
  busco -i ${proteins} -o proteins --out_path . -l ${lineage} -m proteins -c ${task.cpus} \\
      --download_path ${cache} --offline
  cp proteins/short_summary*.txt busco_proteins_short_summary.txt
  rm -rf proteins
  """

  stub:
  """
  touch busco_proteins_short_summary.txt
  """
}

// qc/completeness.tsv from the summaries of the assessments that ran; the
// others are empty placeholder files
process COMPLETENESS_SUMMARY {
  label 'postprocess'
  publishDir "${params.outdir}/qc", mode:'copy', overwrite: true

  input:
    path busco_genome, stageAs: 'in/busco_genome.txt'
    path busco_proteins, stageAs: 'in/busco_proteins.txt'
    path compleasm_genome, stageAs: 'in/compleasm_genome.txt'
    path compleasm_proteins, stageAs: 'in/compleasm_proteins.txt'

  output:
    path "completeness.tsv", emit: tsv

  script:
  """
  completeness_summary.py --out completeness.tsv \\
      --busco-genome ${busco_genome} --busco-proteins ${busco_proteins} \\
      --compleasm-genome ${compleasm_genome} --compleasm-proteins ${compleasm_proteins}
  """

  stub:
  """
  touch completeness.tsv
  """
}
