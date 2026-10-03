nextflow.enable.dsl=2

process CONCAT_PROTEINS {
  tag "concat_proteins"
  label 'container'

  input:
  path proteins, stageAs:  "?/*"

  output:
  path "proteins_concat.faa"

  script:
  def input_str = proteins instanceof List ? proteins.join(" ") : proteins
  """
    : > proteins_concat.faa
    for f in ${input_str}; do
      if [[ "\$f" == *.gz ]]; then
        gunzip -c "\$f" >> proteins_concat.faa
      else
        cat "\$f" >> proteins_concat.faa
      fi
    done
  """

  stub:
  """
  touch proteins_concat.faa
  """
}

// gunzip a FASTA. The output keeps the input name without .gz, so that a
// decompressed genome and decompressed proteins do not collide downstream.
process DECOMPRESS_FASTA {
  tag { infile.name }
  label 'local_only'

  input:
  path infile

  output:
  path "${outName(infile)}"

  script:
  """
  if [[ "${infile.name}" == *.gz ]]; then
    gunzip -c ${infile} > ${outName(infile)}
  else
    ln -sf \$(readlink -f ${infile}) ${outName(infile)}
  fi
  """

  stub:
  """
  touch ${outName(infile)}
  """
}

def outName(infile) {
  def name = infile.name
  return name.toLowerCase().endsWith('.gz') ? name[0..-4] : "decompressed_${name}"
}

process DOWNLOAD_ODB12_PARTITIONS {
  tag "odb12_partitions"
  label 'download'

  input:
  val partitions

  output:
  path "odb12_partitions.faa"

  script:
  """
  set -euo pipefail
  : > odb12_partitions.faa
  for part in ${partitions}; do
      url="https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/\${part}.fa.gz"
      curl -fsSL "\${url}" | gunzip -c >> odb12_partitions.faa
  done
  """

  stub:
  """
  touch odb12_partitions.faa
  """
}


process CONCAT_HINTS {
  publishDir "${params.outdir}", mode: 'copy'
  label 'container'

  input:
    path(prot)
    path(rnaseq), stageAs: 'rnaseq_hints.gff'
    path(isoseq), stageAs: 'isoseq_hints.gff'

  output:
    path "hintsfile.gff", emit: hints

  script:
  """
  cat ${prot} ${ rnaseq ?: "/dev/null" } ${ isoseq ?: "/dev/null" } > hintsfile.gff
  """

  stub:
  """
  touch hintsfile.gff
  """
}

process EMPTY_FILE {
  label 'local_only'

  output:
    path 'empty.txt'

  script:
  """
    touch empty.txt
  """

  stub:
  """
  touch empty.txt
  """
}

// Percent of mapped reads (samtools flagstat). Runs in the tools image; the
// shell is set to pipefail in base.config, so a missing samtools fails the
// task instead of writing an empty rate.
process CALC_ALIGNMENT_RATE {
    tag "${aln_file.simpleName}"
    label 'container'

    input:
    path aln_file

    output:
    tuple path(aln_file), path("alignment_rate.txt")

    script:
    """
    pct=\$(${params.tools.samtools} flagstat ${aln_file} \\
        | awk '/ mapped \\(/ {
              # take 5th field, e.g. (99.79%
              val=\$5
              # remove everything that is not a digit or a dot
              gsub(/[^0-9.]/, "", val)
              print val
              exit
          }')
    echo "\$pct" > alignment_rate.txt
    """

  stub:
  """
  echo 100 > alignment_rate.txt
  """
}


workflow FILTER_ALIGNMENT {
    take:
    aln_ch

    main:
    CALC_ALIGNMENT_RATE(aln_ch)

    CALC_ALIGNMENT_RATE.out
        .map { file, ratefile ->
            def pct = ratefile.text.trim().toFloat()
            tuple(file, pct)
        }
        .set { file_pct_ch }

    def minRate = (params.min_alignment_rate ?: 80) as float
    file_pct_ch
        .filter { _file, pct -> pct < minRate }
        .view { file, pct ->
            "Removing ${file.simpleName} (alignment ${String.format('%.2f', pct)}% < params.min_alignment_rate ${minRate})"
        }

    emit:
    file_pct_ch
        .filter { _file, pct -> pct >= minRate }
        .map { file, _pct -> file }
}
