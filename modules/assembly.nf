nextflow.enable.dsl=2

process STRINGTIE_ASSEMBLE_RNA {
  label 'container'
  input:
    path rnabam
  output:
    path "stringtie_${rnabam.baseName}.gtf", emit: gtf
    path "stringtie_${rnabam.baseName}.gff3", emit: gff3
  script:
  """
  ${params.tools.samtools} sort -@ ${task.cpus} -o sorted.bam ${rnabam}
  ${params.tools.stringtie} -p ${task.cpus} -o "stringtie_${rnabam.baseName}.gtf" sorted.bam
  ${params.tools.transdecoder_gtf2gff} "stringtie_${rnabam.baseName}.gtf" > "stringtie_${rnabam.baseName}.gff3"
  """

  stub:
  """
  touch "stringtie_${rnabam.baseName}.gtf" "stringtie_${rnabam.baseName}.gff3"
  """
}

process STRINGTIE_ASSEMBLE_ISO {
  label 'container'
  input:
    path isobam
  output:
    path "stringtie_${isobam.baseName}.gtf", emit: gtf
    path "stringtie_${isobam.baseName}.gff3", emit: gff3
  script:
  """
  ${params.tools.stringtie} -p ${task.cpus} -o "stringtie_${isobam.baseName}.gtf" -L ${isobam}
  ${params.tools.transdecoder_gtf2gff} "stringtie_${isobam.baseName}.gtf" > "stringtie_${isobam.baseName}.gff3"
  """

  stub:
  """
  touch "stringtie_${isobam.baseName}.gtf" "stringtie_${isobam.baseName}.gff3"
  """
}


process STRINGTIE_MERGE {
    label 'container'
    input:
    path gtfs, stageAs: "?/*"

    output:
    path "stringtie.gtf", emit: gtf
    path "stringtie.gff3", emit: gff3

    script:
    """
    ls ${gtfs} > mergelist.txt

    set +e
    ${params.tools.stringtie} --merge -o stringtie.gtf mergelist.txt
    exitcode=\$?
    set -e

    if [ "\$exitcode" -ne 0 ]; then
        echo "WARNING: stringtie --merge failed with code \$exitcode, using custom merge." >&2

        merge_annotations.py --mode full \\
            ${gtfs} > merged.gff3

        # convert to GTF for TD_ALL
        gffread merged.gff3 -T -o stringtie.gtf
    fi

    # in both cases we now have stringtie.gtf
    # (and we can make a GFF3 if needed)
    ${params.tools.transdecoder_gtf2gff} stringtie.gtf > stringtie.gff3
    """

  stub:
  """
  touch stringtie.gtf stringtie.gff3
  """
}


process STRINGTIE_ASSEMBLE_MIX {
  label 'container'
  input:
    path rnabam,  stageAs: 'rna.bam'
    path isobam,  stageAs: 'isoseq.bam'
  output:
    path "stringtie/stringtie.gtf",  emit: gtf
    path "stringtie/stringtie.gff3", emit: gff3
  script:
  """
  mkdir -p stringtie
  # StringTie needs coordinate-sorted BAMs; user BAMs may be unsorted
  ${params.tools.samtools} sort -@ ${task.cpus} -o rna.sorted.bam rna.bam
  ${params.tools.samtools} sort -@ ${task.cpus} -o isoseq.sorted.bam isoseq.bam
  ${params.tools.stringtie} -p ${task.cpus} -o stringtie/stringtie.gtf --mix rna.sorted.bam isoseq.sorted.bam
  ${params.tools.transdecoder_gtf2gff} stringtie/stringtie.gtf > stringtie/stringtie.gff3
  """

  stub:
  """
  mkdir -p stringtie
  touch stringtie/stringtie.gtf stringtie/stringtie.gff3
  """
}
