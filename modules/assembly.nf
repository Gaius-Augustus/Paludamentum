nextflow.enable.dsl=2

process STRINGTIE_ASSEMBLE_RNA {
  label 'container'
  input:
    path rnabam
  output:
    path "stringtie_${rnabam.baseName}.gtf", emit: gtf
    path "stringtie_${rnabam.baseName}.gff3", emit: gff3
  // The BAMs of RNASEQ_EVIDENCE are coordinate-sorted
  script:
  """
  ${params.tools.stringtie} -p ${task.cpus} -o "stringtie_${rnabam.baseName}.gtf" ${rnabam}
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


// The one StringTie assembly as GTF and as GFF3 (gtf_to_alignment_gff3.pl,
// for TD_REVISED_TO_GENOME). A single input is passed through as it is:
// `stringtie --merge` drops the cov and TPM attributes that FILTER_STRINGTIE
// reads (it stops without them), drops transcripts below TPM 1 (its -T
// default) and renames the transcripts (MSTRG). Several inputs are merged
// (stringtie --merge, merge_annotations.py when that fails).
process STRINGTIE_MERGE {
    label 'container'
    input:
    path gtfs, stageAs: "?/*"

    output:
    path "stringtie.gtf", emit: gtf
    path "stringtie.gff3", emit: gff3

    script:
    """
    set -- ${gtfs}
    if [ \$# -eq 1 ]; then
        cp "\$1" stringtie.gtf
    else
        printf '%s\\n' "\$@" > mergelist.txt

        set +e
        ${params.tools.stringtie} --merge -o stringtie.gtf mergelist.txt
        exitcode=\$?
        set -e

        if [ "\$exitcode" -ne 0 ]; then
            echo "WARNING: stringtie --merge failed with code \$exitcode, using custom merge." >&2

            merge_annotations.py --mode full "\$@" > merged.gff3

            # convert to GTF for TD_ALL
            ${params.tools.gffread} merged.gff3 -T -o stringtie.gtf
        fi
    fi

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
  # the merged BAMs of RNASEQ_EVIDENCE and ISOSEQ_EVIDENCE are coordinate-sorted
  ${params.tools.stringtie} -p ${task.cpus} -o stringtie/stringtie.gtf --mix rna.bam isoseq.bam
  ${params.tools.transdecoder_gtf2gff} stringtie/stringtie.gtf > stringtie/stringtie.gff3
  """

  stub:
  """
  mkdir -p stringtie
  touch stringtie/stringtie.gtf stringtie/stringtie.gff3
  """
}
