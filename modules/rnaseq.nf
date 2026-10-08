nextflow.enable.dsl=2
process HISAT2_BUILD {
  label 'container', 'bigmem'
  // hisat2-build segfaults with very high thread counts on large genomes;
  // cap at 32, but honour a lower params.threads if the user set one.
  cpus { Math.min((params.threads ?: 32) as Integer, 32) }

  input:
    path genome

  output:
    path "hisat2_idx", emit: idxdir

  script:
  """
  mkdir -p hisat2_idx
  ${params.tools.hisat2_build} -p ${task.cpus} ${genome} hisat2_idx/genome
  """

  stub:
  """
  mkdir -p hisat2_idx
  """
}

process HISAT2_MAP_SINGLE {
  label 'container', 'bigmem'
  input:
    path idxdir
    path reads
  output:
    path "${reads.simpleName}.bam", emit: bam

  // No --dta (here and in HISAT2_MAP_PAIRED): its longer junction anchors cost
  // 12-21 % of the intron hints and did not improve gene F1 (pyVARUS
  // benchmark, T. rubripes and B. taurus, 2026-10-03)
  script:
  """
  ${params.tools.hisat2} -x ${idxdir}/genome -U ${reads} -p ${task.cpus} \
    | ${params.tools.samtools} sort -@ ${task.cpus} -o "${reads.simpleName}.bam"
  """

  stub:
  """
  touch "${reads.simpleName}.bam"
  """
}

process HISAT2_MAP_PAIRED {
  label 'container', 'bigmem'
  input:
    path idxdir
    tuple val(sample), path(reads)
  output:
    path "${sample}.bam", emit: bam

  script:
  """
  ${params.tools.hisat2} -x ${idxdir}/genome -1 ${reads[0]} -2 ${reads[1]} -p ${task.cpus} \
    | ${params.tools.samtools} sort -@ ${task.cpus} -o "${sample}.bam"
  """

  stub:
  """
  touch "${sample}.bam"
  """
}

// User BAMs (rnaseq_bam) whose header does not say SO:coordinate. The other
// BAMs of the pipeline are coordinate-sorted when they are made (HISAT2,
// minimap2), so merging, bam2hints and StringTie take them as they are.
process SAMTOOLS_SORT {
  tag "${bam.name}"
  label 'container', 'bigmem'
  input:
    path bam, stageAs: "unsorted/*"

  output:
    path "${bam.name}", emit: bam

  script:
  """
  ${params.tools.samtools} sort -@ ${task.cpus} -o "${bam.name}" "${bam}"
  """

  stub:
  """
  touch "${bam.name}"
  """
}

process SAMTOOLS_MERGE {
  label 'container', 'bigmem'
  input:
    // BAMs of different libraries can carry the same file name (Aligned.out.bam)
    path bams, stageAs: "?/*"

  output:
    path "merged.bam", emit: bam

  script:
  """
  ${params.tools.samtools} merge -f -@ ${task.cpus} merged.bam ${bams}
  """

  stub:
  """
  touch merged.bam
  """
}

// One BAM of all coordinate-sorted BAMs, as a value channel. A single BAM is
// passed on as it is; samtools merge would only copy it.
workflow MERGE_BAMS {
  take:
    bams

  main:
    def all = bams.collect().branch { b ->
        one:  b.size() == 1
        many: true
    }
    def merged = SAMTOOLS_MERGE(all.many).bam.mix(all.one.map { b -> b[0] }).first()

  emit:
    bam = merged
}


// The BAM is coordinate-sorted (MERGE_BAMS).
process BAM2HINTS {
  label 'container', 'bigmem'
  input: path bam; path genome
  output: path "${bam.simpleName}.hints.gff", emit: hints
  script: """
  ${params.tools.bam2hints} --intronsonly --in=${bam} --out=introns.temp
  filterIntronsFindStrand.pl ${genome} introns.temp --score > "${bam.simpleName}.hints.gff"
  """

  stub:
  """
  touch "${bam.simpleName}.hints.gff"
  """
}
