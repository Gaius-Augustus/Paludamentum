// Gene finder independent processes (genome splitting, merging, protein
// extraction) and the Tiberius run. Vipsania processes are in vipsania.nf.
// `prefix` is the name of the gene finder; it names the published files.

process RUN_TIBERIUS {
    label 'gpu', 'container', 'bigmem'
    // Cap concurrent GPU tasks; set params.tiberius.max_parallel to e.g. 1
    // on single-GPU workstations to avoid GPU OOM from parallel chunks.
    maxForks params.tiberius?.max_parallel ? (params.tiberius.max_parallel as Integer) : Integer.MAX_VALUE

    input:
        path genome
        path model_cfg
        // contents of params.tiberius.model_dir, staged under their own names;
        // the weights directory is passed with --model, so Tiberius does not
        // download. Empty list = Tiberius downloads the weights of model_cfg.
        path weights

    output:
        path "tiberius.${genome.name}.gtf"

    script:
    def extra = ''
    if (params.tiberius?.batch_size) extra += " --batch_size ${params.tiberius.batch_size}"
    if (params.tiberius?.seq_len)    extra += " --seq_len ${params.tiberius.seq_len}"
    // Without a configured batch size, cap Tiberius' automatic choice where it
    // would overflow int32 on large GPUs (e.g. 96 GB); see bin/tiberius_batch_size.py.
    def cap = params.tiberius?.batch_size ? '' :
        "BATCH_ARG=\$(tiberius_batch_size.py --model_cfg ${model_cfg}" +
        (params.tiberius?.seq_len ? " --seq_len ${params.tiberius.seq_len}" : '') + ")"
    // With model_dir, run the same model with --model: with --model_cfg Tiberius
    // ignores the staged weights where its own model_weights directory is
    // writable (Docker) and downloads them. See bin/tiberius_model_args.py.
    // The script prints one argument per line, read into an array.
    def model = params.tiberius?.model_dir ?
        "MODEL_OUT=\$(tiberius_model_args.py --model_cfg ${model_cfg}" +
        (params.tiberius?.seq_len ? " --seq_len ${params.tiberius.seq_len}" : '') + "); " +
        "mapfile -t MODEL_ARGS <<< \"\$MODEL_OUT\"" :
        "MODEL_ARGS=(--model_cfg ${model_cfg})"
    """
    ${cap}
    ${model}
    tiberius.py \\
        --genome ${genome} \\
        "\${MODEL_ARGS[@]}" \\
        --out "tiberius.${genome.name}.gtf"${extra} \${BATCH_ARG:-}
    """

    stub:
    """
    touch "tiberius.${genome.name}.gtf"
    """
}

process SPLIT_GENOME {
    label 'container'
    input:
      path genome
      val min_size
      val max_files

    output:
      path "chunks/*.fa", emit: chunks

    script:
    """
    mkdir -p chunks
    split_genome_fasta.py \\
        --genome ${genome} \\
        --outdir chunks \\
        --prefix genome_chunk \\
        --min-size ${min_size} \\
        --max-files ${max_files}
    """

    stub:
    """
    mkdir -p chunks
    cp ${genome} chunks/genome_chunk_1.fa
    cp ${genome} chunks/genome_chunk_2.fa
    """
}

// Ab initio predictions as final result (ab initio mode).
process MERGE_GENEFINDER {
    label 'container'
    publishDir "${params.outdir}/", mode:'copy'

    input:
      val prefix
      path gff_files, stageAs: "?/*"

    output:
      path "${prefix}_ab_initio.gff3"

    script:
    """
    merge_annotations.py --mode full \\
        ${gff_files} > ${prefix}_ab_initio.gff3
    """

    stub:
    """
    touch ${prefix}_ab_initio.gff3
    """
}

// Ab initio predictions as intermediate result (evidence modes).
process MERGE_GENEFINDER_EVI {
    label 'container'
    publishDir "${params.outdir}/intermediate/", mode:'copy'

    input:
      val prefix
      path gff_files, stageAs: "?/*"

    output:
      path "${prefix}_ab_initio.gff3"

    script:
    """
    merge_annotations.py --mode full \\
        ${gff_files} > ${prefix}_ab_initio.gff3
    """

    stub:
    """
    touch ${prefix}_ab_initio.gff3
    """
}

process MERGE_GENEFINDER_TRAIN {
    label 'container'
    publishDir "${params.outdir}/", mode:'copy'

    input:
      val prefix
      path ab_initio
      path traingenes

    output:
      path "${prefix}_evidence.gff3", emit: merged

    script:
    """
    merge_annotations.py --mode full \\
        ${ab_initio} ${traingenes} > ${prefix}_evidence.gff3
    """

    stub:
    """
    touch ${prefix}_evidence.gff3
    """
}

process PROTEIN_FROM_GFF {
  label 'container'

  input:
      val prefix
      path ab_initio
      path genome

  output:
      path "${prefix}_proteins.fa"

  script:
    """
    # gffread malloc()s on genes with many isoforms; cap to a safe number
    # because the resulting protein FASTA is only used by DIAMOND species
    # ranking, so one representative per gene is enough.
    cap_isoforms_per_gene.py --max 10 ${ab_initio} > ab_initio.capped.gff3
    gffread ab_initio.capped.gff3 \\
        -g ${genome} \\
        -y ${prefix}_proteins.fa
    """

  stub:
  """
  touch ${prefix}_proteins.fa
  """
}

process PROTEIN_FROM_GFF_FINAL {
  publishDir "${params.outdir}/", mode:'copy'

  label 'container'

  input:
      val prefix
      path annotation
      path genome

  output:
      path "${prefix}_evidence_proteins.fa"

  script:
    """
    gffread ${annotation} \\
        -g ${genome} \\
        -y ${prefix}_evidence_proteins.fa
    """

  stub:
  """
  touch ${prefix}_evidence_proteins.fa
  """
}
