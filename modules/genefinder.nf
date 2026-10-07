// Gene finder independent processes (genome splitting, merging, protein
// extraction) and the Tiberius run. Vipsania processes are in vipsania.nf.
// `prefix` is the name of the gene finder; it names the published files.

process RUN_TIBERIUS {
    // Runs in the Tiberius image (label 'tiberius'), not the evidence image
    label 'gpu', 'tiberius', 'bigmem'
    // Cap concurrent GPU tasks; set params.tiberius.max_parallel to e.g. 1
    // on single-GPU workstations to avoid GPU OOM from parallel chunks.
    maxForks params.tiberius?.max_parallel ? (params.tiberius.max_parallel as Integer) : Integer.MAX_VALUE

    input:
        path genome
        path model_cfg
        // use_weights: weights holds the extracted weights of model_cfg (the
        // contents of params.tiberius.model_dir, or of DOWNLOAD_TIBERIUS_WEIGHTS);
        // the weights directory is passed with --model, so Tiberius does not
        // download. false: Tiberius downloads the weights of model_cfg itself.
        val use_weights
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
    // With staged weights, run the same model with --model: with --model_cfg Tiberius
    // ignores the staged weights where its own model_weights directory is
    // writable (Docker) and downloads them. See bin/tiberius_model_args.py.
    // The script prints one argument per line, read into an array.
    def model = use_weights ?
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

// Download and extract the weights archive of a Tiberius model configuration
// (its weights_url) once per run, on the submitting host like the other
// downloads (label 'download'). The GPU tasks get the extracted weights staged,
// so they need no internet and do not each download the archive.
process DOWNLOAD_TIBERIUS_WEIGHTS {
    tag "${url.tokenize('/')[-1]}"
    label 'download'

    input:
        val url

    output:
        path "weights/*"

    script:
    // the directory the archive extracts to, as Tiberius and
    // bin/tiberius_model_args.py name it
    def name = url.tokenize('/')[-1].tokenize('.')[0]
    """
    set -euo pipefail
    mkdir -p weights
    curl -fsSL -o weights.archive "${url}"
    tar -xf weights.archive -C weights
    rm weights.archive
    if [ ! -d "weights/${name}" ]; then
        echo "The weights archive ${url} has no directory ${name}, found: \$(ls weights)" >&2
        exit 1
    fi
    """

    stub:
    def name = url.tokenize('/')[-1].tokenize('.')[0]
    """
    mkdir -p "weights/${name}"
    touch "weights/${name}/weights.h5"
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

// Ab initio predictions of the ab initio mode, before the sanity filter; the
// final <prefix>_ab_initio.gff3 in outdir is written by FINALIZE_ANNOTATION.
process MERGE_GENEFINDER {
    label 'container'
    publishDir "${params.outdir}/intermediate/", mode:'copy', overwrite: true

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
    publishDir "${params.outdir}/intermediate/", mode:'copy', overwrite: true

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

// Ab initio predictions merged with the HC genes, before the sanity filter;
// the final <prefix>_evidence.gff3 is written by FINALIZE_ANNOTATION. The HC
// genes enter the annotation here; nothing is trained on them (the 'TRAIN' of
// the name, and the traingenes input, come from BRAKER, where they train
// AUGUSTUS). Both are kept because the process name and its script body are
// part of the Nextflow task hash: renaming them would invalidate -resume for
// this process and everything downstream of it.
process MERGE_GENEFINDER_TRAIN {
    label 'container'
    publishDir "${params.outdir}/intermediate/", mode:'copy', overwrite: true

    input:
      val prefix
      path ab_initio
      path traingenes

    output:
      path "${prefix}_merged.gff3", emit: merged

    script:
    """
    merge_annotations.py --mode full \\
        ${ab_initio} ${traingenes} > ${prefix}_merged.gff3
    """

    stub:
    """
    touch ${prefix}_merged.gff3
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
    ${params.tools.gffread} ab_initio.capped.gff3 \\
        -g ${genome} \\
        -y ${prefix}_proteins.fa
    """

  stub:
  """
  touch ${prefix}_proteins.fa
  """
}
