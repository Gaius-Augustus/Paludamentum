nextflow.enable.dsl=2

// GO terms with FANTASIA-Lite (params.fantasia.run, GPU; docs/postprocessing.md).
// The commands follow BRAKER4 rules/postprocessing/run_fantasia.smk (commit
// 3535ed3). Nextflow runs the FANTASIA-Lite image itself (label fantasia);
// conf/base.config binds bin/fantasia_generate_embeddings.py over the image's
// generate_embeddings.py, which loads ProtT5 from PALUDAMENTUM_HF_MODEL_PATH.

process FANTASIA_ANNOTATE {
  label 'fantasia', 'gpu'
  publishDir "${params.outdir}/qc", mode:'copy', overwrite: true, pattern: 'fantasia/{results.csv,failed_sequences.csv}'

  input:
    path proteins
    path hf_cache
    path lookup_dir

  output:
    path "fantasia/results.csv", emit: results
    path "fantasia/failed_sequences.csv", optional: true, emit: failed

  script:
  def extra = params.fantasia.additional_params ?: ''
  """
  if ! command -v nvidia-smi > /dev/null 2>&1 || ! nvidia-smi -L > /dev/null 2>&1; then
      echo "FANTASIA-Lite needs a CUDA GPU, but no GPU is visible on \$(hostname). Give the label gpu a GPU queue in your config, or set fantasia.run = false." >&2
      exit 1
  fi
  # ProtT5-XL needs about 14 GB of GPU memory
  gpu=\${CUDA_VISIBLE_DEVICES:-0}
  gpu=\${gpu%%,*}
  free=\$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits -i "\$gpu" 2>/dev/null | tr -d ' ' || echo 0)
  if [ "\${free:-0}" -lt 15000 ] 2>/dev/null; then
      echo "GPU \$gpu has \$free MiB free; ProtT5-XL needs at least 15000 MiB." >&2
      exit 1
  fi
  # The snapshot with pytorch_model.bin: torch.load reads it into memory,
  # the safetensors mmap of the other snapshot fails with SIGBUS in
  # Singularity on some clusters
  hf=\$(readlink -f ${hf_cache})/hub/models--Rostlab--prot_t5_xl_uniref50
  snapshot=\$(find "\$hf/snapshots" -maxdepth 2 -name pytorch_model.bin 2>/dev/null | head -n 1 | xargs -r dirname)
  if [ -z "\$snapshot" ]; then
      rev=\$(tr -d '[:space:]' < "\$hf/refs/main" 2>/dev/null || true)
      [ -n "\$rev" ] || { echo "No ProtT5 model (pytorch_model.bin or refs/main) in \$hf" >&2; exit 1; }
      snapshot=\$hf/snapshots/\$rev
  fi
  lookup=\$(readlink -f ${lookup_dir})
  filter_proteins.py --in ${proteins} --out proteins.fa --strip-stop
  mkdir -p fantasia
  export PALUDAMENTUM_HF_MODEL_PATH=\$snapshot TRANSFORMERS_OFFLINE=1 HF_HUB_OFFLINE=1
  python3 /opt/fantasia-lite/src/fantasia_pipeline.py \\
      --serial-models \\
      --embed-models prot_t5 \\
      --device cuda \\
      --venv-dir /opt/venv \\
      --lookup-npz "\$lookup/lookup_table.npz" \\
      --annotations-json "\$lookup/annotations.json" \\
      --accessions-json "\$lookup/accessions.json" \\
      --embeddings-npz fantasia/query_embeddings.npz \\
      --config-yaml fantasia/fantasia_config.yaml \\
      --results-csv fantasia/results.csv \\
      --topgo \\
      --topgo-dir fantasia/topgo \\
      --chunk-dir fantasia/tmp/fasta_chunks \\
      --chunk-embed-dir fantasia/tmp/chunk_embeddings \\
      --chunk-results-dir fantasia/tmp/chunk_results \\
      --chunk-config-dir fantasia/tmp/chunk_configs \\
      --chunk-failure-dir fantasia/tmp/failures \\
      --failure-report fantasia/failed_sequences.csv \\
      ${extra} \\
      proteins.fa
  rm -rf fantasia/tmp fantasia/query_embeddings.npz
  """

  stub:
  """
  mkdir -p fantasia
  touch fantasia/results.csv
  """
}

process FANTASIA_SUMMARY {
  label 'postprocess'
  publishDir "${params.outdir}/qc/fantasia", mode:'copy', overwrite: true

  input:
    path results

  output:
    path "fantasia_summary.txt", emit: summary
    path "fantasia_go_terms.tsv", emit: go_terms
    path "fantasia_go_categories.png", emit: plot

  script:
  """
  fantasia_summary.py --results ${results} --out-dir . --min-score ${params.fantasia.min_score}
  """

  stub:
  """
  touch fantasia_summary.txt fantasia_go_terms.tsv fantasia_go_categories.png
  """
}

// Ontology_term=GO:... on the mRNAs and genes of a GFF3: <name>_go.gff3
process FANTASIA_DECORATE {
  tag "${name}"
  label 'container'
  publishDir "${params.outdir}/", mode:'copy', overwrite: true

  input:
    val name
    path gff3
    path results

  output:
    path "${name}_go.gff3", emit: gff3

  script:
  """
  fantasia_decorate_gff3.py --gff3-in ${gff3} --gff3-out ${name}_go.gff3 --results ${results} \\
      --min-score ${params.fantasia.min_score}
  validate_gff3.sh ${name}_go.gff3
  """

  stub:
  """
  touch ${name}_go.gff3
  """
}
