// Vipsania (https://github.com/Gaius-Augustus/Vipsania) as gene finder.
//
// These processes carry the label 'vipsania' and never the label 'container':
// they run in the Vipsania image, which conf/base.config assigns to the label.

include { truthy } from '../lib_nf/functions.nf'

// Singularity/Apptainer do not run the ENTRYPOINT of the Vipsania image, which
// points TensorFlow to its bundled CUDA libraries. Do the same here.
def vipsaniaEnv() {
    return '''
    TF_CUDA=$(python -c "import site, glob, os; print(':'.join(sorted(glob.glob(os.path.join(site.getsitepackages()[0], 'nvidia', '*', 'lib')))))" 2>/dev/null || true)
    if [ -n "$TF_CUDA" ]; then export LD_LIBRARY_PATH="$TF_CUDA${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"; fi
    export VIPSANIA_CACHE="$PWD/.vipsania_cache"
    export WANDB_MODE=disabled
    '''.stripIndent()
}

// Download a pretrained model (clade name or model id) and the clade table.
// Needs internet access. Offline: set params.vipsania.model_dir instead.
process DOWNLOAD_VIPSANIA_MODEL {
    label 'vipsania', 'local_only'

    input:
        val model

    output:
        path "models/*"

    script:
    """
    ${vipsaniaEnv()}
    vipsania download "${model}" -d models
    """

    stub:
    """
    mkdir -p models/stub_model
    touch models/versions.json models/stub_model/config.json
    """
}

// Annotate a genome, or a chunk of it. With params.vipsania.finetune the model
// is first finetuned on the given FASTA; the pipeline then passes the whole genome.
process RUN_VIPSANIA {
    label 'gpu', 'vipsania', 'bigmem'
    maxForks params.vipsania?.max_parallel ? (params.vipsania.max_parallel as Integer) : Integer.MAX_VALUE
    publishDir "${params.outdir}/intermediate/vipsania", mode: 'copy', pattern: "{*.log,finetuning_*}"

    input:
        path genome
        path model_files, stageAs: 'vip_models/*'
        val model

    output:
        path "vipsania.${genome.name}.gtf", emit: gtf
        path "vipsania.${genome.name}.log", emit: log, optional: true
        path "finetuning_*",               emit: checkpoint, optional: true

    script:
    def v = params.vipsania ?: [:]
    def extra = ''
    if (v.batch_size) extra += " -B ${v.batch_size}"
    if (v.context)    extra += " -T ${v.context}"
    if (truthy(v.finetune)) {
        extra += " --finetune"
        if (v.finetune_epochs) extra += " --finetune_epochs ${v.finetune_epochs}"
        if (v.finetune_B)      extra += " --finetune_B ${v.finetune_B}"
        if (v.finetune_lr)     extra += " --finetune_lr ${v.finetune_lr}"
    }
    if (v.extra_args) extra += " ${v.extra_args}"
    """
    ${vipsaniaEnv()}
    vipsania annotate "${model}" ${genome} \\
        --model_dir vip_models \\
        -o "vipsania.${genome.name}.gtf"${extra}
    """

    stub:
    """
    touch "vipsania.${genome.name}.gtf" "vipsania.${genome.name}.log"
    """
}
