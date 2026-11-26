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
    publishDir "${params.outdir}/intermediate/vipsania", mode: 'copy', overwrite: true, pattern: "{*.log,finetuning_*}"

    input:
        path genome
        path model_files, stageAs: 'vip_models/*'
        val model

    output:
        path "vipsania.${genome.name}.gtf", emit: gtf
        path "vipsania.${genome.name}.log", emit: log, optional: true
        path "finetuning_*",               emit: checkpoint, optional: true
        path "versions.tsv",               emit: versions

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
    printf 'Vipsania\\t%s\\t%s\\n' "\$(python3 -c 'import importlib.metadata as m; print(m.version("vipsania"))' 2>/dev/null || echo unknown)" "${task.container ?: 'none'}" > versions.tsv
    """

    stub:
    """
    touch "vipsania.${genome.name}.gtf" "vipsania.${genome.name}.log" versions.tsv
    """
}

// A genome whose genetic code is not the standard one (params.translation_table):
// the pretrained Vipsania models are table-1 models and the code is part of
// the weights, so a model is trained from scratch on the target genome
// (unsupervised: it reads the sequence alone). The training set are the
// sequences of at least vipsania.train_min_contig bp (Vipsania skips
// sequences below its training window of 20 kb). ckpt/ holds the run
// directory with the weights, MODEL_ID its name.
process TRAIN_VIPSANIA {
    label 'gpu', 'vipsania', 'bigmem'
    publishDir "${params.outdir}/intermediate/vipsania", mode: 'copy', overwrite: true, pattern: "{vipsania_train.log,train_set.tsv}"

    input:
        path genome
        val table

    output:
        path "ckpt",               emit: ckpt
        path "MODEL_ID",           emit: model_id
        path "vipsania_train.log", emit: log
        path "train_set.tsv",      emit: train_set

    script:
    def v = params.vipsania ?: [:]
    def config = v.train_config ?: 'configs/base_10M.json'
    def epochs = (v.train_epochs ?: 300) as Integer
    def steps  = (v.train_steps ?: 800) as Integer
    def minLen = (v.train_min_contig ?: 20000) as Integer
    """
    ${vipsaniaEnv()}
    awk -v min=${minLen} '
        function flush() { if (name != "") { n++; if (length(seq) >= min) { used++; bp += length(seq); print ">" name > "train.fa"; print seq > "train.fa" } } }
        /^>/ { flush(); name = substr(\$1, 2); seq = ""; next }
        { seq = seq \$0 }
        END { flush(); printf "sequences\\tsequences_used\\tbp_used\\n%d\\t%d\\t%d\\n", n, used, bp > "train_set.tsv" }' ${genome}
    if [ ! -s train.fa ]; then
        echo "No sequence of ${genome} has ${minLen} bp or more: Vipsania cannot be trained for translation table ${table}." >&2
        exit 1
    fi
    printf '{"dataset": {"train_paths": ["%s"], "indexed_files": true}, "trainer": {"epochs": %d, "train_steps": %d}}\\n' \\
        "\$PWD/train.fa" ${epochs} ${steps} > train.json
    # a relative configuration is one of those of the image (/opt/vipsania/configs)
    config="${config}"
    case "\$config" in /*) ;; *) [ -f "\$config" ] || config="/opt/vipsania/\$config" ;; esac
    mkdir -p ckpt
    vipsania train "\$config" --translation_table ${table} --checkpoints ckpt -oc train.json > vipsania_train.log 2>&1 \\
        || { tail -n 30 vipsania_train.log >&2; exit 1; }
    run=""
    for d in \$(ls -1t ckpt); do
        if [ -s "ckpt/\$d/latest_checkpoint.weights.h5" ]; then run="\$d"; break; fi
    done
    if [ -z "\$run" ]; then
        echo "vipsania train wrote no checkpoint under ckpt/, see vipsania_train.log" >&2
        exit 1
    fi
    printf '%s\\n' "\$run" > MODEL_ID
    rm -f train.fa train.fa.fai
    """

    stub:
    """
    mkdir -p ckpt/stub_run
    touch ckpt/stub_run/latest_checkpoint.weights.h5 vipsania_train.log train_set.tsv
    echo stub_run > MODEL_ID
    """
}

// Annotate the genome with the model of TRAIN_VIPSANIA. No finetuning: the
// model was trained on this genome.
process RUN_VIPSANIA_TRAINED {
    label 'gpu', 'vipsania', 'bigmem'
    publishDir "${params.outdir}/intermediate/vipsania", mode: 'copy', overwrite: true, pattern: "*.log"

    input:
        path genome
        path ckpt
        path model_id
        val table

    output:
        path "vipsania.${genome.name}.gtf", emit: gtf
        path "vipsania.${genome.name}.log", emit: log, optional: true
        path "versions.tsv",               emit: versions

    script:
    def v = params.vipsania ?: [:]
    def extra = ''
    if (v.batch_size) extra += " -B ${v.batch_size}"
    if (v.context)    extra += " -T ${v.context}"
    if (v.extra_args) extra += " ${v.extra_args}"
    """
    ${vipsaniaEnv()}
    vipsania annotate "\$(cat ${model_id})" ${genome} \\
        --model_dir "\$PWD/${ckpt}" \\
        --translation_table ${table} \\
        -o "vipsania.${genome.name}.gtf"${extra}
    printf 'Vipsania\\t%s\\t%s\\n' "\$(python3 -c 'import importlib.metadata as m; print(m.version("vipsania"))' 2>/dev/null || echo unknown)" "${task.container ?: 'none'}" > versions.tsv
    """

    stub:
    """
    touch "vipsania.${genome.name}.gtf" "vipsania.${genome.name}.log" versions.tsv
    """
}
