# Vipsania as gene finder

[Vipsania](https://github.com/Gaius-Augustus/Vipsania) is an unsupervised
*ab initio* gene finder. It takes no extrinsic evidence itself. In Paludamentum
its predictions play the same role as those of Tiberius: they are merged with
the high-confidence genes that the pipeline derives from proteins and
transcripts.

## Selecting Vipsania

Set `vipsania.run: true` in the params file, or `genefinder: vipsania`.
Only one gene finder can run per pipeline run.

```yaml
genome: /abs/path/genome.fa
proteins: /abs/path/proteins.faa
rnaseq_paired: "/abs/path/rnaseq/*_{1,2}.fastq.gz"
outdir: results
vipsania:
  run: true
  model: Fungi
```

```bash
paludamentum --params_yaml params.yaml --nf_config slurm_generic
```

or entirely from the command line:

```bash
paludamentum --genefinder vipsania --nf_config slurm_generic --genome genome.fa --model Fungi \
    --proteins proteins.faa --rnaseq_paired "/abs/path/rnaseq/*_{1,2}.fastq.gz"
```

`--genefinder vipsania` can be omitted when `--model` is given and
`--model_cfg` is not. Vipsania itself is a git submodule of this repository
(`vipsania/`); the `vipsania` command does not run the pipeline.

## Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `vipsania.run` | `false` | Run Vipsania and merge its predictions with the HC genes. |
| `vipsania.model` | none (required) | Clade name (`Fungi`, `Insecta`, `Vertebrata`, ...) or model id. See the model table in the Vipsania README. |
| `vipsania.model_dir` | none | Directory with models from `vipsania download <model> -d DIR`. No download is attempted then. |
| `vipsania.result` | none | Existing Vipsania prediction (GTF/GFF3) to use instead of running Vipsania. |
| `vipsania.finetune` | `false` | Finetune the model on the target genome before annotating. |
| `vipsania.finetune_epochs`, `finetune_B`, `finetune_lr` | Vipsania defaults | Forwarded to `vipsania annotate`. |
| `vipsania.batch_size` | automatic | `-B`. Automatic sizing needs `nvidia-smi`; set it on CPU-only nodes. |
| `vipsania.context` | Vipsania default | `-T`, genome context length. |
| `vipsania.min_split_size` | `20000000` | Minimal size in bp of a genome chunk. |
| `vipsania.max_files` | `20` | Maximal number of genome chunks. |
| `vipsania.max_parallel` | unlimited | Cap of concurrently running Vipsania tasks. |
| `vipsania.extra_args` | none | Appended verbatim to `vipsania annotate`. |

## Finetuning

Finetuning is off by default. Vipsania recommends it, in particular when the
quality of the repeat masking is uncertain. With `vipsania.finetune: true`
Vipsania trains on the FASTA file that it annotates. Vipsania 1.0.1 cannot
finetune without annotating, so the pipeline then runs **one** Vipsania task on
the whole genome instead of one task per chunk. The checkpoint
(`finetuning_*`) and the log are published in `intermediate/vipsania/`.

Without finetuning the genome is split into chunks that are annotated in
parallel, as for Tiberius.

## Models and offline nodes

Without `vipsania.model_dir` the process `DOWNLOAD_VIPSANIA_MODEL` fetches the
model. It carries the label `local_only` and therefore runs on the submitting
host, which needs internet access. For clusters without it:

```bash
vipsania download Fungi -d /abs/path/vipsania_models
```

and set `vipsania.model_dir: /abs/path/vipsania_models`.

## Container

Vipsania is not part of the Tiberius or the evidence image. The
Vipsania processes carry the label `vipsania`, to which `conf/base.config`
assigns `docker://gaiusaugustus/vipsania:1.0.1`. They must never carry the
label `container`, because a config selector overrides a `container` directive
in the process.

Singularity and Apptainer do not run the ENTRYPOINT of the image, which points
TensorFlow to its bundled CUDA libraries. The processes set `LD_LIBRARY_PATH`
themselves. If Vipsania reports that no GPU was found, check that the process
has `containerOptions = '--nv'` (label `gpu`).

The image is built on NVIDIA's NGC TensorFlow image and supports GPUs up to
and including the Blackwell generation (for example RTX PRO 6000). A
`pip install` of Vipsania outside the container does not support Blackwell
GPUs.

## Notes

- Vipsania drops sequences shorter than 1000 bp.
- Soft masking is read as repeat information. See the Vipsania documentation
  on masking and finetuning.
- Output files are named `vipsania_evidence.gff3`,
  `vipsania_evidence_proteins.fa` and `vipsania_ab_initio.gff3`.
