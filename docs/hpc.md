# Adapting the pipeline to an HPC

This page explains how to write a Nextflow config for your environment
(queues, scratch paths, GPU settings).

Shipped configs in [conf/](../conf):

| Config | Use |
| --- | --- |
| `base.config` | parameter defaults, process labels, container images. Always loaded by the launcher. |
| `local.config` | one machine, no scheduler; tasks are sized to the machine |
| `slurm_generic.config` | starting point for SLURM clusters; set your GPU partition |
| `user_hpc_template.config` | commented template for your own cluster |

The launcher loads `base.config` before your config. `--nf_config` accepts a
path, or the name of a shipped config with or without the `.config` suffix
(`slurm_generic`).

## 1. Create your personal HPC config

Copy the template shipped with the repository to any place:

```bash
mkdir -p ~/nf_configs
cp conf/user_hpc_template.config ~/nf_configs/mycluster.config
```

and pass it with `--nf_config ~/nf_configs/mycluster.config`. The copy needs
no `includeConfig` line: the launcher always loads `conf/base.config` first.
If you copy `local.config` or `slurm_generic.config` instead, remove their
first line `includeConfig 'base.config'`, which only resolves inside `conf/`.

For a machine without a scheduler use `local.config` as it is.

## 2. Edit the config for your cluster

Open `mycluster.config` and adjust:

### Queues / partitions

```bash
process {
  executor = 'slurm'
  queue    = 'batch'  // change to your default queue
  time     = '12h'
  cpus     = 8
  memory   = '32 GB'

  withLabel: gpu {
    queue          = 'gpu'          // your GPU queue
    cpus           = 16
    memory         = '64 GB'
    time           = '24h'
    clusterOptions = '--gpus=1'     // or '--gres=gpu:1' etc.
    containerOptions = '--nv'
  }
}
```

Set the correct queue/partition names and any required `clusterOptions`
(e.g. constraints, accounts, QoS, etc.) for non GPU tasks like Miniprot and for GPU tasks for the gene finder.
Make sure to include `containerOptions = '--nv'` for GPU processes if you want to use the Singularity container.


###  Singularity

If you want to use the Singularity container for the dependencies:
```
singularity {
  enabled    = true
  autoMounts = true
  // envWhitelist = 'CUDA_VISIBLE_DEVICES,NVIDIA_VISIBLE_DEVICES'
}
```
`base.config` already whitelists `CUDA_VISIBLE_DEVICES`, which SLURM needs;
add more variables here if your site requires them. Make sure to include
`containerOptions = '--nv'` in the GPU process section.

Images are pulled once into `singularity.cacheDir`
(`~/.cache/paludamentum/singularity` by default, or `NXF_SINGULARITY_CACHEDIR`).
On a cluster put it on a shared file system that the compute nodes can read.

### Scratch

Both shipped SLURM configs set `scratch = true`: every task runs in a
temporary directory on the compute node (`$TMPDIR`, or a directory made by
`mktemp` when it is unset) and copies only its declared outputs back to the
work directory. Temporary files, such as the sort chunks and intermediate BAMs
of the read mapping, then stay off the shared file system. The nodes need
enough local `$TMPDIR` for the largest task; set `scratch = false` in your
config if they have little local disk. FEELnc tasks always run in their work
directory (`scratch = false` for the label `feelnc` in `base.config`): its
container cannot start in a task directory under `/tmp`.

### Work directory

Nextflow does not delete the task directories in `work/` (or `--work_dir`)
after a run. All final outputs are copied to `outdir`, so the work
directory of a finished run only serves `-resume`, and it grows to hundreds of
GB on a large genome. `paludamentum --cleanup` adds `conf/cleanup.config`
(Nextflow's `cleanup = true`): once the run has completed successfully, its
task directories are deleted. A failed run keeps them for `--resume`; a
cleaned run cannot be resumed. Without the launcher, pass
`-c conf/cleanup.config` to `nextflow run`. Nextflow deletes nothing during
the run, so the work directory still needs the peak space of the whole run.

## Process labels

| Label | Meaning |
| --- | --- |
| `container` | runs in the evidence tools image pinned in `base.config` |
| `tiberius`, `vipsania`, `drusilla`, `hint_rescue` | runs in that image instead (see [containers.md](containers.md)); never combined with `container` |
| `postprocess`, `busco`, `feelnc`, `omark`, `fantasia` | post-processing images: compleasm, gffcompare, tRNAscan-SE, Infernal, pybarrnap, the report; BUSCO; FEELnc; OMArk; FANTASIA-Lite (see [postprocessing.md](postprocessing.md)) |
| `gpu` | needs a GPU; gets `containerOptions = '--nv'` (FANTASIA_ANNOTATE: `--nv` and a bind, set by `withName` in `base.config`) |
| `bigmem` | high memory task, 100 GB by default (also compleasm and BUSCO on the genome, cmscan) |
| `download` | downloads of SRA reads, OrthoDB partitions, model weights, BUSCO lineages, Rfam and the NCBI taxonomy; runs on the submitting host, which needs internet access (4 CPUs, or all of the host's CPUs when it has fewer; 8 GB). The BUSCO and Rfam downloads also carry the image label `postprocess`, the NCBI taxonomy (built with ete3 for OMArk) the image label `omark`. |
| `local_only` | tiny task that runs on the submitting host |

Tasks get `params.threads` CPUs (48) unless the site config sets
`process.cpus`, with these exceptions in `base.config`: the single-threaded
post-processing tasks (sanity filter, UTRs, final files, statistics, report,
summaries, GFF3 conversions) and the NCBI taxonomy download get one CPU,
FANTASIA_ANNOTATE gets 8 and the GPU probe FANTASIA_GPU_CHECK one CPU, 2 GB
and 1 h (so it is scheduled quickly at the start of the run; keep it on the
same queue as FANTASIA_ANNOTATE). These are `withName` settings, which beat
`process.cpus` and `withLabel:` of a site config; change them with a
`withName:` of your own.

Override resources per label (`withLabel:`) or per process (`withName:`) in
your config. Memory is a closure of `task.attempt` in `base.config`: a task
that the scheduler killed for memory or time (exit codes 137, 140, 143, 247)
is retried up to twice with twice and three times the memory; any other
error ends the run after the running tasks finish.
