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

## Process labels

| Label | Meaning |
| --- | --- |
| `container` | runs in the evidence tools image pinned in `base.config` |
| `tiberius`, `vipsania`, `drusilla`, `hint_rescue` | runs in that image instead (see [containers.md](containers.md)); never combined with `container` |
| `gpu` | needs a GPU; gets `containerOptions = '--nv'` |
| `bigmem` | high memory task, 100 GB by default |
| `download` | downloads of SRA reads, OrthoDB partitions and model weights; runs on the submitting host, which needs internet access (4 CPUs, 8 GB) |
| `local_only` | tiny task that runs on the submitting host |

Override resources per label (`withLabel:`) or per process (`withName:`) in
your config. Memory is a closure of `task.attempt` in `base.config`: a task
that the scheduler killed for memory or time (exit codes 137, 140, 143, 247)
is retried up to twice with twice and three times the memory; any other
error ends the run after the running tasks finish.
