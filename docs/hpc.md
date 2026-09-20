# Adapting the pipeline to an HPC

This page explains how to write a Nextflow config for your environment
(queues, scratch paths, GPU settings).

Shipped configs in [conf/](../conf):

| Config | Use |
| --- | --- |
| `base.config` | parameter defaults, process labels, container images. Always loaded by the launcher. |
| `local.config` | one machine, no scheduler |
| `slurm_generic.config` | starting point for SLURM clusters |
| `greifswald_hpc.config` | the Greifswald cluster |
| `user_hpc_template.config` | commented template for your own cluster |

Each of them includes `base.config`. `--nf_config` accepts a path, or the name
of a shipped config with or without the `.config` suffix (`slurm_generic`).

---

---

## 1. Create your personal HPC config

Copy the template shipped with the repository:

```bash
mkdir -p ~/nf_configs
cp conf/user_hpc_template.config ~/nf_configs/mycluster.config
```

Otherwise, use one of the example configs `local.config` for non HPC usage and `slurm_generic` a generic example for a slurm HPC.

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
For SLURM you need to set `envWhitelist = 'CUDA_VISIBLE_DEVICES'`
Make sure to include `containerOptions = '--nv'` in the GPU process section.

If you copy a shipped config out of `conf/`, replace its first line
`includeConfig 'base.config'` by the absolute path of `conf/base.config`, or
remove the line. The launcher always loads `base.config` before your config.

## Process labels

| Label | Meaning |
| --- | --- |
| `container` | runs in the tools image pinned in `base.config` |
| `gpu` | needs a GPU; gets `containerOptions = '--nv'` |
| `bigmem` | high memory task, 100 GB by default |
| `local_only` | tiny task that runs on the submitting host |

Override resources per label (`withLabel:`) or per process (`withName:`) in
your config.
