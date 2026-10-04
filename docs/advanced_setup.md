# Advanced setup

The supported way to run Paludamentum is with Singularity or Apptainer, as
described in the [README](../README.md). This page is for experts who need
another setup.

## Running without containers

Without containers you install every tool and the gene finder yourself.

1. Install the tools of the steps that your run uses (HISAT2, minimap2,
   SAMtools, StringTie, miniprot, miniprot-boundary-scorer, miniprothint,
   DIAMOND, TransDecoder, AUGUSTUS `bam2hints`, GffRead, SRA Toolkit). They
   must be in your `PATH`, or be configured under `tools:` in the params
   file, see [parameters.md](parameters.md#parameters). The versions that the
   pipeline was benchmarked with are in
   [docker/evidence/Dockerfile](../docker/evidence/Dockerfile).
2. Install the gene finder. Vipsania: `pip install ./vipsania`. Tiberius runs
   from the checkout of the submodule (`tiberius/tiberius.py`) and needs its
   Python dependencies. A `pip install` of Vipsania does not support
   Blackwell GPUs; the container image does.
3. Switch the containers off in your Nextflow config (see [hpc.md](hpc.md)):

   ```groovy
   singularity {
     enabled = false
   }
   ```

4. Launch with `--skip_singularity_check`, so that the launcher does not
   require a `singularity` or `apptainer` executable, and with
   `--check_tools`, which verifies that all tools of the run are found:

   ```bash
   paludamentum --params_yaml params.yaml --nf_config mycluster.config \
       --skip_singularity_check --check_tools --dry_run
   ```

The [Drusilla flow](drusilla_flow.md) additionally needs Drusilla, LightGBM
and, for the hint rescue, the Tiberius and bricks2marble branches of the hint
rescue image, see [containers.md](containers.md).

## Location of the checkout

`pip install -e .` installs the launcher (`paludamentum`, also
`python -m paludamentum`) linked to the checkout; the pipeline (`main.nf`,
`conf/`, `bin/`) runs from the checkout. If you install without `-e`, or move
the checkout, point the launcher at it:

```bash
export PALUDAMENTUM_ROOT=/path/to/Paludamentum
```

The gene finder submodules do not have to be installed when the pipeline
runs with containers. The Tiberius checkout is used to resolve model
configuration names (`--model_cfg diatoms`) and, for runs without containers,
to find `tiberius.py`.

## Image cache and own images

Images are pulled once into `~/.cache/paludamentum/singularity`
(`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir` change the directory).
To use images that you pulled or built yourself, set the `container` of the
process labels in your config, see [containers.md](containers.md) and
[hpc.md](hpc.md).
