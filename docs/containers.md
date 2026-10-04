# Containers

| Processes | Image | Built from |
| --- | --- | --- |
| evidence tools (StringTie, HISAT2, minimap2, miniprot, DIAMOND, TransDecoder, ...) | `docker://gaiusaugustus/paludamentum-evidence:0.1.0` | [docker/evidence/Dockerfile](../docker/evidence/Dockerfile) |
| Tiberius | `docker://gaiusaugustus/tiberius:<version>` | `Dockerfile` in the Tiberius repository |
| Vipsania | `docker://gaiusaugustus/vipsania:<version>` | `Dockerfile` in the Vipsania repository |
| Drusilla flow (ORFs, LightGBM filter) | `docker://gaiusaugustus/drusilla:<version>` | [docker/drusilla/Dockerfile](../docker/drusilla/Dockerfile) |
| hint rescue of the Drusilla flow | `docker://gaiusaugustus/paludamentum-hint-rescue:0.2.0` | [docker/hint_rescue/Dockerfile](../docker/hint_rescue/Dockerfile) |

The evidence image is Ubuntu 24.04 with every tool pinned: Ubuntu packages
at fixed versions, release archives checked against their SHA-256, and tools
built from source at fixed commits. The versions are those of the image the
pipeline was benchmarked with. Its tag is versioned on its own; bump it
whenever the Dockerfile changes.

The Drusilla image is Drusilla at the commit the `drusilla/` submodule
points to, on the same NGC TensorFlow base as the Tiberius image, with
`lightgbm` 4.7.0 (the version the released LightGBM model was converted
with). Its tag is the Drusilla version.

The hint rescue image is the Tiberius image of the submodule version with the
Tiberius branch `hint_integration` and the bricks2marble branch `intron_hints`
in place of the released versions, plus samtools. It is used for the hint rescue only and will be
dropped when both branches are released.

The images are pinned in [conf/base.config](../conf/base.config) through the
process labels `container` (evidence image), `tiberius`, `vipsania`, `drusilla`
and `hint_rescue`; a process carries at most one of them (processes that run
on the submitting host carry none). The image
tag of a gene finder must match the version of its submodule; the launcher
warns when they differ. The pipeline
scripts in `bin/` are not part of an image. Nextflow adds `bin/` to the
`PATH` of every task and mounts it into the container.

Images are pulled once into `~/.cache/paludamentum/singularity` and shared by
all runs. Set `NXF_SINGULARITY_CACHEDIR`, or `singularity.cacheDir` in your
config, to use another directory. Compressed sizes on Docker Hub: Tiberius
9.4 GB, Vipsania 9.4 GB, Drusilla 9.4 GB, hint rescue 9.6 GB, evidence tools
0.25 GB. Most of a GPU image is its NGC TensorFlow base, so a run of the
Drusilla flow pulls about 29 GB (Tiberius, Drusilla, hint rescue, evidence).

The Tiberius, Vipsania and Drusilla images are built on NVIDIA's NGC
TensorFlow image 25.02 (TensorFlow 2.17, CUDA 12.8) and run on NVIDIA GPUs up
to and including the Blackwell generation; Tiberius and Vipsania were run on
RTX PRO 6000 Blackwell cards. A `pip install` of Vipsania outside the
container does not support Blackwell GPUs. See the Vipsania container
documentation.
