# Paludamentum

The paludamentum was the cloak worn by Roman emperors and commanders. This
repository is the cloak around the emperors of the
[Gaius-Augustus](https://github.com/Gaius-Augustus) gene finder family: a
Nextflow pipeline that prepares extrinsic evidence (proteins, RNA-Seq,
Iso-Seq), derives high-confidence genes from it, and integrates them with the
ab initio predictions of a deep learning gene finder.

The gene finders are git submodules of this repository. Paludamentum runs
them; they do not depend on Paludamentum.

| Submodule | Role | Status |
| --- | --- | --- |
| [Tiberius](https://github.com/Gaius-Augustus/Tiberius) (`tiberius/`) | gene finder | supported |
| [Vipsania](https://github.com/Gaius-Augustus/Vipsania) (`vipsania/`) | gene finder | supported, see [docs/vipsania.md](docs/vipsania.md) |
| [Drusilla](https://github.com/Gaius-Augustus/Drusilla) (`drusilla/`) | ORF annotator for assembled transcripts | high-confidence genes for vertebrate models, see [Drusilla flow](#drusilla-flow) |

> **Status (v0.4.0).** Paludamentum is the entry point: `paludamentum`
> launches the pipeline with Tiberius or Vipsania. Earlier versions were a
> submodule of Tiberius and were launched by `tiberius.py`; that direction is
> reversed since v0.3.0. v0.4.0 adds the [Drusilla flow](#drusilla-flow) for
> vertebrate models.

## Table of contents

- [What the pipeline does](#what-the-pipeline-does)
- [Installation](#installation)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Inputs and modes](#inputs-and-modes)
- [Outputs](#outputs)
- [Choosing the gene finder](#choosing-the-gene-finder)
- [Drusilla flow](#drusilla-flow)
- [Containers](#containers)
- [Running on an HPC](#running-on-an-hpc)
- [Repository layout](#repository-layout)
- [For maintainers](#for-maintainers)
- [Testing](#testing)
- [Roadmap](#roadmap)
- [License and citation](#license-and-citation)
- [Funding](#funding)

## What the pipeline does

1. **Ab initio prediction.** The genome is split into chunks, the gene finder
   runs on each chunk on a GPU, and the chunk predictions are merged.
2. **Protein evidence.** Proteins (your FASTA files and/or OrthoDB v12
   partitions) are aligned with miniprot, scored with
   miniprot-boundary-scorer, and converted to hints with miniprothint. For
   very large protein databases the ab initio proteins are used to select the
   most relevant source species first (DIAMOND).
3. **Transcript evidence.** Short reads are aligned with HISAT2, Iso-Seq reads
   with minimap2. Libraries with a low alignment rate are dropped. Transcripts
   are assembled with StringTie, and intron hints are extracted.
4. **High-confidence genes.** Assembled transcripts get ORFs from TransDecoder.
   ORFs that are supported by protein homology (DIAMOND) and by the scored
   protein alignments become the high-confidence (HC) gene set. TD2 can
   replace TransDecoder (`transdecoder: td2`); a comparison is in
   [docs/orf_finder_comparison.md](docs/orf_finder_comparison.md).
5. **Integration.** The HC genes are merged with the ab initio predictions into
   the final annotation, and its protein sequences are extracted.

For vertebrate gene finder models, runs with transcripts use the
[Drusilla flow](#drusilla-flow) in steps 4 and 5 instead.

Without any evidence input the pipeline runs step 1 only. That is useful to
parallelize a gene finder over several GPUs.

## Installation

```bash
git clone --recursive https://github.com/Gaius-Augustus/Paludamentum
cd Paludamentum
pip install -e .
```

`--recursive` checks out the submodules `tiberius/`, `vipsania/` and
`drusilla/`. In an existing clone run `git submodule update --init` after
`git pull`. `pip install -e .` installs the launcher (`paludamentum`, also
`python -m paludamentum`) linked to the checkout; the pipeline (`main.nf`,
`conf/`, `bin/`) runs from the checkout. If you install without `-e`, or
move the checkout, point the launcher at it with
`export PALUDAMENTUM_ROOT=/path/to/Paludamentum`.

The pipeline runs the gene finders and all tools in containers, so the
submodules do not have to be installed. The Tiberius checkout is used to
resolve model configuration names (`--model_cfg diatoms`) and, for runs
without containers, to find `tiberius.py`.

## Requirements

On the machine that launches the pipeline:

- Nextflow 25.04 or newer
- Java 17 or newer (required by Nextflow)
- Singularity or Apptainer (all tools run in containers, see [Containers](#containers))
- Python 3.9 or newer with `pyyaml`
- an NVIDIA GPU on the nodes that run the gene finder step

You do not need to install HISAT2, miniprot, StringTie and the other tools
when you use the containers. To run without containers, all tools must be in
your `PATH`, or be configured under `tools:` in the params file, and the gene
finder must be installed (`pip install ./vipsania`; Tiberius runs from the
checkout). The launcher option `--check_tools` verifies this.

## Quick start

```bash
# Tiberius, ab initio only, parallelized over GPUs by Nextflow
paludamentum --nf_config slurm_generic --genome genome.fa --model_cfg eudicotyledons

# Tiberius with evidence, inputs from a params file
paludamentum --params_yaml params.yaml --nf_config slurm_generic

# Vipsania with evidence, inputs on the command line
paludamentum --genefinder vipsania --nf_config local --genome genome.fa --model Fungi \
    --proteins proteins.faa --rnaseq_paired "/abs/path/rnaseq/*_{1,2}.fastq.gz"
```

`--nf_config` takes a path or the name of a config in `conf/`
(`local`, `slurm_generic`, ...); `local` sizes the tasks to the machine it
runs on, `slurm_generic` needs your GPU partition (see
[Running on an HPC](#running-on-an-hpc)). Evidence can be given on the command line
(`--proteins`, `--odb12Partitions`, `--rnaseq_paired`, `--rnaseq_single`,
`--isoseq`, `--rnaseq_sra_paired`, `--rnaseq_sra_single`, `--isoseq_sra`) or
in the params file; command line values override the file. Useful options:
`--dry_run` writes the params file and validates inputs and executables
without starting Nextflow, `--resume` continues a previous run, `--work_dir`
sets the Nextflow work directory. Arguments after `--` go to Nextflow.
`paludamentum --help` lists everything.

### Params file

Start from [conf/parameters.yaml](conf/parameters.yaml). Minimal example:

```yaml
genome: /abs/path/genome.fa
proteins: /abs/path/proteins.faa
rnaseq_paired: "/abs/path/rnaseq/*_{1,2}.fastq.gz"
outdir: results
tiberius:
  run: true
  model_cfg: eudicotyledons
```

The launcher expands `~` and environment variables and resolves relative
paths against the launch directory, not against the params file; it writes
the merged parameters of a run with absolute paths to `<outdir>/params.yaml`.
Paths must not contain spaces. All parameters are documented in
[docs/parameters.md](docs/parameters.md).

## Inputs and modes

| Parameter | Content |
| --- | --- |
| `genome` | genome FASTA, optionally gzipped (required) |
| `proteins` | one protein FASTA or a list; lists are concatenated |
| `odb12Partitions` | OrthoDB v12 partitions to download: `Metazoa`, `Vertebrata`, `Viridiplantae`, `Arthropoda`, `Fungi`, `Alveolata`, `Stramenopiles`, `Amoebozoa`, `Euglenozoa`, `Eukaryota` |
| `rnaseq_paired`, `rnaseq_single` | short read FASTQ files (glob or list) |
| `rnaseq_bam` | aligned short reads, for example from VARUS |
| `isoseq` | Iso-Seq FASTQ files |
| `rnaseq_sra_paired`, `rnaseq_sra_single`, `isoseq_sra` | SRA run accessions, downloaded by the pipeline |
| `min_alignment_rate` | libraries whose alignment rate (percent mapped) is below this value are dropped; default 80 |

The mode is inferred from the inputs and can be forced with `mode`
(`--mode`). A forced mode without its inputs, or transcript evidence without
proteins, is an error:

| Mode | Inputs |
| --- | --- |
| `abinitio` (historic name: `tiberius`) | genome only |
| `proteins` | proteins |
| `rnaseq` | proteins and short reads |
| `isoseq` | proteins and Iso-Seq |
| `mixed` | proteins, short reads and Iso-Seq |

Transcript evidence requires protein evidence, because HC genes need both.

## Outputs

All files are written to `outdir` (default `<genefinder>_results`). `<tool>`
is `tiberius` or `vipsania`.

| File | Content |
| --- | --- |
| `<tool>_evidence.gff3` | final annotation: ab initio predictions merged with HC genes |
| `<tool>_evidence_proteins.fa` | protein sequences of the final annotation |
| `<tool>_ab_initio.gff3` | ab initio predictions (in `intermediate/` when evidence is used) |
| `intermediate/hc.gff3` | high-confidence genes derived from the evidence (TransDecoder) |
| `intermediate/drusilla_orfs.gtf` | Drusilla ORFs, the HC genes of the Drusilla flow |
| `intermediate/<tool>_lgb_filtered.gtf` | ab initio predictions kept by the LightGBM filter (Drusilla flow) |
| `intermediate/<tool>_lgb_scores.tsv` | LightGBM class probabilities of all ab initio transcripts (Drusilla flow) |
| `intermediate/hint_rescue.gtf` | ab initio genes predicted again with protein hints (Drusilla flow; missing if the rescue was skipped) |
| `hintsfile.gff` | protein, RNA-Seq and Iso-Seq hints |
| `sra_downloads/` | reads downloaded from SRA |
| `params.yaml` | the merged parameters of this run, written by the launcher |

Nextflow's timeline, trace and report files are written as well.

## Choosing the gene finder

`--genefinder tiberius|vipsania` selects the gene finder. Without it, the
launcher takes `genefinder` from the params file, else the block with
`run: true`, else Vipsania if `--model` is given, else Tiberius. Only one
gene finder runs per pipeline run.

### Tiberius block

```yaml
tiberius:
  run: true
  model_cfg: eudicotyledons   # name in tiberius/model_cfg of the submodule, or a path
  model_dir: null             # directory with the extracted weights, for offline nodes
  result: null                # reuse an existing prediction instead of running
  min_split_size: 20000000    # minimal chunk size in bp
  max_files: 20               # maximal number of chunks
  max_parallel: null          # cap concurrent GPU tasks, e.g. 1 on a workstation
  batch_size: null
  seq_len: null
```

### Vipsania block

```yaml
vipsania:
  run: true
  model: Fungi                # clade name or model id
  model_dir: null             # directory with downloaded models, for offline nodes
  finetune: false             # see below
  finetune_epochs: null
  batch_size: null            # required on nodes without nvidia-smi
  context: null
  result: null
  max_parallel: null
```

The pipeline downloads the weights of the model configuration
(`weights_url` in the YAML) once at the start of the run, on the submitting
host like the other downloads, and stages them into every Tiberius task, so
the GPU nodes need no internet. The same holds for the Tiberius model of the
hint rescue. To skip the download, for example when the submitting host has
no internet either, download and extract the archive once and set
`model_dir` to the directory that holds the extracted weights directory: the
archive name without `.tar.gz`, for example `vertebrates_weights` or `fungi`.
Tiberius runs with `--model` on these weights and never downloads; a task
fails if the directory is missing or empty.

Vipsania finetuning is **off by default**. With `finetune: true` (or
`--finetune`) Vipsania first trains on the target genome and then annotates
it. Vipsania 1.0.1 cannot finetune without annotating, so in this case the
genome is processed in a single GPU task instead of chunks. All Vipsania
parameters, the model download and offline use are described in
[docs/vipsania.md](docs/vipsania.md).

## Drusilla flow

For vertebrate gene finder models, Paludamentum replaces the TransDecoder
high-confidence genes with Drusilla ORFs and filters the ab initio
predictions. The flow is derived from the one behind the Tiberius evidence
results on the GCB 2026 poster
([doi:10.13140/RG.2.2.24444.91521](https://doi.org/10.13140/RG.2.2.24444.91521)),
which used an earlier Drusilla model and had no stop/start codon fix, no hint
rescue and no merge of overlapping genes. With this flow and the Tiberius
`vertebrates` model, the gene-level F1 is 80.36 on *Takifugu rubripes* and
79.46 on *Bos taurus*. It runs when all of these hold:

- the run has transcripts (mode `rnaseq`, `isoseq` or `mixed`),
- the gene finder model is a Tiberius model whose `target_species` is
  `Vertebrata` or `Mammalia` (`vertebrates`, `mammalia*`), or Vipsania
  `Vertebrata` (`etb1go6q`); `drusilla.run: true` forces the flow for other
  vertebrate models,
- `drusilla.lgb_model` is not `null` (the default is the released model).

Steps:

1. One StringTie assembly of all reads (short reads, Iso-Seq with `-L`, or
   both with `--mix`).
2. Transcripts are kept if length >= 300, coverage >= 3 and TPM >= 1 (TPM >= 0.5
   for transcripts of 3000 nt or longer).
3. `drusilla annotate` predicts the ORFs of the kept transcripts, including
   ORFs truncated at the 3' or 5' end of a transcript. The released model
   (`vertebrates`: weights of the training run `cnn_lstm_vertebrates_run009`)
   is downloaded once on the submitting host and checked against the
   `weights_sha256` of its manifest in `drusilla/model_cfg/`, so the GPU nodes
   need no internet for it. With `cache_dir`, Drusilla takes the model from
   that cache instead (`drusilla models download vertebrates` with
   `DRUSILLA_CACHE_DIR` set to it fills it once).
4. Stop codon fix: ORFs with an early stop and truncated ORFs are extended to
   a stop codon supported by a miniprot alignment. Start codon fix: ORFs
   without an upstream in-frame stop are extended to a start codon hint of
   miniprothint. Isoforms whose CDS is a subsequence of another one are then
   collapsed.
5. A LightGBM model classifies each ab initio transcript as wrong, partial or
   correct, from miniprot alignment and miniprothint hint features. Transcripts
   with P(partial) + P(correct) >= 0.5 whose most likely class is `correct`
   are kept. The model (`drusilla_lgb_3class_v1`, a LightGBM text model) is
   downloaded by Nextflow at the start of the run and checked against its
   sha256, so the compute nodes need no internet for it.
6. Hint rescue: loci of `partial` transcripts without a kept transcript are
   predicted again by Tiberius with the hints of the best protein chain of the
   locus; loci without a protein chain are predicted ab initio, without hints.
   With `rescue_orf_filter: true`, loci where a Drusilla ORF already has all
   introns of the chain are skipped; the filter is off by default, and the
   F1 values above were measured without it. This needs a Tiberius with
   `--hints` (branch `hint_integration`, with bricks2marble from its branch
   `intron_hints`). The step therefore runs in its own image (see
   [Containers](#containers)). `rescue_tiberius`
   is another `tiberius.py` to use instead. If the Tiberius has no `--hints`,
   the step is skipped with a warning in the Nextflow log, and no GPU task is
   started.
7. The kept and rescued ab initio transcripts and the Drusilla ORFs are merged
   into `<tool>_evidence.gff3`.

```yaml
drusilla:
  run: auto              # auto: vertebrate models only; true: always; false: never
  model: vertebrates     # released Drusilla model
  weights: null          # a local .weights.h5 instead (needs config)
  config: null
  cache_dir: null        # model cache (DRUSILLA_CACHE_DIR); null = downloaded once on the submitting host
  batch_size: null       # null = Drusilla sizes the batch from the GPU memory
  shards: 1              # parallel Drusilla processes; e.g. 24 on a 48-core CPU node
  min_coding_length: 200 # minimal CDS length of a Drusilla ORF
  lgb_model: https://bioinf.uni-greifswald.de/bioinf/drusilla/models/drusilla_lgb_3class_v1.tar.gz
                         # LightGBM model of the ab initio filter (required): the
                         # archive (URL or file), its unpacked directory, or a .txt with its .json
  lgb_model_sha256: d5bf3a9774914c6ca6421507c4cb956ef813b0bd2305bf45717af5f5eb42f1da
                         # sha256 of the archive; null for another model
  lgb_threshold: 0.5
  lgb_keep: correct      # classes kept after the threshold, e.g. "correct,partial"
  fix_stop: true         # stop codon fix of the ORFs
  fix_start: true        # start codon fix of the ORFs (needs fix_stop)
  rescue: true           # hint rescue of partial ab initio genes
  rescue_tiberius: null  # tiberius.py with --hints; null = the one of the hint rescue image
  rescue_model_cfg: null # model file or name; null = the Tiberius model of the run, vertebrates for Vipsania
  rescue_flank: 25000
  rescue_hint_weight: 2.5
  rescue_orf_filter: false # skip loci where a Drusilla ORF has all introns of the chain
  rescue_seq_len: 99990   # Tiberius seq_len of the rescue loci
  min_length: 300        # StringTie pre-filter
  min_cov: 3
  min_tpm: 1
  long_length: 3000
  min_tpm_long: 0.5
```

Drusilla runs on a GPU (label `gpu`) in the Drusilla image, which also runs
the LightGBM filter (`lightgbm` 4.7.0, `pyfaidx`, `pandas`). Without a
GPU, send `DRUSILLA_ANNOTATE` to CPU nodes in your site config and set
`shards`: one Drusilla process uses only one to two cores, so `shards`
splits the transcripts by gene into parts that run in parallel. On Bos taurus
24 shards used 26 of 48 cores and 40 GB.

## Containers

| Processes | Image | Built from |
| --- | --- | --- |
| evidence tools (StringTie, HISAT2, minimap2, miniprot, DIAMOND, TransDecoder, ...) | `docker://gaiusaugustus/paludamentum-evidence:0.1.0` | [docker/evidence/Dockerfile](docker/evidence/Dockerfile) |
| Tiberius | `docker://gaiusaugustus/tiberius:<version>` | `Dockerfile` in the Tiberius repository |
| Vipsania | `docker://gaiusaugustus/vipsania:<version>` | `Dockerfile` in the Vipsania repository |
| Drusilla flow (ORFs, LightGBM filter) | `docker://gaiusaugustus/drusilla:<version>` | [docker/drusilla/Dockerfile](docker/drusilla/Dockerfile) |
| hint rescue of the Drusilla flow | `docker://gaiusaugustus/paludamentum-hint-rescue:0.2.0` | [docker/hint_rescue/Dockerfile](docker/hint_rescue/Dockerfile) |

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

The images are pinned in [conf/base.config](conf/base.config) through the
process labels `container` (evidence image), `tiberius`, `vipsania`, `drusilla`
and `hint_rescue`; a process carries at most one of them (processes that run
on the submitting host carry none). The image
tag of a gene finder must match the version of its submodule; the launcher
warns when they differ. The pipeline
scripts in `bin/` are not part of an image. Nextflow adds `bin/` to the
`PATH` of every task and mounts it into the container.

Images are pulled once into `~/.cache/paludamentum/singularity` and shared by
all runs. Set `NXF_SINGULARITY_CACHEDIR`, or `singularity.cacheDir` in your
config, to use another directory (the Tiberius image is about 11 GB).

Vipsania requires `tensorflow<2.20` and therefore does not support Blackwell
GPUs. See the Vipsania container documentation.

## Running on an HPC

Copy [conf/user_hpc_template.config](conf/user_hpc_template.config), set your
queues, GPU options and scratch paths, and pass it with `--nf_config`.
`conf/slurm_generic.config` is a starting point for SLURM,
`conf/local.config` runs everything on one machine. Every user config includes
`conf/base.config`. `--nf_config` accepts a path or the name of a config in
`conf/`. Details are in [docs/hpc.md](docs/hpc.md).

GPU processes carry the label `gpu`. Give them your GPU queue and
`containerOptions = '--nv'`. On SLURM keep
`singularity.envWhitelist = 'CUDA_VISIBLE_DEVICES'`. Downloads (SRA reads,
OrthoDB partitions, Tiberius weights, Vipsania models) carry the label `download` or
`local_only` and run on the submitting host, which needs internet access.
Tasks that the scheduler kills for memory or time are retried twice with more
memory; any other error stops the run.

## Repository layout

```text
main.nf                entry workflow
nextflow.config        manifest, includes conf/base.config
lib_nf/                shared Groovy functions
modules/               Nextflow processes
subworkflows/          inputs, protein, RNA-Seq, Iso-Seq, HC genes, gene finder
bin/                   scripts called by processes
conf/                  base config, site configs, parameters.yaml, blosum62.csv
paludamentum/          Python launcher (paludamentum, python -m paludamentum)
tiberius/              submodule: Tiberius (gene finder)
vipsania/              submodule: Vipsania (gene finder)
drusilla/              submodule: Drusilla (ORF annotator for transcripts)
docs/                  parameters.md, hpc.md, vipsania.md, orf_finder_comparison.md
tests/                 launcher tests and Nextflow stub runs
```

## For maintainers

- **Submodule pinning.** Each submodule is pinned to the release whose
  container image `conf/base.config` runs. A new gene finder release means:
  bump the submodule (`git -C tiberius checkout <tag>`, `git add tiberius`)
  and the image tag in `conf/base.config` in one commit.
  `tests/test_launcher.py` checks that the two agree.
- **Stable interfaces.** Users rely on the published file names, on the
  `tiberius.*` and `vipsania.*` parameter blocks, and on `conf/<name>.config`.
  Do not change them without a deprecation path.
- **Adding a gene finder.** Add the repository as a submodule, a module with
  the run process, a label with its container in `conf/base.config`, a
  parameter block, a branch in `subworkflows/genefinder.nf`, and entries in
  the launcher's `GENEFINDER_CLI`, `SUBMODULES` and `GENEFINDER_CLI_KEYS`
  tables. The process takes a genome FASTA and emits GTF or GFF3.
  `bin/merge_annotations.py` renumbers gene IDs during merging, writes
  transcripts with a CDS as `mRNA` and marks their genes
  `gene_biotype=protein_coding`, as in NCBI/Ensembl GFF3, so that
  [Annotrieve](https://genome.crg.es/annotrieve/) reports them the same way.
- Renaming a process invalidates `-resume` for runs in progress. Mention it in
  the release notes.

## Testing

```bash
pip install -e .[test]
pytest tests --ignore=tests/test_stub_run.py   # launcher and scripts, no Nextflow needed
pytest tests/test_stub_run.py                  # needs nextflow, or NEXTFLOW_BIN=/path/to/nextflow
nextflow lint main.nf modules subworkflows lib_nf
```

CI runs the same on Python 3.9 and 3.12, and the stub runs on the oldest
supported Nextflow (25.04.0) and the latest stable release.

The stub runs execute `nextflow run main.nf -stub-run -c tests/stub.config` for
every mode on the tiny inputs in `tests/data`. They check the wiring and the
published file names without tools, containers or a GPU. Every process has a
`stub:` block for this purpose; keep it in sync when you change outputs. Tests
that need the Tiberius submodule are skipped when it is not checked out. Real
smoke tests use `tiberius/test_data/Panthera_pardus` and
`vipsania/docs/example/aspergillus_fumigatus_chr7.fa`.

## Roadmap

- [x] v0.1.0: copy of the pipeline from Tiberius, launcher, stub blocks, tests, CI
- [x] v0.2.0: gene finder abstraction and Vipsania processes
- [x] v0.3.0: Paludamentum imports the gene finders (submodules `tiberius/`,
      `vipsania/`, `drusilla/`) and is launched by `paludamentum`; Tiberius
      no longer runs the pipeline
- [x] v0.4.0: Drusilla flow for vertebrate models (Drusilla ORFs as HC genes,
      LightGBM filter of the ab initio predictions)
- [ ] `vipsania annotate --finetune_only`, so finetuning can be combined with chunked annotation

## Known issues

Drusilla flow:

- The LightGBM model was trained on Drusilla ORFs of 48 vertebrates (listed
  in the `.json` of the archive) and is applied to gene finder predictions.
- The hint rescue needs the Tiberius branch `hint_integration` and the
  bricks2marble branch `intron_hints`, which are not part of a release. It
  runs in a separate image until they are. With a site config that runs the
  rescue in the released Tiberius image, it is skipped with a warning.
- The scripts of the flow (`bin/filter_stringtie_gtf.py`,
  `compute_orf_features.py`, `apply_lgb_model_gtf.py`, `fix_stop_by_miniprot.py`,
  `prepare_hint_rescue_loci.py`, `filter_and_merge_rescue_gtf.py`,
  `chainedHints.py`) are copies from tiberius_orf_finder and Tiberius.

## License and citation

Artistic License 1.0, see [LICENSE](LICENSE). Scripts in `bin/` copied from
Tiberius say so in their header and stay under its MIT License, see
[LICENSE-Tiberius](LICENSE-Tiberius).

If you use the pipeline, cite the gene finder you ran and this repository
([CITATION.cff](CITATION.cff)). The references are in
the READMEs of [Tiberius](https://github.com/Gaius-Augustus/Tiberius) and
[Vipsania](https://github.com/Gaius-Augustus/Vipsania). Please also cite the
tools that the pipeline runs on your data: miniprot, miniprot-boundary-scorer,
miniprothint, HISAT2, minimap2, StringTie, TransDecoder, DIAMOND, SAMtools,
BEDTools, gffread, and AUGUSTUS (bam2hints).

If the [Drusilla flow](#drusilla-flow) ran, please also cite:

- the poster the flow and its LightGBM filter are derived from. The hint
  rescue is not on the poster and is not published yet; it is Lars Gabriel's
  work (tiberius_orf_finder scripts, Tiberius branch `hint_integration`).
  Gabriel L, Hoff KJ. Annotating Eukaryotic Genomes by Combining Deep Learning
  with extrinsic evidence. Poster, GCB 2026.
  [doi:10.13140/RG.2.2.24444.91521](https://doi.org/10.13140/RG.2.2.24444.91521)
- [Drusilla](https://github.com/Gaius-Augustus/Drusilla) (repository; there is
  no publication yet)
- LightGBM: Ke G, Meng Q, Finley T, Wang T, Chen W, Ma W, Ye Q, Liu T-Y.
  LightGBM: A Highly Efficient Gradient Boosting Decision Tree. Advances in
  Neural Information Processing Systems 30 (NIPS 2017).
  [proceedings](https://proceedings.neurips.cc/paper_files/paper/2017/hash/6449f44a102fde848669bdd9eb6b76fa-Abstract.html)
- Tiberius, also when Vipsania is the gene finder, because the hint rescue
  runs Tiberius

## Funding

Funded by the Deutsche Forschungsgemeinschaft (DFG, German Research
Foundation), project number 552910312: "AI-GUSTUS: Eine cloud-native Pipeline
für genaue Genom-Annotation" (AI-GUSTUS: a cloud-native pipeline for accurate
genome annotation).
