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
| [Drusilla](https://github.com/Gaius-Augustus/Drusilla) (`drusilla/`) | ORF annotator for assembled transcripts | imported; use in the high-confidence gene step is planned, see [Roadmap](#roadmap) |

> **Status (v0.3.0).** Paludamentum is the entry point: `paludamentum`
> launches the pipeline with Tiberius or Vipsania. Earlier versions were a
> submodule of Tiberius and were launched by `tiberius.py`; that direction is
> reversed since v0.3.0.

## Table of contents

- [What the pipeline does](#what-the-pipeline-does)
- [Installation](#installation)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Inputs and modes](#inputs-and-modes)
- [Outputs](#outputs)
- [Choosing the gene finder](#choosing-the-gene-finder)
- [Containers](#containers)
- [Running on an HPC](#running-on-an-hpc)
- [Repository layout](#repository-layout)
- [For maintainers](#for-maintainers)
- [Testing](#testing)
- [Roadmap](#roadmap)
- [Known issues](#known-issues)
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
   protein alignments become the high-confidence (HC) gene set.
5. **Integration.** The HC genes are merged with the ab initio predictions into
   the final annotation, and its protein sequences are extracted.

Without any evidence input the pipeline runs step 1 only. That is useful to
parallelize a gene finder over several GPUs.

## Installation

```bash
git clone --recursive https://github.com/Gaius-Augustus/Paludamentum
cd Paludamentum
pip install .
```

`--recursive` checks out the submodules `tiberius/`, `vipsania/` and
`drusilla/`. In an existing clone run `git submodule update --init` after
`git pull`. `pip install .` installs the launcher (`paludamentum`, also
`python -m paludamentum`); the pipeline runs from the checkout.

The pipeline runs the gene finders and all tools in containers, so the
submodules do not have to be installed. The Tiberius checkout is used to
resolve model configuration names (`--model_cfg diatoms`) and, for runs
without containers, to find `tiberius.py`.

## Requirements

On the machine that launches the pipeline:

- Nextflow
- Java 11 or newer
- Singularity or Apptainer (all tools run in containers, see [Containers](#containers))
- Python 3 with `pyyaml`
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
(`local`, `slurm_generic`, ...). Evidence can be given on the command line
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

Use absolute paths. Nextflow does not expand `~`, and relative paths are
resolved against the launch directory, not the location of the params file.
All parameters are documented in [docs/parameters.md](docs/parameters.md).
The launcher writes the merged parameters of a run to `<outdir>/params.yaml`.

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

The mode is inferred from the inputs and can be forced with `mode`:

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
| `intermediate/hc.gff3` | high-confidence genes derived from the evidence |
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

Vipsania finetuning is **off by default**. With `finetune: true` (or
`--finetune`) Vipsania first trains on the target genome and then annotates
it. Vipsania 1.0.0 cannot finetune without annotating, so in this case the
genome is processed in a single GPU task instead of chunks. All Vipsania
parameters, the model download and offline use are described in
[docs/vipsania.md](docs/vipsania.md).

## Containers

Paludamentum does not ship its own image.

| Processes | Image | Built from |
| --- | --- | --- |
| all evidence tools and Tiberius | `docker://larsgabriel23/tiberius:<version>` | `Dockerfile` in the Tiberius repository |
| Vipsania | `docker://gaiusaugustus/vipsania:<version>` | `Dockerfile` in the Vipsania repository |

The images are pinned in [conf/base.config](conf/base.config) through the
process labels `container` and `vipsania`. The image tag must match the
version of the submodule; the launcher warns when they differ. The pipeline
scripts in `bin/` are not part of an image. Nextflow adds `bin/` to the
`PATH` of every task and mounts it into the container.

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
`singularity.envWhitelist = 'CUDA_VISIBLE_DEVICES'`.

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
docs/                  parameters.md, hpc.md, vipsania.md
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
pytest tests/test_launcher.py      # launcher, no Nextflow needed
pytest tests/test_stub_run.py      # needs nextflow, or NEXTFLOW_BIN=/path/to/nextflow
nextflow lint main.nf modules subworkflows
```

The stub runs execute `nextflow run main.nf -stub-run -c tests/stub.config` for
every mode on the tiny inputs in `tests/data`. They check the wiring and the
published file names without tools, containers or a GPU. Every process has a
`stub:` block for this purpose; keep it in sync when you change outputs. Tests
that need the Tiberius submodule are skipped when it is not checked out. Real
smoke tests use `tiberius/test_data/Panthera_pardus` and
`vipsania/docs/example/aspergillus_fumigatus_chr7.fa`.

## Roadmap

The full plan is in [MIGRATION_PLAN.md](MIGRATION_PLAN.md).

- [x] v0.1.0: copy of the pipeline from Tiberius, launcher, stub blocks, tests, CI
- [x] v0.2.0: gene finder abstraction and Vipsania processes
- [x] v0.3.0: Paludamentum imports the gene finders (submodules `tiberius/`,
      `vipsania/`, `drusilla/`) and is launched by `paludamentum`; Tiberius
      no longer runs the pipeline
- [ ] Drusilla as alternative to TransDecoder in the HC gene step
- [ ] `vipsania annotate --finetune_only`, so finetuning can be combined with chunked annotation

## Known issues

Carried over unchanged from Tiberius and tracked for a later fix:

- The automatic mode inference has operator precedence slips. Set `mode`
  explicitly if the inferred mode is not what you expect.
- `restart` and `prothint_conflict_filter` are declared but unused.

## License and citation

Artistic License 1.0, see [LICENSE](LICENSE). Scripts in `bin/` copied from
Tiberius say so in their header and stay under its MIT License, see
[LICENSE-Tiberius](LICENSE-Tiberius).

If you use the pipeline, cite the gene finder you ran. The references are in
the READMEs of [Tiberius](https://github.com/Gaius-Augustus/Tiberius) and
[Vipsania](https://github.com/Gaius-Augustus/Vipsania). Please also cite the
tools that the pipeline runs on your data: miniprot, miniprot-boundary-scorer,
miniprothint, HISAT2, minimap2, StringTie, TransDecoder, DIAMOND, SAMtools,
BEDTools, gffread, and AUGUSTUS (bam2hints).

## Funding

Funded by the Deutsche Forschungsgemeinschaft (DFG, German Research
Foundation), project number 552910312: "AI-GUSTUS: Eine cloud-native Pipeline
für genaue Genom-Annotation" (AI-GUSTUS: a cloud-native pipeline for accurate
genome annotation).
