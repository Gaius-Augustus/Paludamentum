# Paludamentum

The paludamentum was the cloak worn by Roman emperors and commanders. This
repository is the cloak around the emperors of the
[Gaius-Augustus](https://github.com/Gaius-Augustus) gene finder family: a
Nextflow pipeline that prepares extrinsic evidence (proteins, RNA-Seq,
Iso-Seq), derives high-confidence genes from it, and integrates them with the
ab initio predictions of a deep learning gene finder.

Supported gene finders:

| Gene finder | Status |
| --- | --- |
| [Tiberius](https://github.com/Gaius-Augustus/Tiberius) | supported, launched by `tiberius.py` |
| [Vipsania](https://github.com/Gaius-Augustus/Vipsania) | supported in the pipeline, see [docs/vipsania.md](docs/vipsania.md); launch through `vipsania annotate` is planned |

> **Status (v0.2.0).** The pipeline runs with Tiberius and with Vipsania as
> gene finder. Tiberius uses this repository as a submodule. Launching the
> pipeline through the Vipsania command line is the remaining step, see the
> [Roadmap](#roadmap). Until then, run Vipsania through `python -m paludamentum`.

## Table of contents

- [What the pipeline does](#what-the-pipeline-does)
- [How you get Paludamentum](#how-you-get-paludamentum)
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

## How you get Paludamentum

You normally do not install Paludamentum yourself. It is a **git submodule**
of the gene finder repositories and is checked out at `<repo>/paludamentum/`:

```bash
git clone --recursive https://github.com/Gaius-Augustus/Tiberius
# existing clone:
git submodule update --init --recursive
```

The gene finder's own command line launches the pipeline. Nothing about the
Tiberius command line changes through the migration.

Standalone use is possible for development:

```bash
git clone https://github.com/Gaius-Augustus/Paludamentum
cd Paludamentum
python -m paludamentum --params_yaml params.yaml --nf_config conf/local.config
```

## Requirements

On the machine that launches the pipeline:

- Nextflow
- Java 11 or newer
- Singularity or Apptainer (all tools run in containers, see [Containers](#containers))
- Python 3 with `pyyaml`
- an NVIDIA GPU on the nodes that run the gene finder step

You do not need to install HISAT2, miniprot, StringTie and the other tools
when you use the containers. To run without containers, all tools must be in
your `PATH`, or be configured under `tools:` in the params file. The launcher
option `--check_tools` verifies this.

## Quick start

### Tiberius

```bash
# evidence pipeline, inputs from a params file
python tiberius.py --params_yaml params.yaml --nf_config conf/slurm_generic.config

# ab initio only, parallelized over GPUs by Nextflow
python tiberius.py --nf_config conf/local.config --genome genome.fa --model_cfg eudicotyledons
```

Evidence can also be given on the command line (`--proteins`,
`--odb12Partitions`, `--rnaseq_paired`, `--rnaseq_single`, `--isoseq`,
`--rnaseq_sra_paired`, `--rnaseq_sra_single`, `--isoseq_sra`). Useful launcher
options: `--dry_run` validates inputs and executables without starting
Nextflow, `--resume` continues a previous run, `--work_dir` sets the Nextflow
work directory.

### Vipsania

```bash
python -m paludamentum --genefinder vipsania --params_yaml params.yaml --nf_config slurm_generic
```

with `vipsania: {run: true, model: Fungi}` in the params file, see
[docs/vipsania.md](docs/vipsania.md). Planned:

```bash
vipsania annotate Fungi genome.fa --params_yaml params.yaml --nf_config slurm_generic
```

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

All files are written to `outdir`. `<tool>` is `tiberius` or `vipsania`.

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

The gene finder is selected by which block has `run: true`, or explicitly with
`genefinder: tiberius|vipsania`. The launcher of each gene finder sets this
for you.

### Tiberius block

```yaml
tiberius:
  run: true
  model_cfg: eudicotyledons   # name from Tiberius/model_cfg or a path
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

Vipsania finetuning is **off by default**. With `finetune: true` Vipsania
first trains on the target genome and then annotates it. Vipsania 1.0.0 cannot
finetune without annotating, so in this case the genome is processed in a
single GPU task instead of chunks. All Vipsania parameters, the model download
and offline use are described in [docs/vipsania.md](docs/vipsania.md).

## Containers

Paludamentum does not ship its own image.

| Processes | Image | Built from |
| --- | --- | --- |
| all evidence tools and Tiberius | `docker://larsgabriel23/tiberius:<version>` | `Dockerfile` in the Tiberius repository |
| Vipsania | `docker://gaiusaugustus/vipsania:<version>` | `Dockerfile` in the Vipsania repository |

The images are pinned in [conf/base.config](conf/base.config) through the
process labels `container` and `vipsania`. The pipeline scripts in `bin/` are
not part of an image. Nextflow adds `bin/` to the `PATH` of every task and
mounts it into the container.

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
paludamentum/          Python launcher used by tiberius.py and vipsania
docs/                  parameters.md, hpc.md, vipsania.md
tests/                 launcher tests and Nextflow stub runs
```

## For maintainers

- **Submodule pinning.** Tiberius and Vipsania pin Paludamentum to a tag.
  After a Paludamentum release, bump the submodule pointer in both
  repositories.
- **Image tag coupling.** `conf/base.config` pins the Tiberius image tag,
  while `tiberius.py --singularity` derives its tag from the installed
  Tiberius version. A new Tiberius image therefore needs a Paludamentum commit
  and a submodule bump.
- **Stable interfaces.** Tiberius users rely on the published file names, on
  the `tiberius.*` parameter block, and on `conf/<name>.config`. Do not change
  them without a deprecation path.
- **Adding a gene finder.** Add a module with the run process, a label with
  its container in `conf/base.config`, a parameter block, a branch in
  `subworkflows/genefinder.nf`, and an entry in the launcher's gene finder
  table. The process takes a genome FASTA and emits GTF or GFF3.
  `bin/merge_annotations.py` renumbers gene IDs during merging.
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
`stub:` block for this purpose; keep it in sync when you change outputs. Real smoke tests use
`Tiberius/test_data/Panthera_pardus` and
`Vipsania/docs/example/aspergillus_fumigatus_chr7.fa`.

## Roadmap

The full plan is in [MIGRATION_PLAN.md](MIGRATION_PLAN.md).

- [x] v0.1.0: copy of the pipeline from Tiberius, launcher, stub blocks, tests, CI
- [x] Tiberius uses the submodule; the originals are removed from Tiberius
- [x] v0.2.0: gene finder abstraction and Vipsania processes
- [ ] Vipsania launches the pipeline with `vipsania annotate --params_yaml/--nf_config`
- [ ] Follow-up: `vipsania annotate --finetune_only`, so finetuning can be combined with chunked annotation

## Known issues

Carried over unchanged from Tiberius and tracked for a later fix:

- The automatic mode inference has operator precedence slips. Set `mode`
  explicitly if the inferred mode is not what you expect.
- `rnaseq_bam` alone is not reliably detected as RNA-Seq input by the mode
  inference.
- `restart` and `prothint_conflict_filter` are declared but unused.

## License and citation

MIT, see [LICENSE](LICENSE).

If you use the pipeline, cite the gene finder you ran. The references are in
the READMEs of [Tiberius](https://github.com/Gaius-Augustus/Tiberius) and
[Vipsania](https://github.com/Gaius-Augustus/Vipsania). Please also cite the
tools that the pipeline runs on your data: miniprot, miniprot-boundary-scorer,
miniprothint, HISAT2, minimap2, StringTie, TransDecoder, DIAMOND, SAMtools,
BEDTools, gffread, and AUGUSTUS (bam2hints).
