# Development

## Submodules and status

The gene finders are git submodules of this repository. Paludamentum runs
them; they do not depend on Paludamentum.

| Submodule | Role | Status |
| --- | --- | --- |
| [Tiberius](https://github.com/Gaius-Augustus/Tiberius) (`tiberius/`) | gene finder | supported |
| [Vipsania](https://github.com/Gaius-Augustus/Vipsania) (`vipsania/`) | gene finder | supported, see [vipsania.md](vipsania.md) |
| [Drusilla](https://github.com/Gaius-Augustus/Drusilla) (`drusilla/`) | ORF annotator for assembled transcripts | high-confidence genes for vertebrate models, see [Drusilla flow](drusilla_flow.md) |

> **Status (v0.5.0).** Paludamentum is the entry point: `paludamentum`
> launches the pipeline with Tiberius or Vipsania. Earlier versions were a
> submodule of Tiberius and were launched by `tiberius.py`; that direction is
> reversed since v0.3.0. v0.4.0 adds the [Drusilla flow](drusilla_flow.md) for
> vertebrate models, v0.5.0 the [post-processing and quality
> control](postprocessing.md) of the final annotation.

## Repository layout

```text
main.nf                entry workflow
nextflow.config        manifest, includes conf/base.config
lib_nf/                shared Groovy functions
modules/               Nextflow processes
subworkflows/          inputs, protein, RNA-Seq, Iso-Seq, HC genes, gene finder, post-processing
bin/                   scripts called by processes (bin/gff3_lib.py: GFF3 model of the post-processing scripts)
conf/                  base config, site configs, parameters.yaml, blosum62.csv
paludamentum/          Python launcher (paludamentum, python -m paludamentum)
tiberius/              submodule: Tiberius (gene finder)
vipsania/              submodule: Vipsania (gene finder)
drusilla/              submodule: Drusilla (ORF annotator for transcripts)
docs/                  parameters, Vipsania, Drusilla flow, containers, HPC, comparisons, development
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

The test environment is managed with [pixi](https://pixi.sh) (`pixi.toml`,
no root access needed). It has the Python packages of the tests (the `test`
extra of `pyproject.toml`), samtools, GffRead, GenomeTools and Nextflow:

```bash
pixi run -e test test-scripts   # launcher and scripts, no Nextflow needed
pixi run -e test test           # everything, with the Nextflow stub runs
pixi run -e test nextflow lint main.nf modules subworkflows lib_nf
```

Without pixi, `pip install -e .[test]` installs the Python packages; samtools,
GffRead and GenomeTools (`gt`, or `GT_BIN`) then have to be on the `PATH`, and
the tests that need them are skipped otherwise (`pytest -rs` lists the
reasons). `pytest tests -m "not nextflow"` runs the tests without Nextflow,
`pytest tests -m nextflow` the stub runs (`NEXTFLOW_BIN=/path/to/nextflow`
for another Nextflow).

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
      LightGBM filter of the *ab initio* predictions)
- [x] v0.5.0: post-processing and quality control (docs/postprocessing.md):
  - [x] subworkflow `POSTPROCESS` and its parameter blocks
  - [x] sanity filter, GFF3 contract and validator, GTF, CDS FASTA, longest isoform
  - [x] completeness with compleasm and BUSCO
  - [x] UTRs from the StringTie assemblies
  - [x] hint support, gene set statistics, software versions, `report.html`
  - [x] ncRNA: pybarrnap, tRNAscan-SE, Infernal/Rfam, FEELnc
  - [x] OMArk and gffcompare
  - [x] GO terms with FANTASIA-Lite (GPU)
- [ ] BUSCO rescue of *ab initio* genes dropped by the Drusilla flow (decide with the completeness numbers)
- [ ] `vipsania annotate --finetune_only`, so finetuning can be combined with chunked annotation
