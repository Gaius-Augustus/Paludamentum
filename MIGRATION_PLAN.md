# Migration plan: evidence pipeline from Tiberius to Paludamentum, integration of Vipsania

Status: revised 2026-09-25. Steps 1 to 4 of the first plan are done
(Paludamentum v0.2.0). The dependency direction was then reversed: Paludamentum
imports the gene finders, not the other way around. Steps 5 to 8 below
implement the reversal; step 5 is implemented and unreleased.

## Context

The Nextflow evidence integration pipeline lived inside Tiberius
(`tiberius/main.nf`, `tiberius/modules/`, `tiberius/subworkflows/`,
`tiberius/scripts/`, `conf/`, launched by `tiberius/evidence_pipeline_wrapper.py`).
It moved to its own repository, Paludamentum, so that a second gene finder,
Vipsania, can use the same pipeline.

## Design

Paludamentum is the entry point and imports the gene finders. Tiberius and
Vipsania are gene finders only; nothing in them depends on Paludamentum.

1. Paludamentum has **git submodules** `tiberius/`, `vipsania/` and `drusilla/`
   (Drusilla is the ORF annotator for assembled transcripts, imported for a
   later step of the HC gene workflow). Each submodule is pinned to the release
   whose container image `conf/base.config` runs; the launcher warns when the
   versions differ.
2. The pipeline is launched by `paludamentum` (`python -m paludamentum`):
   `paludamentum --nf_config <config> --genome ... --model_cfg ...` (Tiberius),
   `paludamentum --genefinder vipsania --nf_config <config> --genome ... --model ...`,
   or `paludamentum --params_yaml params.yaml --nf_config <config>`. Command
   line values override the params file; the merged params go to
   `<outdir>/params.yaml`. `tiberius.py` no longer runs the pipeline; its
   pipeline options exit with a pointer to Paludamentum.
3. **Tiberius keeps the fat Dockerfile** (`larsgabriel23/tiberius:<ver>`, all
   evidence tools). Paludamentum ships no image and mounts its `bin/` scripts
   into the tasks. Vipsania processes use `docker://gaiusaugustus/vipsania:<ver>`
   through a dedicated process label.
4. Vipsania **finetuning is off by default**, switchable with
   `vipsania.finetune: true` / `--finetune`.
5. What the launcher takes from the submodules: `tiberius/model_cfg/` to
   resolve model configuration names to files (the Nextflow process stages
   the file), and `tiberius/tiberius.py` on the `PATH` of the Nextflow process
   for runs without containers. Vipsania and Drusilla are pinned sources and
   documentation; the pipeline runs their images.

History: the first plan (2026-09-20) made Paludamentum a submodule of
Tiberius and of Vipsania, launched by `tiberius.py --nf_config` and
`vipsania annotate --nf_config`. That was reversed on 2026-09-25 before
Vipsania was wired, because the pipeline is the integrator and must own the
version pins of the gene finders, not the reverse.

## Target layout of Paludamentum

```text
README.md  LICENSE  .gitignore  .gitmodules  pyproject.toml (hatchling, packages=["paludamentum"], script paludamentum)
main.nf
nextflow.config            # manifest + includeConfig 'conf/base.config'
lib_nf/functions.nf        # truthy(), resolveGenefinder(), inferMode()
modules/                   # evidence processes + genefinder.nf + vipsania.nf
subworkflows/              # inputs, protein, RNA-Seq, Iso-Seq, HC genes, genefinder, ab_initio
bin/                       # pipeline scripts, mode 100755
conf/                      # base/local/slurm_generic/greifswald_hpc/user_hpc_template .config, parameters.yaml, blosum62.csv
paludamentum/              # __init__.py, launcher.py, cli.py, __main__.py
tiberius/  vipsania/  drusilla/   # git submodules
docs/                      # parameters.md, hpc.md, vipsania.md
figures/evi_wflow.png
tests/                     # test_launcher.py, test_stub_run.py, stub.config, data/
.github/workflows/test.yml # checkout with submodules
```

## Rollout

| Step | Repository | Branch | What happens | Ends with |
| --- | --- | --- | --- | --- |
| 1 | Paludamentum | `main` | pipeline copied in, launcher, stubs, tests, CI | done, tag `v0.1.0` |
| 2 | Tiberius | `paludamentum` | submodule added at `v0.1.0`, shim, config shims | done (superseded by step 6) |
| 3 | Tiberius | `paludamentum` | original pipeline files deleted | done |
| 4 | Paludamentum | `main` | gene finder abstraction, Vipsania processes | done, tag `v0.2.0` |
| 5 | Paludamentum | `main` | submodules `tiberius/`, `vipsania/`, `drusilla/`; `paludamentum` command line; model_cfg resolution; docs | commit, push, tag `v0.3.0` |
| 6 | Tiberius | `paludamentum` | submodule, shims and pipeline options removed; docs point to Paludamentum | commit, push, merge to `main`, release |
| 7 | Vipsania | `paludamentum` | docs point to Paludamentum (`docs/pipeline.md`, README section) | commit, push, merge to `main` |
| 8 | Paludamentum | `main` | submodule pointers bumped to the merged Tiberius and Vipsania commits | commit |

Steps 5, 6 and 7 are independent of each other. Step 8 needs the commits of
steps 6 and 7 on GitHub. Until step 8, `tiberius/` points at Tiberius `main`
(v2.0.7 plus the Blackwell Dockerfile change), which still carries the old
in-tree pipeline; that does not affect Paludamentum, which uses only
`model_cfg/` and `tiberius.py` from the checkout.

### Step 5 [Paludamentum]. v0.3.0: import the gene finders

- `git submodule add` Tiberius, Vipsania, Drusilla at `tiberius/`, `vipsania/`, `drusilla/`.
- `paludamentum/launcher.py`: `SUBMODULES`, `submodule_root()`, `submodule_version()`,
  `container_tags()`, `version_mismatches()` (warning at launch), `tiberius_checkout()`,
  `nextflow_env()` (Tiberius checkout appended to `PATH`), `resolve_model_cfg()`
  (name to file via `tiberius/model_cfg/`, superseded configs with a warning),
  `infer_genefinder()`, `cli_overrides()`, `build_params()` (defaults, params file,
  command line, written to `<outdir>/params.yaml`). `--check_tools` accepts
  `tiberius.py` from the checkout.
- `paludamentum/cli.py`: full command line (Nextflow options, inputs, gene finder
  options, Tiberius `--model_cfg --seq_len`, Vipsania `--model --model_dir --finetune
  --finetune_epochs --context`); `--params_yaml` optional. `pyproject.toml`: console
  script `paludamentum`, sdist excludes the submodules.
- CI checks out with `submodules: true`. Tests: submodule table, version/tag
  consistency, model_cfg resolution, `build_params`, CLI dry run, CLI stub run.
- README, docs, MIGRATION_PLAN rewritten for the reversed direction.
- Commit, push, tag `v0.3.0`.

### Step 6 [Tiberius]. Remove the reverse dependency

- `git submodule deinit paludamentum`, `git rm paludamentum .gitmodules`, `git rm -r conf`
  (the shims), `git rm tiberius/evidence_pipeline_wrapper.py tests/unit_tests/pipeline/`.
- `tiberius.py`: the `nextflow` mode, `DEFAULT_PARAMS`, `collect_cli_params`,
  `ensure_params_yaml`, `merge_dicts` removed. `reject_removed_pipeline_flags()` exits
  with a pointer to Paludamentum when one of the removed options is used.
  `--params_yaml` stays: it fills `--genome` and `--model_cfg` from a params file.
- `tiberius/tiberius_args.py`: groups "Nextflow Pipeline" and "Nextflow Params" removed.
- `Dockerfile`: plain clone, no `paludamentum/bin` on `PATH`. The image keeps all
  evidence tools. `pyproject.toml`, `pytest.ini`, CI: no submodule references.
- README: the Nextflow sections replaced by "Multi-GPU runs and extrinsic evidence
  (Paludamentum)". New test `tests/unit_tests/launcher/test_removed_pipeline_flags.py`.
- Commit, push, merge to `main`.

### Step 7 [Vipsania]. Documentation only

- README section "Evidence integration (Paludamentum)", `docs/pipeline.md`, row in
  `docs/README.md`. No code, no submodule, no launcher flags. Commit, push, merge.

### Step 8 [Paludamentum]. Bump the submodule pointers

- `git -C tiberius fetch && git -C tiberius checkout <merged commit or tag>`; same for
  `vipsania`; `git add tiberius vipsania`. Keep the image tags in `conf/base.config` in
  step with the versions. Commit.

### Follow-ups, not part of this plan

- Drusilla as alternative to TransDecoder in `subworkflows/hc_genes.nf` (process with
  a `drusilla` label and image `larsgabriel23/drusilla`).
- `vipsania annotate --finetune_only`, to allow finetuning once plus chunked annotation.
- Fix the mode-inference quirks.

## Verification

| Check | Paludamentum | Tiberius | Vipsania |
| --- | --- | --- | --- |
| A. pytest | steps 5, 8 | step 6 | |
| B. launcher dry run | steps 5, 8 | | |
| C. Nextflow stub runs | steps 5, 8 | | |
| D. GPU smoke test | step 8 | step 6 (direct mode unchanged) | |

- **A. pytest (no Nextflow).** Paludamentum `tests/test_launcher.py`: submodule
  table, version/tag consistency, `resolve_model_cfg`, `build_params`, command
  assembly, CLI dry run. Tiberius `tests/unit_tests/launcher/`: the removed options
  exit with the pointer; no pipeline files left.
- **B. Dry run.** `paludamentum --nf_config slurm_generic --genome <genome.fa>
  --model_cfg mammalia_softmasking_v2 --dry_run --skip_singularity_check`; assert that
  `tiberius_results/params.yaml` has `model_cfg` under `tiberius/model_cfg/`.
- **C. Nextflow stub runs.** `pytest tests/test_stub_run.py` with `NEXTFLOW_BIN`; per
  mode and per gene finder, plus one run through the `paludamentum` command line.
- **D. GPU smoke tests.** Run B and C before any cluster submission. Tiberius: ab
  initio on `tiberius/test_data/Panthera_pardus` with a small `min_split_size` to
  force 2+ chunks; check `.command.run` for the `bin` bind. Vipsania:
  `vipsania/docs/example/aspergillus_fumigatus_chr7.fa` (2 Mb); compare the gene
  count with `docs/example/vipsania_fh1kg88z.gff`; repeat with `--finetune
  --finetune_epochs 1` and with a small protein FASTA.

## Risks

- `base.config` pins `tiberius:2.0.8`; the submodule pins the matching Tiberius commit. Each
  Tiberius release needs one Paludamentum commit that bumps both.
- Existing Tiberius users who call `tiberius.py --nf_config` get an error with the
  new command; their params files work unchanged with `paludamentum --params_yaml`.
- The `bin/` bind under Singularity `--contain` relies on Nextflow's automatic mount
  of the project `bin/`; the image no longer holds a copy of the scripts.
- `DOWNLOAD_VIPSANIA_MODEL` needs internet on the executing node; `vipsania.model_dir`
  is the offline path.
- Vipsania's image does not support Blackwell GPUs; CPU-only runs need `vipsania.batch_size`.
