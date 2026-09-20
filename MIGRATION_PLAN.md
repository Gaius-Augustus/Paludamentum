# Migration plan: evidence pipeline from Tiberius to Paludamentum, integration of Vipsania

Status: approved 2026-09-20, implementation not started.

## Context

The Nextflow evidence integration pipeline lives inside Tiberius today
(`tiberius/main.nf`, `tiberius/modules/`, `tiberius/subworkflows/`,
`tiberius/scripts/`, `conf/`, launched by `tiberius/evidence_pipeline_wrapper.py`).
It moves to its own repository, Paludamentum, so that a second gene finder,
Vipsania, can use the same pipeline. Tiberius usage must not change: same
commands, same flags, same output file names. Vipsania has never run the
pipeline.

Design decisions:

1. Tiberius and Vipsania obtain Paludamentum as a **git submodule** at `<repo>/paludamentum/`.
2. **Tiberius keeps the fat Dockerfile** (`larsgabriel23/tiberius:<ver>`, all evidence tools). Paludamentum ships no image. Vipsania is not in that image, so Vipsania processes use `docker://gaiusaugustus/vipsania:1.0.0` through a dedicated process label.
3. The Vipsania launcher **mirrors Tiberius**: `vipsania annotate --params_yaml/--nf_config ...`.
4. Vipsania **finetuning is off by default**, switchable with `vipsania.finetune: true` / `--finetune`.

Submodule pointers need a pushed Paludamentum commit. Paludamentum is therefore
pushed and tagged at the end of step 1 (before Tiberius is wired) and at the
end of step 4 (before Tiberius is bumped and Vipsania is wired).

Facts that shape the design (checked against Tiberius@66c228b and Vipsania@413f756):

- The wrapper returns only an exit code; Tiberius reads nothing back. The only path coupling is `pipeline_paths()` (`evidence_pipeline_wrapper.py:70-78`), `tiberius.py:50` (blosum path), `tiberius.py:482` (default base.config), and `conf/base.config:30` (`${projectDir}/../conf/blosum62.csv`).
- All 13 files in `tiberius/scripts/` are pipeline-only. They reach PATH only through `Dockerfile:64,85`. In Paludamentum they go to `bin/`, which Nextflow adds to PATH and bind-mounts automatically (`singularity.autoMounts = true` is already set).
- The user configs (`local`, `slurm_generic`, `greifswald_hpc`, `user_hpc_template`) each `includeConfig 'base.config'`, so `conf/` moves as a unit.
- The gene finder seam is `RUN_TIBERIUS` (`modules/tiberius.nf:1-25`), called from `protein_evidence.nf:31-44` and `tiberius_only.nf`. `merge_annotations.py` renumbers gene IDs, so Vipsania's per-sequence ID restart is harmless, and it accepts Vipsania GTF.
- Nextflow config selectors override in-process directives. Vipsania processes must therefore not carry the `container` label.
- Vipsania 1.0.0 has no finetune-only mode. With finetune on, the pipeline runs one whole-genome `vipsania annotate --finetune` task instead of chunking.
- `vipsania annotate` has required positionals and already uses `-p`. The Nextflow flags are long-only, positionals become optional, and extra Nextflow args are a quoted `--nextflow_args` string.
- No pipeline tests and no Nextflow CI exist anywhere today.

## Target layout of Paludamentum

```text
README.md  LICENSE  .gitignore  pyproject.toml (hatchling, packages=["paludamentum"], dep pyyaml)
main.nf
nextflow.config            # manifest + includeConfig 'conf/base.config'
lib_nf/functions.nf        # truthy(), resolveGenefinder(), inferMode()
modules/                   # 9 unchanged + genefinder.nf (ex tiberius.nf) + vipsania.nf (new)
subworkflows/              # 5 unchanged + genefinder.nf (new) + ab_initio.nf (ex tiberius_only.nf)
bin/                       # ex tiberius/scripts, all mode 100755
conf/                      # base/local/slurm_generic/greifswald_hpc/user_hpc_template .config, parameters.yaml, blosum62.csv
paludamentum/              # __init__.py, launcher.py (ex evidence_pipeline_wrapper.py), cli.py, __main__.py
docs/                      # parameters.md + hpc.md (ex conf/README.md), vipsania.md
figures/evi_wflow.png
tests/                     # test_launcher.py, test_stub_run.py, params_*.yaml, data/tiny.fa
.github/workflows/test.yml
```

## Rollout

Tiberius keeps working after every step. Every step touches exactly one
repository. All file paths inside a step are relative to that repository.

| Step | Repository | Branch | What happens | Ends with |
| --- | --- | --- | --- | --- |
| 1 | Paludamentum | `main` | pipeline copied in, launcher, stubs, tests, CI | commit, push, tag `v0.1.0` |
| 2 | Tiberius | `paludamentum` | submodule added at `v0.1.0`, shim, config shims, Dockerfile, README | commit |
| 3 | Tiberius | `paludamentum` | original pipeline files deleted | commit, merge to `main` |
| 4 | Paludamentum | `main` | gene finder abstraction, Vipsania processes | commit, push, tag `v0.2.0` |
| 5 | Tiberius | `main` | submodule pointer bumped to `v0.2.0` | commit |
| 6 | Vipsania | `paludamentum` | submodule added at `v0.2.0`, launcher flags, docs | commit |

Steps 2 and 6 cannot start before the tag of the preceding Paludamentum step
is pushed, because the submodule pointer needs that commit on GitHub. Vipsania
is not touched before step 6.

### Step 1 [Paludamentum]. v0.1.0: pure copy, Tiberius only

- Copy from Tiberius@66c228b per the layout above and name the source commit in the commit message. Set the exec bit on every `bin/*`.
- `conf/base.config:30` becomes `scoring_matrix = "${projectDir}/conf/blosum62.csv"`. Declare the currently undeclared `tiberius.min_split_size`, `tiberius.max_files`, `rnaseq_bam`.
- `paludamentum/launcher.py`: generalize the wrapper.
  - `pipeline_paths(root)` resolves `<root>/main.nf` and `<root>/conf/base.config`; root defaults to `Path(__file__).parent.parent`.
  - New `resolve_nf_config(value, root)`: existing path, else `root/conf/<basename>`, else `root/conf/<value>.config`, else an error listing the search.
  - `run_nextflow_pipeline(args, genefinder="tiberius", pipeline_root=None)`; `GENEFINDER_CLI = {"tiberius": "tiberius.py", "vipsania": "vipsania"}` replaces the hard-coded check at old line 296.
  - Keep the command shape `nextflow run main.nf -params-file P -c base.config -c user.config ...`.
- `paludamentum/cli.py`: minimal `python -m paludamentum --params_yaml ... --nf_config ... [--genefinder]`.
- Add `stub:` blocks (touch the declared outputs) to every process so `-stub-run` works without tools or GPU.
- Add docs, tests, CI (pytest + `nf-core/setup-nextflow` stub runs).
- Commit, push, tag `v0.1.0`.

### Step 2 [Tiberius]. Wire Tiberius to the submodule (old files still present, unused)

- `.gitmodules`: remove the stale `learnMSA` stanza; `git submodule add https://github.com/Gaius-Augustus/Paludamentum paludamentum`, pinned to `v0.1.0`.
- `tiberius/evidence_pipeline_wrapper.py`: replace with a ~25-line shim. It inserts `<repo>/paludamentum` at `sys.path[0]` (avoids the namespace-package shadowing by the outer directory), imports `paludamentum.launcher` lazily, and exits with a `git submodule update --init --recursive` hint if `main.nf` is missing. It exposes `pipeline_paths`, `resolve_nf_config`, and `run_nextflow_pipeline(args)`, which calls the launcher with `genefinder="tiberius"`.
- `tiberius.py`: add `PALUDAMENTUM_ROOT = SCRIPT_ROOT / "paludamentum"`; the blosum path (line 50) and the default base.config (line 482) point under it; resolve `args.nf_config` through `resolve_nf_config` before `ensure_params_yaml`. `collect_cli_params`, `ensure_params_yaml`, `resolve_model_cfg`, and `tiberius_args.py:89-139` stay untouched except the help text at line 126.
- `conf/`: keep the directory with one-line shims, `includeConfig '../paludamentum/conf/<name>.config'`, for the 5 configs. This keeps `--nf_config conf/slurm_generic.config` and users' own configs in `conf/` working unchanged. Drop `blosum62.csv` and `parameters.yaml`; reduce `conf/README.md` to a pointer.
- `Dockerfile`: line 60 `git clone --recursive`; line 64 add `chmod +x paludamentum/bin/*`; line 85 PATH `/opt/Tiberius/paludamentum/bin` (fallback if the host bin mount is ever absent). No image rebuild is needed for the move itself.
- `pyproject.toml`: `[tool.hatch.build.targets.wheel] packages = ["tiberius"]`; sdist excludes `paludamentum`, `test_data`.
- `README.md`: `git clone --recursive`, and `git submodule update --init` for existing clones; re-point the links at lines 218, 231, 245 into `paludamentum/`; fix the wrong default outdir at line 249 (`tiberius_results`).
- `.github/workflows/test.yml`: checkout with `submodules: recursive`. Add `tests/unit_tests/test_pipeline_shim.py`.
- Commit.

### Step 3 [Tiberius]. Delete the originals

Remove `tiberius/main.nf`, `tiberius/modules/`, `tiberius/subworkflows/`,
`tiberius/scripts/`, `figures/evi_wflow.png`. Re-run verification A to C.
Commit and merge `paludamentum` into `main`. The next Tiberius image build
picks up the Dockerfile change.

### Step 4 [Paludamentum]. v0.2.0: gene finder abstraction plus Vipsania processes

- `lib_nf/functions.nf`: `truthy()` replaces three copies; `resolveGenefinder(params)` returns `params.genefinder`, else `vipsania` if `vipsania.run`, else `tiberius`; it errors if both run flags are set. Move `inferMode` verbatim. Do not fix the known quirks (precedence at old `main.nf:17-18`, `rnaseq_bam?.size`) in this step; they are listed in the README as known issues.
- `modules/genefinder.nf` (ex `tiberius.nf`): `RUN_TIBERIUS` unchanged. `SPLIT_GENOME` takes `min_size`, `max_files` as inputs. Merge and protein processes take `val prefix` and emit `${prefix}_ab_initio.gff3`, `${prefix}_evidence.gff3`, `${prefix}_evidence_proteins.fa`. With prefix `tiberius` all published names equal today's.
- `modules/vipsania.nf`, label `vipsania` (never `container`), with a shared preamble that exports `LD_LIBRARY_PATH` from `site-packages/nvidia/*/lib` (Singularity bypasses the image ENTRYPOINT), `VIPSANIA_CACHE=$PWD/.vipsania_cache`, `WANDB_MODE=disabled`:
  - `DOWNLOAD_VIPSANIA_MODEL` (`local_only`): `vipsania download <model> -d models`, emits `models/*`.
  - `RUN_VIPSANIA` (`gpu`, `bigmem`, `maxForks` from `vipsania.max_parallel`): `vipsania annotate <model> <genome> --model_dir vip_models -o vipsania.<name>.gtf [-B] [-T] [--finetune ...] [extra_args]`. It also emits the `.log` sidecar and, with finetune, `finetuning_*/` to `intermediate/`.
- `subworkflows/genefinder.nf`: `GENEFINDER(genome, params, publish_top)`. Reuse `<tool>.result` if given. Tiberius: split, `RUN_TIBERIUS`, merge. Vipsania: models from `vipsania.model_dir` or the download process; finetune off means split, per-chunk `RUN_VIPSANIA`, merge; finetune on means one whole-genome `RUN_VIPSANIA --finetune`, then merge. Emits `gff`, `tool`.
- `protein_evidence.nf:26-47` and `ab_initio.nf` call `GENEFINDER`; `main.nf:56` accepts modes `tiberius|abinitio|vipsania`; `main.nf:99-103` uses the prefixed merge processes.
- `conf/base.config`: add `genefinder = null`, a `vipsania { run=false; model; model_dir; result; finetune=false; finetune_epochs; finetune_B; finetune_lr; batch_size; context; max_parallel; min_split_size; max_files; extra_args }` block, and `withLabel: vipsania { container = 'docker://gaiusaugustus/vipsania:1.0.0' }`. All `.nf` access is null-safe (`params.vipsania?.x ?: default`).
- Launcher: `default_params(genefinder)` and `write_params_yaml()` helpers for Vipsania.
- Regression gate: stub runs with Tiberius produce exactly the old published names. Release note: process renames invalidate `-resume` for runs in progress.
- Commit, push, tag `v0.2.0`.

### Step 5 [Tiberius]. Bump the submodule to v0.2.0

- `git -C paludamentum fetch --tags && git -C paludamentum checkout v0.2.0`, then stage the new pointer.
- No other Tiberius file changes. Re-run verification A to C; the published file names must be unchanged.
- Commit.

### Step 6 [Vipsania]. Pipeline integration

- `git submodule add ... paludamentum`, pinned to `v0.2.0`. `pyproject.toml`: add `pyyaml`; sdist excludes `paludamentum`.
- `vipsania/cli/annotate.py`: `model`, `fasta` get `nargs="?"` (lines 281-288). After line 488 add the group "Nextflow pipeline (Paludamentum)": `--params_yaml --nf_config --profile --nextflow_bin --resume --work_dir --check_tools --skip_singularity_check --dry_run --nextflow_args`; and the group "Nextflow params": `--outdir --threads --proteins --odb12Partitions --rnaseq_single --rnaseq_paired --rnaseq_sra_single --rnaseq_sra_paired --isoseq --isoseq_sra --mode --scoring_matrix --vipsania_result --max_parallel`. `run()` (line 491) dispatches to the pipeline if `--nf_config` or `--params_yaml` is given; otherwise it requires the positionals, so direct mode is unchanged.
- New `vipsania/cli/pipeline.py`: `find_paludamentum()` (`$PALUDAMENTUM_HOME`, then `<repo>/paludamentum`, else a clone hint); `build_params(args)` merges launcher defaults, YAML, and CLI into `genome`, `vipsania: {run: true, model, model_dir, finetune, ...}`, `tiberius: {run: false}`; writes `<outdir>/params.yaml` (default outdir `vipsania_results`); calls the launcher with `genefinder="vipsania"`.
- Docs: `docs/pipeline.md` and a README section. The Vipsania Dockerfile is unchanged.
- Commit.

### Follow-ups, not part of this plan

- `vipsania annotate --finetune_only`, to allow finetuning once plus chunked annotation.
- Let the Tiberius launcher override the container tag pinned in `base.config`.
- Fix the mode-inference quirks.

## Verification

Which check runs in which repository, and after which step:

| Check | Paludamentum | Tiberius | Vipsania |
| --- | --- | --- | --- |
| A. pytest | steps 1, 4 | steps 2, 3, 5 | step 6 |
| B. launcher dry run | | steps 2, 3, 5 | step 6 |
| C. Nextflow stub runs | steps 1, 4 | steps 2, 3, 5 (through the submodule) | step 6 (through the submodule) |
| D. GPU smoke test | | steps 3, 5 | step 6 |

- **A. pytest (no Nextflow).**
  - Paludamentum, `tests/test_launcher.py`: path resolution, `resolve_nf_config`, input validation globs, command assembly via monkeypatched `subprocess.run`, gene finder CLI check.
  - Tiberius, `tests/unit_tests/test_pipeline_shim.py`: paths exist, a missing submodule gives a clean exit.
  - Vipsania: parser test that direct mode still requires the positionals, and a `build_params` snapshot.
- **B. Dry run.** Needs Nextflow and Java on the machine.
  - Tiberius: `python tiberius.py --nf_config conf/slurm_generic.config --genome <Panthera genome.fa> --model_cfg mammalia_softmasking_v2 --dry_run --skip_singularity_check`; assert that `tiberius_results/params.yaml` has `scoring_matrix` under `paludamentum/conf/`.
  - Vipsania: `vipsania annotate Fungi docs/example/aspergillus_fumigatus_chr7.fa --nf_config local --dry_run`.
- **C. Nextflow stub runs.** In Paludamentum, or in `paludamentum/` of a gene finder checkout: `nextflow run main.nf -stub-run -params-file tests/params_<mode>.yaml -c conf/local.config` with Singularity disabled, per mode (abinitio, proteins, rnaseq, mixed) and per gene finder; assert the published file names.
- **D. GPU smoke tests.** Run cheap static checks and B/C before any cluster submission.
  - Tiberius: ab initio on `test_data/Panthera_pardus` with a small `min_split_size` to force 2+ chunks. Check `.command.run` for the `bin` bind.
  - Vipsania: `docs/example/aspergillus_fumigatus_chr7.fa` (2 Mb); compare the gene count with `docs/example/vipsania_fh1kg88z.gff`; repeat with `--finetune --finetune_epochs 1` and with a small protein FASTA. Check Vipsania's log for a detected GPU.

## Risks

- `base.config` pins `tiberius:2.0.7` inside Paludamentum while Tiberius derives its own tag from the package version. Each Tiberius release therefore needs a Paludamentum commit and a submodule bump.
- Existing clones that `git pull` without `git submodule update --init` hit the shim's error message. This is documented in the README.
- The `bin/` bind under Singularity `--contain` is expected to work through `autoMounts` but is not yet verified; the image PATH fallback covers it.
- `DOWNLOAD_VIPSANIA_MODEL` needs internet on the executing node; `vipsania.model_dir` is the offline path.
- Vipsania's image does not support Blackwell GPUs; CPU-only runs need `vipsania.batch_size`.
