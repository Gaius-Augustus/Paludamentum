# Paludamentum pre-release audit

Audited 2026-09-29 (`main` at `18110c1`, `drusilla` at `ea68c20`). Updated
2026-09-30 after the overnight fixes and again 2026-10-01: fixed items are
removed; what remains is open. The fixed items are summarised at the end.

State 2026-10-01: the `main` fixes are committed and pushed (`29a00c4`). The
`drusilla` worktree `../Paludamentum_drusilla` is at `3513075` (pushed) with
the script fixes of 2026-09-30 and the model loader of 2026-10-01 (D1, D5)
still **uncommitted** (six `bin/` scripts, `lib_nf/functions.nf`,
`docs/hint_rescue_comparison.md`, `tests/test_stub_run.py`, new
`tests/test_drusilla_scripts.py`).

Scope: launcher and packaging, Nextflow code and configs, `bin/` scripts, the
Drusilla branch diff, CI, docs and licences. The submodules were only checked
where Paludamentum depends on them.

Tags:

- **[repro]** marks a failure reproduced on a toy input or in a stub run.
- All other items were checked by reading the code.

What works (checked 2026-09-30 on the working trees):

- Python tests pass: 91 on `main` (launcher, modules, scripts), 7 script
  tests on `drusilla`.
- Nextflow stub runs pass (32 on `main`, plus 10 new ones for the guards
  added tonight) with Nextflow 25.04.6.
- `nextflow lint main.nf modules subworkflows lib_nf` is clean on both trees.
- `conf/local.config` runs a stub of the ab initio mode on an 8-CPU laptop.
- The Tiberius and Vipsania images are pinned (`tiberius:2.0.8`,
  `vipsania:1.0.0`) and can be pulled publicly.
- No secrets are in the repo.
- No shell injection: `subprocess` is always called with argument lists and
  YAML is read with `safe_load`.
- CI passed on `29a00c4` (`main`) in all four jobs: Python 3.9 and 3.12,
  stub runs on Nextflow 25.04.0 and latest-stable. CI has not run on
  `drusilla`: its own workflow file triggers only on `main` and PRs, so it
  runs with the merge PR.

---

## Decisions needed before the release

1. **LightGBM model for the Drusilla flow.** RESOLVED 2026-10-01: the
   pickle was converted to a native LightGBM text model plus a JSON sidecar
   (48 training species, script names, labels, hyperparameters; no paths) and
   published at
   `https://bioinf.uni-greifswald.de/bioinf/drusilla/models/drusilla_lgb_3class_v1.tar.gz`
   (sha256 `d5bf3a97…2f1da`, checked after the upload). The drusilla worktree
   downloads it by default and loads it without the pickle (D1, D5), in the
   image `gaiusaugustus/drusilla:0.1.0` with lightgbm 4.7.0 (L7, D2).
2. **HC gene logic.** RESOLVED 2026-09-29: the intrinsic stage never changed
   the output and is removed (`5eada8e`). Tonight the dead intrinsic code was
   deleted from `bin/hc_module.py` as well. `choose_one_isoform` keeps one ORF
   per StringTie transcript, not per gene; changing that would drop about
   4,100 HC transcripts per genome and needs a benchmark.
3. **Internal files:**
   - `MIGRATION_PLAN.md` is still tracked and the README Roadmap links to
     it. It mentions "storm" and internal steps. Remove it, or rewrite it
     for the public.
   - `conf/greifswald_hpc.config`: keep as a labelled example (docs/hpc.md
     now calls it that), or drop.
4. **Benchmark again** (T. rubripes and Bos taurus) before the release: the
   merge (B5, B6) and the stop/start fix (D4) now change the output. Expected
   effect: fewer wrong isoforms and no CDS with an in-frame stop; gene-level
   F1 should move little, but it must be measured.

---

## 1. Release logistics

| # | Item | Where |
|---|------|-------|
| L1 | The GitHub repo is still **private**. Make it public last. | GitHub |
| L2 | `MIGRATION_PLAN.md` (see decision 3). Everything else with internal names (test docstring, `docs/orf_finder_comparison.md`, `docs/hint_rescue_comparison.md`) was rewritten. | `MIGRATION_PLAN.md:81`, README Roadmap |
| L6 | Submodules are not pinned to releases, although README "Submodule pinning" says they are. Checked 2026-09-30: Tiberius' newest tag is still `v2.0.7` (submodule at `v2.0.7-11-g9734138`, image 2.0.8); Vipsania and Drusilla have no tags. Ask for tags upstream, then point the submodules at them. | `.gitmodules`, README "For maintainers" |
| L8 | Version bump to 0.4.0 with the Drusilla merge. The version now has one source, `paludamentum/__init__.py` (pyproject reads it); `nextflow.config`, the README status line and `CITATION.cff` must follow, which `test_version_is_the_same_everywhere` enforces. | `paludamentum/__init__.py`, `nextflow.config:11`, README:19, `CITATION.cff` |
| L10 | `CITATION.cff` and `pyproject.toml` `authors` name Katharina J. Hoff only. **Confirm the author list** (Lars Gabriel wrote the Tiberius pipeline the scripts come from). The drusilla README does not cite Drusilla, LightGBM or the hint rescue. | `CITATION.cff`, `pyproject.toml`, README (drusilla) |

---

## 2. Drusilla flow (branch `drusilla`)

The `hint_rescue` image wiring in `README.md`, `conf/base.config`,
`modules/drusilla.nf` and `subworkflows/drusilla.nf` is committed (`3513075`)
and the image `gaiusaugustus/paludamentum-hint-rescue:0.1.0` is on Docker
Hub, so these files are no longer held. The items below in them are still
open.

### BLOCKER

| # | Item | Where |
|---|------|-------|
| D3 | **The README claim about the GCB 2026 poster is wrong.** The poster used `epoch_74`, no stop/start fix, no rescue, a `cat` merge, and `cds_length_nt` filled with zeros. The released `vertebrates` model is run009. Say "derived from" and give the benchmark numbers of this flow (T. rubripes 80.26, Bos 79.46; to be re-measured, decision 4). | README:236-238 |

### SHOULD

| # | Item | Where |
|---|------|-------|
| D4 | Done: a corrected CDS is now written only if it has no in-frame stop codon (`has_internal_stop`, tested). Still open: the splice sites of extension exons are not checked. | `bin/fix_stop_by_miniprot.py` |
| D6 | The ORF-agreement filter of the hint rescue is switched off (benchmarked that way). The script docstrings and `docs/hint_rescue_comparison.md` now say so. Still wrong: the module comments ("each with the hints of its best protein chain", "the predictions that agree with the hints") and the README sentence "hints of the best protein chain". Loci without a chain run ab initio without hints. | `modules/drusilla.nf:202`, `:254`, README:265, `conf/base.config:99` |
| D8 | `DRUSILLA_ANNOTATE` and `HINT_RESCUE_TIBERIUS` have no memory label, so they get the 30 GB default (now growing with the attempt). Measured: 40 GB for Bos, more than 90 GB for the rescue at seq_len 400k. Give them `bigmem`, or their own label. | `modules/drusilla.nf:40`, `:258`, `conf/base.config` |
| D9 | With `shards > 1` and `cache_dir: null`, all shards download and extract into the same `drusilla_cache` at once; the Drusilla registry `rmtree`s before it extracts, so the shards race. Download the model once before the fork. The weights URL has no sha256 in the manifest; which run the weights are (run009) is recorded nowhere. Nodes without internet fail. | `modules/drusilla.nf:56`, `drusilla/model_cfg/vertebrates.yaml` |
| D10 | Remaining gaps in `hcMethod()` after tonight (`run: "on"` is now true; `tiberius.result` without `model_cfg` now warns and uses TransDecoder): <br>• Forcing `run: true` with `tiberius.result` and no model sets the rescue model to the string `"null"`. <br>• A custom model path counts as a vertebrate model by its file name (`/my/vertebrates.yaml`); read `target_species` from the YAML instead, or document it. | `subworkflows/drusilla.nf:43-44`, `lib_nf/functions.nf` drusillaModelEligible |
| D12 | Tests added: `gff_to_cds_gtf`, the in-frame stop check, feature/model column alignment, the empty feature table; on 2026-10-01 the model loader (archive, directory, .txt, both checksums, pickle refused), a missing feature, the columns of `compute_orf_features.py` against the 26 features of the released model, and the filter end to end on a tiny model. These need lightgbm and pandas, which CI does not install (`.[test]` is pytest only), so CI skips them. Still without tests: `compute_orf_features.py`, `filter_stringtie_gtf.py`, `prepare_hint_rescue_loci.py`, `filter_and_merge_rescue_gtf.py`, and stub runs of `mammalia*`, `Vertebrata`, `shards > 1`, `fix_stop: false`, `lgb_keep`, `tiberius.result` with Drusilla. | `tests/test_drusilla_scripts.py`, `tests/test_stub_run.py` |

### NIT

| # | Item |
|---|------|
| D13 | The `drusilla:` config snippet in the README leaves out `batch_size` and `min_coding_length` (README:275-303). |
| D14 | `--proteins-fasta` is never passed, so `best_protein_coverage` is always 0. Check with Lars whether the model was trained with this feature (`modules/drusilla.nf`). |

### Merge notes

Merging `main` into `drusilla` (2026-10-03, state `d3027cf` into `d89a27d`)
conflicts in `main.nf`, `conf/base.config`, `README.md`,
`tests/test_stub_run.py` and the `drusilla` submodule. Resolution: the
include lines of both sides without `MERGE_GENEFINDER_TRAIN_PRIO` (removed on
`main`); `HC_GENES` with three inputs (its `scored_gff` input was removed on
`main`); the `drusilla` params block before the `mode` comment of `main`; both
test blocks; "Known issues" with the Drusilla items only (`main` fixed the
others); the submodule at `e3c5cf5`, the commit of the image. `hcMethod()`
gets the normalised mode. The `drusilla` params block is now in
`conf/parameters.yaml` and `docs/parameters.md`.

The headers of Lars Gabriel's scripts in `bin/` of the drusilla branch say
"Artistic License 1.0, see LICENSE", but the branch still has the MIT
`LICENSE` and the README says MIT: both come right with the merge of `main`.

---

## 3. Launcher, packaging, configs, docs

### NIT

- `params.yaml` is written before validation, so a failed `--dry_run` leaves
  it behind (`cli.py`). Harmless; restructuring `build_params` to validate
  before writing is the fix.
- `--nf_config` still accepts a bare name that matches a file in the launch
  directory (`./local.config` wins over `conf/local.config`); a test locks
  this in. A mistyped *path* now errors instead of falling back.

---

## 4. Nextflow workflow and configs

### SHOULD

| # | Item | Where |
|---|------|-------|
| N10 | `tiberius.model_dir` (added tonight) stages the extracted weights into the task directory, where Tiberius' `download_weights()` finds `<model>_weights` and skips the download. Checked against the Tiberius code and in a stub run only: **run it once on a GPU node without internet** before documenting it as supported. | `modules/genefinder.nf`, `subworkflows/genefinder.nf`, `tiberius/tiberius/main.py:357-395` |
| N14 | The retry rule now retries only exit codes 137, 140, 143 and 247 with memory × attempt, and ends the run on anything else (`finish`). SLURM reports OOM kills as 137 in most setups, but some sites use other codes; check the first real failure on your cluster. | `conf/base.config` |
| N15 | `greifswald_hpc.config` holds site partitions (`snowball,pinky,batch,vision`) and is listed as a shipped config (decision 3). | `conf/greifswald_hpc.config`, `docs/hpc.md` |

### NIT

- N19: the shell blocks still do not quote file variables. The launcher now
  rejects paths with whitespace, and the README says so; a params file that
  bypasses the launcher is not protected.
- `main.nf` has no `-resume`-safe guard against renamed processes; mention
  renames in release notes (already in README "For maintainers").

---

## 5. `bin/` scripts (main)

### Inherited HC logic

Resolved (decisions 2): `training.gff` = `choose_one_isoform(P)`, where P is
the strict DIAMOND set. The intrinsic stage and its helpers were deleted from
`hc_module.py` tonight (they were dead code since `5eada8e`), which removed
the items about `hc_module.py:654` (overlap test), `:440-441` (stop in 5'
UTR) with them. Kept by decision, benchmark needed to change:

- `hc_module.py` `choose_one_isoform` groups by transcript, not by gene, so
  HC training genes overlap (about 4,100 more HC transcripts per genome reach
  the final merge).

### NIT

- `extend_cds_with_stop_codon.py`: a stop codon split by an intron extends
  the CDS across the intron. `check_stop_codons` removes these genes later,
  so the effect is a lost gene, not a wrong one.
- `merge_annotations.py --mode priority` and
  `shorten_incomplete.py --keep-non-incomplete` are unused by the pipeline
  (kept: they work and are cheap).
- Upstream, in the Tiberius Dockerfile: TransDecoder master is cloned without
  a pin, next to v5.7.1.

---

## Fixed 2026-09-30 (for the release notes)

Launcher (`paludamentum/`): the pipeline root is found from the checkout or
`PALUDAMENTUM_ROOT`, with a clear error after a plain `pip install .` (P1);
`--genefinder` overrides the params file and switches the other finder off
(P2), two `run: true` blocks are an error (P4), `--genefinder` with
`run: false` switches the block on and says so (P5); a mistyped `--nf_config`
path errors (P6); `tiberius.py` of the checkout is called through an
executable wrapper (P7); `~`, `$VARS` and relative paths are written absolute
(P8); `rnaseq_bam`, `*.result`, `*.model_dir`, `tiberius.model_cfg` are
validated and paths with whitespace rejected (P9, N19); `--mode` has choices
(P10); Java 17 (P13); `apptainer` satisfies the container check (P14);
`--threads 0` is a value, three `--rnaseq_paired` files and non-string paths
are errors, no more unused arguments or commented code (NITs). Single version
source in `pyproject.toml` (L8), `authors` set.

Nextflow: `nextflowVersion >= 25.04.0` (N1); decompressed files keep their
names (N2); `slurm_generic.config` books a GPU (N3); `CALC_ALIGNMENT_RATE`
and the `CONCAT_*` processes run in the image, and all tasks run with
`pipefail` (N4, N9); a run without any assembly, a forced mode without its
inputs, a wrong mode, transcript evidence without proteins, and a missing
`result` or `model_cfg` stop the run early (N5, N7, P9, P10);
`min_alignment_rate` is a parameter (N5); same-named BAMs and read files no
longer collide (N6); the dead `memory '180 GB'` is gone (N8);
`tiberius.model_dir` (N10); downloads carry the label `download` and run on
the submitting host (N11); `singularity.cacheDir` defaults to
`~/.cache/paludamentum/singularity` (N12); the HPC template needs no
`includeConfig` and `docs/hpc.md` says so (N13); retries only after a memory
or time kill, with memory × attempt (N14); `params.tools.*` is used
everywhere (N16); the `grep -c` count is right (N17); `restart`,
`prothint_conflict_filter`, `MERGE_GENEFINDER_TRAIN_PRIO`, `EMPTY_TSV`,
`OUT_CH` and the commented `STRINGTIE_MERGE` are removed (N18); the
`local.config` sizes tasks to the machine (P3); `parameters.yaml` works with
only the genome filled in (P11); the defaults tables agree (P12);
`STRINGTIE_ASSEMBLE_MIX` sorts its BAMs (D11).

Scripts: FASTA headers with descriptions keep their intron hints (B1); the
`STRG.` crash path is gone with the intrinsic stage (B2); genome chunks are
balanced under the `max_files` cap (B3); `hc.gff3` keeps gene and transcript
lines (B4); a gene nested in an intron is its own gene (B5); the same CDS
with and without UTRs is one transcript (B6); species ranking only for
OrthoDB ids, full database otherwise or when no species passes (B7); gene and
mRNA lines have distinct IDs (B8); an empty merge is an empty GFF3, not an
error (B9); `--transdecoder_util` is honoured, dead code removed, modules
without shebang are not executable, `generate_hc_genes.py` logs and has a
main guard, docstring names fixed (NITs). Drusilla branch: in-frame stop
check in the stop/start fix (D4), empty feature table guard (D15), batch
mode removed (D16), private paths and names removed from docstrings and docs
(L2, L3), `run: "on"` and `tiberius.result` without model in `hcMethod`
(D10, part), first script tests (D12, part).

Fixed 2026-10-01 (`main`): the overnight fixes are committed (`29a00c4`) and
pushed, and the first CI run with the new jobs (Python 3.9, Nextflow 25.04.0)
passed (L9). Drusilla branch: the hint rescue runs in its own image,
`gaiusaugustus/paludamentum-hint-rescue:0.1.0`, built from
`docker/hint_rescue/Dockerfile` and pushed (`ecde717`, `3513075`).

Fixed 2026-10-01 (drusilla worktree): the released LightGBM model is the
default `drusilla.lgb_model` (URL, downloaded once by Nextflow, sha256 in
`drusilla.lgb_model_sha256`) (D1); `apply_lgb_model_gtf.py` loads the
LightGBM text model and its JSON instead of the joblib pickle, checks the
sha256 of the archive and of the text model, the classes and the feature
count, and stops on a feature missing from the table instead of filling it
with 0 (D5); Artistic License 1.0 and copyright lines in the headers of the
six scripts from tiberius_orf_finder; the README describes the model and its
download (D13, part).

Tests: `tests/test_bin_scripts.py` (9), 15 new launcher tests, 10 new stub
runs, `tests/test_drusilla_scripts.py` (7). CI runs all of them (L9).

2026-10-01/03 (L7, D2): the Drusilla image of the flow is
`gaiusaugustus/drusilla:0.1.0` (on Docker Hub, accepted by Katharina on
2026-10-03), built from
`docker/drusilla/Dockerfile` (Drusilla `e3c5cf5`, lightgbm 4.7.0, Keras 3,
NGC TensorFlow stack guarded); the `drusilla` submodule points at the same
commit; `version_mismatches()` and `test_container_tags_of_base_config` cover
Drusilla.

---

## Suggested order

1. Commit the drusilla worktree (D1, D5 and the script fixes). Decisions 3
   and 4.
2. Merge `drusilla` into `main` (see "Merge notes"), fix D3, D6, D8, D9,
   D10, D13. Run CI through a PR.
3. Run T. rubripes and Bos taurus again (decision 4); test `tiberius.model_dir`
   on a GPU node (N10).
4. Get tags upstream and pin the submodules (L6), bump to v0.4.0 (L8),
   confirm the authors (L10).
5. Remove or rewrite `MIGRATION_PLAN.md` (L2). Tag, then make the repo
   public (L1).
