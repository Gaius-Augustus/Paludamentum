# Paludamentum pre-release audit

Open items only; fixed items are deleted, the git history has them. Last
checked 2026-10-03 on `main` (v0.4.0, after the Drusilla merge, PR #2).

The submodules were only checked where Paludamentum depends on them.

---

## Before the release

1. **Benchmark again** (T. rubripes and Bos taurus): the annotation merge
   (exon-overlap clustering, CDS deduplication) and the in-frame stop check
   of the stop/start fix change the output. Gene-level F1 should move little,
   but it must be measured. T. rubripes with the v0.4.0 images: 80.36
   (80.26 before). Update the numbers in README "Drusilla flow"
   (80.36 / 79.46) with the results.

---

## 1. Release logistics

| # | Item | Where |
|---|------|-------|
| L1 | The GitHub repo is still **private**. Make it public last. | GitHub |
| L6 | Submodules are not pinned to releases, although README "Submodule pinning" says they are. Tiberius' newest tag is `v2.0.7` (submodule at `v2.0.7-11-g9734138`, image 2.0.8); Vipsania and Drusilla have no tags. Tags requested 2026-10-03: [Tiberius#120](https://github.com/Gaius-Augustus/Tiberius/issues/120) (`v2.0.8`), [Vipsania#8](https://github.com/Gaius-Augustus/Vipsania/issues/8) (`v1.0.1`; 1.0.1 is only on the branches `translation_table` and `blackwell-translation-table`, ; the Docker Hub tags `1.0.1` and `latest` held Vipsania 1.0.0 until 2026-10-03 and now point at the 1.0.1 build `1.0.1-tt`), [Drusilla#1](https://github.com/Gaius-Augustus/Drusilla/issues/1) (`v0.1.0` on `e3c5cf5`). When they exist, check out each tag in its submodule, set `gaiusaugustus/vipsania:1.0.1` in `conf/base.config`, and commit both. | `.gitmodules`, README "For maintainers" |
| L11 | This file mentions internal checks on brain. Move it to Paludamentum-UG or delete it before the repo goes public. | `RELEASE_AUDIT.md` |

---

## 2. Drusilla flow

### SHOULD

| # | Item | Where |
|---|------|-------|
| D4 | The splice sites of the extension exons of the stop/start fix are not checked. | `bin/fix_stop_by_miniprot.py` |
| D6 | The ORF-agreement filter of the hint rescue is switched off (benchmarked that way), but the module comments ("each with the hints of its best protein chain", "the predictions that agree with the hints") and the README sentence "hints of the best protein chain" still describe it. Loci without a chain run ab initio without hints. | `modules/drusilla.nf:202`, `:254`, README:282, `conf/base.config:103` |
| D8 | `DRUSILLA_ANNOTATE` and `HINT_RESCUE_TIBERIUS` have no memory label, so they get the 30 GB default (growing with the attempt). Measured: 40 GB for Bos, more than 90 GB for the rescue at seq_len 400k. Give them `bigmem`, or their own label. | `modules/drusilla.nf:40`, `:258`, `conf/base.config` |
| D9 | With `shards > 1` and `cache_dir: null`, all shards download and extract into the same `drusilla_cache` at once; the Drusilla registry `rmtree`s before it extracts, so the shards race. Download the model once before the fork. The weights URL has no sha256 in the manifest; which run the weights are (run009) is recorded nowhere. Nodes without internet fail. | `modules/drusilla.nf:56`, `drusilla/model_cfg/vertebrates.yaml` |
| D10 | Gaps in `hcMethod()`: <br>• Forcing `run: true` with `tiberius.result` and no model sets the rescue model to the string `"null"`. <br>• A custom model path counts as a vertebrate model by its file name (`/my/vertebrates.yaml`); read `target_species` from the YAML instead, or document it. | `subworkflows/drusilla.nf:43-44`, `lib_nf/functions.nf` drusillaModelEligible |
| D12 | Without tests: `compute_orf_features.py` (only its columns are tested), `filter_stringtie_gtf.py`, `prepare_hint_rescue_loci.py`, `filter_and_merge_rescue_gtf.py`, and stub runs of `mammalia*`, `Vertebrata`, `shards > 1`, `fix_stop: false`, `lgb_keep`, `tiberius.result` with Drusilla. The script tests need lightgbm and pandas, which CI does not install, so CI skips them. | `tests/test_drusilla_scripts.py`, `tests/test_stub_run.py` |

### NIT

| # | Item |
|---|------|
| D13 | The `drusilla:` config snippet in the README leaves out `batch_size` and `min_coding_length` (README:292-320). |
| D14 | `--proteins-fasta` is never passed, so `best_protein_coverage` is always 0. Check with Lars whether the model was trained with this feature (`modules/drusilla.nf`). |

---

## 3. Launcher, packaging, configs, docs

### NIT

- `params.yaml` is written before validation, so a failed `--dry_run` leaves
  it behind (`paludamentum/launcher.py`). Validate before writing.

---

## 4. Nextflow workflow and configs

### SHOULD

| # | Item | Where |
|---|------|-------|
| N10 | `tiberius.model_dir` stages the extracted weights into the task directory, where Tiberius' `download_weights()` finds `<model>_weights` and skips the download. Checked against the Tiberius code and in a stub run only: **run it once on a GPU node without internet** before documenting it as supported. | `modules/genefinder.nf`, `subworkflows/genefinder.nf`, `tiberius/tiberius/main.py:357-395` |

### NIT

- N19: the shell blocks do not quote file variables, so a path with
  whitespace in a params file run without the launcher breaks the run.

---

## 5. `bin/` scripts

### NIT

- `extend_cds_with_stop_codon.py`: a stop codon split by an intron extends
  the CDS across the intron. `check_stop_codons` removes these genes later,
  so the effect is a lost gene, not a wrong one.
- Upstream, in the Tiberius Dockerfile: TransDecoder master is cloned without
  a pin, next to v5.7.1.

---

## Suggested order

1. Fix D6, D8, D9, D10, D13 on `main`. Run CI through a PR.
2. Run T. rubripes and Bos taurus again ("Before the release"); test `tiberius.model_dir`
   on a GPU node (N10).
3. Pin the submodules to the requested tags (L6).
4. Move this file out (L11). Tag, then make the repo public (L1).
