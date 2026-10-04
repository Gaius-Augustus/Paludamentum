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
| L11 | This file mentions internal checks on brain. Move it to Paludamentum-UG or delete it before the repo goes public. | `RELEASE_AUDIT.md` |

---

## 2. Drusilla flow

### SHOULD

| # | Item | Where |
|---|------|-------|
| D9 | With `shards > 1` and `cache_dir: null`, all shards download and extract into the same `drusilla_cache` at once; the Drusilla registry `rmtree`s before it extracts, so the shards race. Download the model once before the fork. The weights URL has no sha256 in the manifest; which run the weights are (run009) is recorded nowhere. Nodes without internet fail. | `modules/drusilla.nf:56`, `drusilla/model_cfg/vertebrates.yaml` |

---

## 3. Launcher, packaging, configs, docs

### NIT

- `params.yaml` is written before validation, so a failed `--dry_run` leaves
  it behind (`paludamentum/launcher.py`). Validate before writing.

---

## 4. Nextflow workflow and configs

No open items.

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

1. Fix D9 on `main`. Run CI through a PR.
2. Run T. rubripes and Bos taurus again ("Before the release").
3. Move this file out (L11). Tag, then make the repo public (L1).
