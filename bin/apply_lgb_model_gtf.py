#!/usr/bin/env python3
# From tiberius_orf_finder/scripts/apply_lgb_model_gtf.py (Lars Gabriel). Changed in
# Paludamentum: cds_length_nt added to the features (it was filled with 0), the batch
# mode over a directory of species removed, no scoring of an empty feature table, the
# native LightGBM text model of the release instead of the joblib pickle, and a
# missing feature is an error instead of a column of 0.
# Copyright (c) 2026 Lars Gabriel. Artistic License 1.0, see LICENSE.
"""Score the ORF features of gene finder transcripts with the saved LightGBM
model of the Drusilla flow and filter the GTF.

Inputs:
  1. the model: the released archive (drusilla_lgb_3class_v1.tar.gz), its
     unpacked directory, or the .txt model with its .json next to it. The
     .txt is a LightGBM text model (Booster.save_model); the .json lists the
     feature names in model order, the classes and the sha256 of the .txt.
  2. the orf_features.tsv of compute_orf_features.py
  3. the GTF whose transcripts are scored (CDS lines with transcript_id)

The feature matrix follows the feature order of the .json. A feature that
cannot be built from the feature table is an error.

Filtering criterion
-------------------
A transcript is KEPT when:
    P(correct) + P(partial) >= --threshold   (default 0.5)
Use --score-col prob_correct to apply a stricter correct-only filter.
Every kept line gets lgb_class and lgb_prob_* attributes; the scores of all
transcripts go to <out-gtf without suffix>.scores.tsv.

Usage
-----
apply_lgb_model_gtf.py \\
  --model     drusilla_lgb_3class_v1.tar.gz \\
  --sha256    <sha256 of the archive> \\
  --features  orf_features.tsv \\
  --in-gtf    ab_initio.gtf \\
  --out-gtf   ab_initio.lgb_scored.gtf
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import tarfile
from pathlib import Path

import lightgbm as lgb
import numpy as np
import pandas as pd


TID_RE = re.compile(r'transcript_id\s+"([^"]+)"')

NUMERIC_FEATURES = [
    "n_exons", "cds_length_nt", "dist_upstream_stop_nt", "n_upstream_atgs",
    "n_overlapping_alignments", "best_identity", "best_norm_bitscore",
    "best_target_coverage", "best_protein_coverage", "frac_introns_supported",
    "cds_length_pct", "n_overlapping_alignments_pct",
    "protein_extends_5prime_codons", "protein_extends_3prime_codons",
]
BINARY_FEATURES = [
    "has_protein_support", "has_start_hint", "has_stop_hint",
    "has_conflict", "has_upstream_partner", "has_downstream_partner",
]
CATEGORICAL_FEATURES = {
    "lorf_class":    ["LORF_UPSTOP", "LORF_NOUPSTOP", "sORF_UPSTOP", "sORF_NOUPSTOP", "upLORF"],
    "support_level": ["fullSupport", "anySupport", "noSupport"],
}


def build_feature_matrix(df: pd.DataFrame, expected_names: list[str]) -> np.ndarray:
    """Build the same feature matrix layout the model was trained on."""
    parts, names = [], []
    for col in NUMERIC_FEATURES:
        if col in df.columns:
            parts.append(pd.to_numeric(df[col], errors="coerce").fillna(0).rename(col))
            names.append(col)
    for col in BINARY_FEATURES:
        if col in df.columns:
            parts.append(pd.to_numeric(df[col], errors="coerce").fillna(0).rename(col))
            names.append(col)
    for col, cats in CATEGORICAL_FEATURES.items():
        if col not in df.columns:
            continue
        for cat in cats:
            name = f"{col}__{cat}"
            parts.append((df[col] == cat).astype(float).rename(name))
            names.append(name)

    X = pd.concat(parts, axis=1)

    # training feature order; a missing feature would silently score as 0
    missing = [n for n in expected_names if n not in X.columns]
    if missing:
        sys.exit(f"ERROR: the feature table has no column for {len(missing)} model "
                 f"feature(s): {', '.join(missing)}")
    return X[expected_names].values


CLASS_LABELS = {0: "wrong", 1: "partial", 2: "correct"}


def sha256_of(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_model_files(path: Path) -> tuple[bytes, bytes, str]:
    """(text model, json, name) of a .tar.gz, a directory, or a .txt with its .json."""
    if path.is_file() and tarfile.is_tarfile(path):
        # read the members, nothing is extracted
        with tarfile.open(path) as tar:
            files = [m for m in tar.getmembers() if m.isfile()]
            jsons = [m for m in files if m.name.endswith(".json")]
            if len(jsons) != 1:
                sys.exit(f"ERROR: {path} must hold exactly one .json, it holds {len(jsons)}")
            meta_bytes = tar.extractfile(jsons[0]).read()
            model_file = json.loads(meta_bytes)["model_file"]
            txts = [m for m in files if Path(m.name).name == model_file]
            if len(txts) != 1:
                sys.exit(f"ERROR: {path} has no {model_file} next to its .json")
            return tar.extractfile(txts[0]).read(), meta_bytes, path.name
    if path.is_dir():
        jsons = sorted(path.glob("*.json"))
        if len(jsons) != 1:
            sys.exit(f"ERROR: {path} must hold exactly one .json, it holds {len(jsons)}")
        js = jsons[0]
    elif path.suffix == ".txt":
        js = path.with_suffix(".json")
    else:
        sys.exit(f"ERROR: {path} is not a .tar.gz, a directory or a .txt model "
                 "(the joblib pickle is not supported any more)")
    if not js.is_file():
        sys.exit(f"ERROR: {js} not found")
    meta_bytes = js.read_bytes()
    txt = js.parent / json.loads(meta_bytes)["model_file"]
    if not txt.is_file():
        sys.exit(f"ERROR: {txt} not found")
    return txt.read_bytes(), meta_bytes, path.name


def load_model(path: Path, sha256: str | None = None) -> tuple[lgb.Booster, list[str], dict]:
    """Booster, feature names and metadata of the model at path.

    sha256 is checked against the archive (or the .txt) given as path; the
    sha256 of the text model recorded in the .json is always checked."""
    if sha256:
        got = sha256_of(path.read_bytes()) if path.is_file() else None
        if got != sha256.lower():
            sys.exit(f"ERROR: sha256 of {path} is {got}, expected {sha256}. "
                     "For another model, set drusilla.lgb_model_sha256 to its sha256, or to null.")
    model_bytes, meta_bytes, name = read_model_files(path)
    meta = json.loads(meta_bytes)
    if sha256_of(model_bytes) != meta.get("model_sha256"):
        sys.exit(f"ERROR: the text model in {name} does not match the model_sha256 of its .json")
    classes = {int(k): v for k, v in meta.get("classes", {}).items()}
    if classes != CLASS_LABELS:
        sys.exit(f"ERROR: the model classes are {classes}, expected {CLASS_LABELS}")
    feat_names = list(meta["feature_names"])
    booster = lgb.Booster(model_str=model_bytes.decode())
    if booster.num_feature() != len(feat_names):
        sys.exit(f"ERROR: the model has {booster.num_feature()} features, "
                 f"its .json lists {len(feat_names)}")
    if booster.num_model_per_iteration() != len(CLASS_LABELS):
        sys.exit(f"ERROR: the model has {booster.num_model_per_iteration()} classes, "
                 f"expected {len(CLASS_LABELS)}")
    return booster, feat_names, meta


def score_features(features_tsv: Path, booster, feat_names: list[str],
                   threshold: float, score_col: str
                   ) -> tuple[set[str], pd.DataFrame, dict[str, str]]:
    """Returns (passing_tids, scores_df, tid→gtf_attr_string)."""
    df = pd.read_csv(features_tsv, sep="\t", low_memory=False)
    if df.empty:
        # no gene finder transcript: nothing to score, nothing passes
        print("  No transcripts in the feature table.", file=sys.stderr)
        return set(), pd.DataFrame(columns=["transcript_id"]), {}
    X  = build_feature_matrix(df, feat_names)
    proba = booster.predict(X)   # (n, 3) class probabilities

    scores = df[["transcript_id"]].copy()
    scores["prob_wrong"]     = proba[:, 0]
    scores["prob_partial"]   = proba[:, 1]
    scores["prob_correct"]   = proba[:, 2]
    scores["prob_not_wrong"] = proba[:, 1] + proba[:, 2]
    pred_class = np.argmax(proba, axis=1)
    scores["lgb_class"] = [CLASS_LABELS[c] for c in pred_class]

    # build per-transcript GTF attribute string
    attr_map: dict[str, str] = {}
    for row in scores.itertuples(index=False):
        attr_map[row.transcript_id] = (
            f' lgb_class "{row.lgb_class}";'
            f' lgb_prob_correct "{row.prob_correct:.4f}";'
            f' lgb_prob_partial "{row.prob_partial:.4f}";'
            f' lgb_prob_wrong "{row.prob_wrong:.4f}";'
        )

    passing_tids = set(scores.loc[scores[score_col] >= threshold, "transcript_id"])
    return passing_tids, scores, attr_map


def filter_gtf(in_gtf: Path, passing_tids: set[str],
               attr_map: dict[str, str], out_gtf: Path) -> tuple[int, int]:
    """Write kept GTF lines, appending LGB attributes to the 9th column."""
    kept, total = set(), set()
    with in_gtf.open() as fin, out_gtf.open("w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            m = TID_RE.search(line)
            if not m:
                continue
            tid = m.group(1)
            total.add(tid)
            if tid in passing_tids:
                kept.add(tid)
                # append attributes before the trailing newline
                extra = attr_map.get(tid, "")
                fout.write(line.rstrip("\n") + extra + "\n")
    return len(kept), len(total)


def process_one(features_tsv: Path, in_gtf: Path, out_gtf: Path,
                booster, feat_names: list[str], threshold: float,
                score_col: str, label: str) -> None:
    passing, scores, attr_map = score_features(features_tsv, booster, feat_names,
                                               threshold, score_col)
    n_kept, n_total = filter_gtf(in_gtf, passing, attr_map, out_gtf)
    pct = 100 * n_kept / n_total if n_total else 0
    print(f"  {label:<40} {n_kept:>6}/{n_total:<6} ({pct:.1f}%)  → {out_gtf.name}")

    scores_out = out_gtf.with_suffix(".scores.tsv")
    scores.to_csv(scores_out, sep="\t", index=False)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--model",     required=True, type=Path,
                   help="released model archive (.tar.gz), its unpacked directory, "
                        "or the .txt text model with its .json next to it")
    p.add_argument("--sha256", default=None,
                   help="expected sha256 of --model (the archive)")
    p.add_argument("--threshold", type=float, default=0.5)
    p.add_argument("--score-col", default="prob_not_wrong",
                   choices=["prob_not_wrong", "prob_correct", "prob_partial"],
                   help="Column to threshold on (default: prob_not_wrong = P(partial)+P(correct))")
    p.add_argument("--features", type=Path, required=True,
                   help="orf_features.tsv of compute_orf_features.py")
    p.add_argument("--in-gtf",   type=Path, required=True, help="GTF to score and filter")
    p.add_argument("--out-gtf",  type=Path, required=True, help="Filtered GTF with lgb_* attributes")
    args = p.parse_args()

    print(f"Loading model from {args.model} …", flush=True)
    booster, feat_names, meta = load_model(args.model, args.sha256)
    print(f"  {meta.get('name', args.model.name)}: {len(feat_names)} features, "
          f"threshold={args.threshold} on '{args.score_col}'", flush=True)

    if not args.in_gtf.exists():
        sys.exit(f"ERROR: --in-gtf {args.in_gtf} not found")
    print(f"\nScoring {args.features.name} …", flush=True)
    process_one(args.features, args.in_gtf, args.out_gtf,
                booster, feat_names, args.threshold, args.score_col,
                args.features.parent.name)
    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
