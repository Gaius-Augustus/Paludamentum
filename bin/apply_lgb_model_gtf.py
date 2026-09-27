#!/usr/bin/env python3
# From tiberius_orf_finder/scripts/apply_lgb_model_gtf.py (Lars Gabriel), with one fix: cds_length_nt
# is a feature of lgb_3class_model.pkl but was missing below and filled with 0.
"""Score a new orf_features.tsv with a saved LGB model and filter the GTF.

This is the stand-alone inference script. It requires only:
  1. A trained model saved by train_orf_lgb_3class.py  (lgb_3class_model.pkl)
  2. An orf_features.tsv for the target species/run
  3. The corresponding ORF GTF to filter

No re-training needed. The model and feature names are loaded from the pkl.

Filtering criterion
-------------------
A transcript is KEPT when:
    P(correct) + P(partial) >= --threshold   (default 0.5)
Use --score-col prob_correct to apply a stricter correct-only filter.

Usage
-----
# Single species
python scripts/apply_lgb_model_gtf.py \\
  --model     /projects/AI-GUSTUS/tiberius_orf_finder/results/filter_analysis/lgb_3class_model.pkl \\
  --features  <species_dir>/annotate_.../orf_features.tsv \\
  --in-gtf    <species_dir>/annotate_.../orfs.filtered.gtf \\
  --out-gtf   <species_dir>/annotate_.../orfs_lgb3_filtered.gtf

# All species under a base-dir (same layout as training)
python scripts/apply_lgb_model_gtf.py \\
  --model     /projects/AI-GUSTUS/tiberius_orf_finder/results/filter_analysis/lgb_3class_model.pkl \\
  --base-dir  /projects/AI-GUSTUS/tiberius_orf_finder/results/vertebrates_test \\
  --annot-tag annotate_epoch_74_filt_tpm1cov3len300_lorf \\
  --in-gtf    orfs.filtered.gtf \\
  --out-gtf   orfs_lgb3_filtered.gtf
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import joblib
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

    # align to training feature order; fill missing columns with 0
    missing = [n for n in expected_names if n not in X.columns]
    if missing:
        print(f"  Warning: {len(missing)} features not in input, filling with 0: {missing[:5]}…",
              file=sys.stderr)
    for m in missing:
        X[m] = 0.0
    return X[expected_names].values


CLASS_LABELS = {0: "wrong", 1: "partial", 2: "correct"}


def score_features(features_tsv: Path, clf, feat_names: list[str],
                   threshold: float, score_col: str
                   ) -> tuple[set[str], pd.DataFrame, dict[str, str]]:
    """Returns (passing_tids, scores_df, tid→gtf_attr_string)."""
    df = pd.read_csv(features_tsv, sep="\t", low_memory=False)
    X  = build_feature_matrix(df, feat_names)
    proba = clf.predict_proba(X)   # (n, 3) for 3-class or (n, 2) for binary

    scores = df[["transcript_id"]].copy()
    if proba.shape[1] == 3:
        scores["prob_wrong"]     = proba[:, 0]
        scores["prob_partial"]   = proba[:, 1]
        scores["prob_correct"]   = proba[:, 2]
        scores["prob_not_wrong"] = proba[:, 1] + proba[:, 2]
        pred_class = np.argmax(proba, axis=1)
        scores["lgb_class"] = [CLASS_LABELS[c] for c in pred_class]
    else:
        scores["prob_wrong"]     = proba[:, 0]
        scores["prob_correct"]   = proba[:, 1]
        scores["prob_not_wrong"] = proba[:, 1]
        scores["lgb_class"] = np.where(proba[:, 1] >= 0.5, "correct", "wrong")

    # build per-transcript GTF attribute string
    attr_map: dict[str, str] = {}
    for row in scores.itertuples(index=False):
        tid = row.transcript_id
        if proba.shape[1] == 3:
            attr = (
                f' lgb_class "{row.lgb_class}";'
                f' lgb_prob_correct "{row.prob_correct:.4f}";'
                f' lgb_prob_partial "{row.prob_partial:.4f}";'
                f' lgb_prob_wrong "{row.prob_wrong:.4f}";'
            )
        else:
            attr = (
                f' lgb_class "{row.lgb_class}";'
                f' lgb_prob_correct "{row.prob_correct:.4f}";'
                f' lgb_prob_wrong "{row.prob_wrong:.4f}";'
            )
        attr_map[tid] = attr

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
                clf, feat_names: list[str], threshold: float,
                score_col: str, label: str) -> None:
    passing, scores, attr_map = score_features(features_tsv, clf, feat_names,
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
                   help="lgb_3class_model.pkl saved by train_orf_lgb_3class.py")
    p.add_argument("--threshold", type=float, default=0.5)
    p.add_argument("--score-col", default="prob_not_wrong",
                   choices=["prob_not_wrong", "prob_correct", "prob_partial"],
                   help="Column to threshold on (default: prob_not_wrong = P(partial)+P(correct))")

    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--base-dir", type=Path,
                      help="Root with one subdir per species (batch mode)")
    mode.add_argument("--features", type=Path,
                      help="Single orf_features.tsv (single-species mode)")

    p.add_argument("--annot-tag", default="annotate_epoch_74_filt_tpm1cov3len300_lorf",
                   help="Subdir inside each species dir (batch mode only)")
    p.add_argument("--in-gtf",   default="orfs.filtered.gtf",
                   help="Input GTF filename (batch mode) or full path (single mode)")
    p.add_argument("--out-gtf",  default="orfs_lgb3_filtered.gtf",
                   help="Output GTF filename (batch mode) or full path (single mode)")
    args = p.parse_args()

    print(f"Loading model from {args.model} …", flush=True)
    clf, feat_names = joblib.load(args.model)
    print(f"  Features: {len(feat_names)}, threshold={args.threshold} on '{args.score_col}'",
          flush=True)

    if args.features:
        # single-species mode
        in_gtf  = Path(args.in_gtf)
        out_gtf = Path(args.out_gtf)
        if not in_gtf.exists():
            sys.exit(f"ERROR: --in-gtf {in_gtf} not found")
        print(f"\nScoring {args.features.name} …", flush=True)
        process_one(args.features, in_gtf, out_gtf,
                    clf, feat_names, args.threshold, args.score_col,
                    args.features.parent.name)
    else:
        # batch mode
        print(f"\n{'species':<42} {'kept/total':>14}  output")
        print("-" * 75)
        for sp_dir in sorted(args.base_dir.iterdir()):
            if not sp_dir.is_dir():
                continue
            tag_dir  = sp_dir / args.annot_tag
            feat_tsv = tag_dir / "orf_features.tsv"
            in_gtf   = tag_dir / args.in_gtf
            out_gtf  = tag_dir / args.out_gtf
            if not feat_tsv.exists() or not in_gtf.exists():
                continue
            process_one(feat_tsv, in_gtf, out_gtf,
                        clf, feat_names, args.threshold, args.score_col,
                        sp_dir.name)

    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
