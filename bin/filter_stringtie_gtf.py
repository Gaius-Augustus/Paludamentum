#!/usr/bin/env python3
# From tiberius_orf_finder/scripts/filter_stringtie_gtf.py (Lars Gabriel), used unchanged in the Drusilla flow.
"""Pre-filter a StringTie GTF by transcript-level cov/TPM/length thresholds.

Drops whole transcripts (both their ``transcript`` line and all their
``exon`` lines). Other feature types are passed through unchanged. The
output GTF preserves StringTie's gene_id / transcript_id, so the
downstream ``annotate.py`` path is unaffected.

Filtering rule
--------------
A transcript is kept iff::

    length >= min_length
    AND  cov   >= min_cov
    AND  TPM   >= tpm_min(length)

where TPM is the StringTie ``TPM`` attribute, ``cov`` is the per-base
read coverage, and ``length`` is the sum of exon spans (inclusive). The
TPM threshold is length-aware to avoid penalising long, evenly-covered
transcripts that TPM (a per-kb quantity) naturally pushes downward::

    tpm_min(length) = min_tpm_long  if length >= long_length
    tpm_min(length) = min_tpm       otherwise

Defaults (TPM>=1, cov>=3, length>=300; relax to TPM>=0.5 when
length>=3000) match the discussion in plans/eval_observations.md for
the vertebrates_test precision pass.

Outputs
-------
* ``--out-gtf``     filtered GTF.
* ``--out-tsv``     per-transcript decision log
  (tx_id, gene_id, length, cov, tpm, kept, reason).
* stderr           kept/dropped counts and the breakdown by reason.

CLI::

    python scripts/filter_stringtie_gtf.py \\
        --in-gtf  stringtie.gtf \\
        --out-gtf stringtie.filtered.gtf \\
        --out-tsv stringtie.filtered.decisions.tsv
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path


def _parse_attr(attr_col: str, key: str) -> str | None:
    for chunk in attr_col.strip().strip(";").split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        parts = chunk.split(None, 1)
        if len(parts) == 2 and parts[0] == key:
            return parts[1].strip().strip('"')
    return None


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-gtf", type=Path, required=True)
    ap.add_argument("--out-gtf", type=Path, required=True)
    ap.add_argument("--out-tsv", type=Path, default=None,
                    help="Per-transcript decision log. Defaults to "
                         "<out-gtf>.decisions.tsv.")
    ap.add_argument("--min-length", type=int, default=300)
    ap.add_argument("--min-cov", type=float, default=3.0)
    ap.add_argument("--min-tpm", type=float, default=1.0)
    ap.add_argument("--long-length", type=int, default=3000,
                    help="Transcripts at or above this exon-summed length "
                         "get the relaxed TPM threshold.")
    ap.add_argument("--min-tpm-long", type=float, default=0.5,
                    help="Relaxed TPM threshold for long transcripts. TPM "
                         "is length-normalised, so long real transcripts "
                         "drift downward in TPM even at the same gene-level "
                         "expression as a short one.")
    return ap.parse_args(argv)


def _collect_transcripts(in_gtf: Path) -> dict[str, dict]:
    """Return ``{tx_id: {gene_id, length, cov, tpm}}`` for every transcript.

    ``length`` is summed exon span (inclusive). ``cov`` and ``tpm`` come
    from the ``transcript`` feature line (StringTie writes them there).
    A transcript missing from the dict's final view (no exons) is treated
    as length 0 by the caller.
    """
    info: dict[str, dict] = {}
    with in_gtf.open() as fh:
        for raw in fh:
            if not raw or raw.startswith("#"):
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            feature = cols[2]
            start = int(cols[3])
            end = int(cols[4])
            attrs = cols[8]
            tx_id = _parse_attr(attrs, "transcript_id")
            if tx_id is None:
                continue
            entry = info.setdefault(tx_id, {
                "gene_id": _parse_attr(attrs, "gene_id"),
                "length": 0,
                "cov": None,
                "tpm": None,
            })
            if entry["gene_id"] is None:
                entry["gene_id"] = _parse_attr(attrs, "gene_id")
            if feature == "exon":
                entry["length"] += (end - start + 1)
            elif feature == "transcript":
                cov = _parse_attr(attrs, "cov")
                tpm = _parse_attr(attrs, "TPM")
                if cov is not None:
                    try:
                        entry["cov"] = float(cov)
                    except ValueError:
                        pass
                if tpm is not None:
                    try:
                        entry["tpm"] = float(tpm)
                    except ValueError:
                        pass
    return info


def _decide(
    info: dict,
    min_length: int,
    min_cov: float,
    min_tpm: float,
    long_length: int,
    min_tpm_long: float,
) -> tuple[bool, str]:
    """Return ``(keep, reason)``. ``reason`` names the first failing rule
    or ``'pass'`` / ``'pass_long_relaxed'`` on success."""
    length = info["length"]
    cov = info["cov"]
    tpm = info["tpm"]
    if length < min_length:
        return False, "length"
    if cov is None or cov < min_cov:
        return False, "cov"
    is_long = length >= long_length
    tpm_threshold = min_tpm_long if is_long else min_tpm
    if tpm is None or tpm < tpm_threshold:
        return False, "tpm_long" if is_long else "tpm"
    return True, "pass_long_relaxed" if is_long else "pass"


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    out_tsv = args.out_tsv or args.out_gtf.with_suffix(args.out_gtf.suffix + ".decisions.tsv")

    args.out_gtf.parent.mkdir(parents=True, exist_ok=True)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    info = _collect_transcripts(args.in_gtf)

    keep: dict[str, str] = {}
    reasons: Counter[str] = Counter()
    for tx_id, entry in info.items():
        ok, reason = _decide(
            entry, args.min_length, args.min_cov, args.min_tpm,
            args.long_length, args.min_tpm_long,
        )
        reasons[reason] += 1
        if ok:
            keep[tx_id] = reason

    with out_tsv.open("w") as fh:
        fh.write("tx_id\tgene_id\tlength\tcov\ttpm\tkept\treason\n")
        for tx_id, entry in info.items():
            kept = "1" if tx_id in keep else "0"
            ok, reason = _decide(
                entry, args.min_length, args.min_cov, args.min_tpm,
                args.long_length, args.min_tpm_long,
            )
            fh.write(
                f"{tx_id}\t{entry['gene_id']}\t{entry['length']}\t"
                f"{entry['cov']}\t{entry['tpm']}\t{kept}\t{reason}\n"
            )

    kept_lines = 0
    total_lines = 0
    with args.in_gtf.open() as fin, args.out_gtf.open("w") as fout:
        for raw in fin:
            if raw.startswith("#"):
                fout.write(raw)
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 9:
                fout.write(raw)
                continue
            total_lines += 1
            tx_id = _parse_attr(cols[8], "transcript_id")
            if tx_id is None or tx_id in keep:
                fout.write(raw)
                kept_lines += 1

    n_total = len(info)
    n_kept = len(keep)
    n_dropped = n_total - n_kept
    pct = (100.0 * n_kept / n_total) if n_total else 0.0
    print(
        f"transcripts total={n_total} kept={n_kept} ({pct:.1f}%) "
        f"dropped={n_dropped}",
        file=sys.stderr,
    )
    for reason in ("pass", "pass_long_relaxed", "length", "cov", "tpm", "tpm_long"):
        if reasons[reason]:
            print(f"  {reason}: {reasons[reason]}", file=sys.stderr)
    print(
        f"gtf lines kept={kept_lines}/{total_lines}",
        file=sys.stderr,
    )
    print(f"out_gtf={args.out_gtf}", file=sys.stderr)
    print(f"out_tsv={out_tsv}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
