#!/usr/bin/env python3
"""
Collect BUSCO and compleasm results of the genome and of the predicted
proteome into one TSV (qc/completeness.tsv):

    assessment  genome_or_proteome  complete  single  duplicated  fragmented  missing  n

one row per given input, in the order BUSCO genome, BUSCO proteome, compleasm
genome, compleasm proteome. Percentages are written as plain numbers (no '%'),
n as an integer. Every input is optional; a missing or empty file (the
pipeline passes placeholder files of size 0) is skipped.

Inputs:
  - BUSCO short summary (short_summary*.txt): the line
    'C:95.2%[S:93.1%,D:2.1%],F:1.2%,M:3.6%,n:255', with BUSCO 6 possibly
    followed by ',E:..%' (share of complete BUSCOs with internal stop codons,
    miniprot genome mode); E is written to a comment line, not to the table.
  - compleasm summary.txt: lines 'S:93.33%, 238', 'D:..', 'F:..', 'I:..',
    'M:..', 'N:255'. complete = S + D. compleasm splits fragmented genes in
    F (subclass 1: only part of the gene aligns) and I (subclass 2: parts of
    the gene align to different places of the assembly); both are fragments
    in the BUSCO sense, so 'fragmented' = F + I. Protein mode has no I line.

The lineage dataset ('# The lineage dataset is: ...' in BUSCO, '## lineage:
...' in compleasm) goes into a comment line '# lineage <assessment>
<genome_or_proteome>: <lineage>' before the header.

A non-empty input without a result line is an error (exit 1).

Usage:
    completeness_summary.py --out completeness.tsv [--busco-genome F] [--busco-proteins F]
                            [--compleasm-genome F] [--compleasm-proteins F]
"""
from __future__ import annotations

import argparse
import os
import re
import sys

HEADER = ["assessment", "genome_or_proteome", "complete", "single", "duplicated", "fragmented", "missing", "n"]

BUSCO_RE = re.compile(
    r"C:(?P<C>[\d.]+)%\[S:(?P<S>[\d.]+)%,D:(?P<D>[\d.]+)%\],F:(?P<F>[\d.]+)%,M:(?P<M>[\d.]+)%,n:(?P<n>\d+)"
    r"(?:,E:(?P<E>[\d.]+)%)?")
BUSCO_LINEAGE_RE = re.compile(r"^#\s*The lineage dataset is:\s*(\S+)")
COMPLEASM_LINE_RE = re.compile(r"^([SDFIM]):\s*([\d.]+)%\s*,\s*(\d+)\s*$")
COMPLEASM_N_RE = re.compile(r"^N:\s*(\d+)\s*$")
COMPLEASM_LINEAGE_RE = re.compile(r"^##\s*lineage:\s*(\S+)")


class SummaryError(Exception):
    pass


def fmt(value: float) -> str:
    """A percentage as a plain number: 2 decimals at most, no trailing zeros."""
    return f"{round(value, 2):g}"


def usable(path: str | None) -> bool:
    return bool(path) and os.path.isfile(path) and os.path.getsize(path) > 0


def parse_busco(path: str) -> dict:
    """Percentages, n, lineage and E (or None) of a BUSCO short summary."""
    lineage = None
    found = None
    with open(path) as fh:
        for line in fh:
            m = BUSCO_LINEAGE_RE.match(line.strip())
            if m and lineage is None:
                lineage = m.group(1)
            if found is None:
                found = BUSCO_RE.search(line)
    if found is None:
        raise SummaryError(f"{path}: no BUSCO result line 'C:..%[S:..%,D:..%],F:..%,M:..%,n:..'")
    g = found.groupdict()
    return {
        "complete": float(g["C"]), "single": float(g["S"]), "duplicated": float(g["D"]),
        "fragmented": float(g["F"]), "missing": float(g["M"]), "n": int(g["n"]),
        "lineage": lineage, "E": float(g["E"]) if g["E"] is not None else None,
    }


def parse_compleasm(path: str) -> dict:
    """Percentages, n and lineage of a compleasm summary.txt (fragmented = F + I)."""
    lineage = None
    pct: dict[str, float] = {}
    n = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            m = COMPLEASM_LINEAGE_RE.match(line)
            if m:
                lineage = m.group(1)
                continue
            m = COMPLEASM_LINE_RE.match(line)
            if m:
                pct[m.group(1)] = float(m.group(2))
                continue
            m = COMPLEASM_N_RE.match(line)
            if m:
                n = int(m.group(1))
    missing = [k for k in "SDFM" if k not in pct]
    if missing or n is None:
        need = ", ".join(missing + ([] if n is not None else ["N"]))
        raise SummaryError(f"{path}: no compleasm line(s) {need}")
    return {
        "complete": pct["S"] + pct["D"], "single": pct["S"], "duplicated": pct["D"],
        "fragmented": pct["F"] + pct.get("I", 0.0), "missing": pct["M"], "n": n,
        "lineage": lineage, "E": None,
    }


def summarize(inputs: list[tuple[str, str, str | None]]) -> tuple[list[str], list[list[str]]]:
    """inputs: (assessment, genome_or_proteome, path). Returns (comment lines, rows)."""
    comments: list[str] = []
    rows: list[list[str]] = []
    for assessment, kind, path in inputs:
        if not usable(path):
            if path:
                print(f"completeness_summary: {assessment} {kind}: {path} missing or empty, skipped",
                      file=sys.stderr)
            continue
        res = parse_busco(path) if assessment == "BUSCO" else parse_compleasm(path)
        if res["lineage"]:
            comments.append(f"# lineage {assessment} {kind}: {res['lineage']}")
        if res["E"] is not None:
            comments.append(f"# {assessment} {kind}: E (complete BUSCOs with internal stop codons) {fmt(res['E'])}%")
        rows.append([assessment, kind] + [fmt(res[k]) for k in HEADER[2:7]] + [str(res["n"])])
    return comments, rows


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", required=True, help="output TSV")
    parser.add_argument("--busco-genome", help="BUSCO short summary, genome mode")
    parser.add_argument("--busco-proteins", help="BUSCO short summary, protein mode")
    parser.add_argument("--compleasm-genome", help="compleasm summary.txt, genome mode")
    parser.add_argument("--compleasm-proteins", help="compleasm summary.txt, protein mode")
    args = parser.parse_args(argv)

    inputs = [("BUSCO", "genome", args.busco_genome), ("BUSCO", "proteome", args.busco_proteins),
              ("compleasm", "genome", args.compleasm_genome), ("compleasm", "proteome", args.compleasm_proteins)]
    try:
        comments, rows = summarize(inputs)
    except SummaryError as err:
        print(f"completeness_summary: {err}", file=sys.stderr)
        return 1
    with open(args.out, "w") as out:
        for line in comments:
            out.write(line + "\n")
        out.write("\t".join(HEADER) + "\n")
        for row in rows:
            out.write("\t".join(row) + "\n")
    print(f"completeness_summary: {len(rows)} row(s) written to {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
