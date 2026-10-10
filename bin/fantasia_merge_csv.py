#!/usr/bin/env python3
"""
Join the CSV files of the FANTASIA_ANNOTATE tasks (results.csv or
failed_sequences.csv of each protein chunk) into one, in the order given.

The header of the first file and the rows of all files, line by line and
with their line ends (failed_sequences.csv ends lines with \\r\\n), as
fantasia_pipeline.py of FANTASIA-Lite joins the results of its own chunks
(append_csv): one row per protein, GO term and reference, nothing computed
across proteins. 300 A. thaliana proteins in one task and in three chunks
(CPU, brain, 2026-10-10) gave the same 1,032 rows; distances differed by at
most 3.6e-7 (float32 embeddings of other batches), the row order differed
(FANTASIA-Lite does not keep the input order). Empty files are skipped; a
file with another header stops the merge.

Usage:
    fantasia_merge_csv.py -o results.csv chunk1/results.csv chunk2/results.csv ...
"""

import argparse
import sys


def merge(paths, out):
    """Write the header of the first non-empty file and the rows of all files to out."""
    header = None
    for path in paths:
        with open(path, newline="") as fh:
            first = fh.readline()
            if not first:
                continue
            if header is None:
                header = first
                out.write(first)
            elif first != header:
                sys.exit(f"fantasia_merge_csv.py: {path}: header {first.rstrip()!r} differs from "
                         f"{header.rstrip()!r} of the first file")
            for line in fh:
                out.write(line)
    if header is None:
        sys.exit("fantasia_merge_csv.py: no rows: all input files are empty")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--out", required=True, help="joined CSV")
    parser.add_argument("csv", nargs="+", help="CSV files in chunk order")
    args = parser.parse_args()
    with open(args.out, "w", newline="") as out:
        merge(args.csv, out)


if __name__ == "__main__":
    main()
