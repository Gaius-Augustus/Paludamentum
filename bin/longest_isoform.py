#!/usr/bin/env python3
"""
Longest coding isoform per gene of a GFF3, for BUSCO and compleasm on the
proteome: one protein per gene avoids duplicates that are only alternative
isoforms. Per gene the mRNA with the largest summed CDS length is kept (the
first in the file on ties), with its children; genes without a coding
transcript are dropped. Port of GALBA2 scripts/get_longest_isoform.py to GFF3.

Usage:
    longest_isoform.py --gff3 in.gff3 --out out.gff3
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gff3_lib import Annotation, GFF3Error, read_gff3, update_spans, write_gff3  # noqa: E402


def longest_isoforms(ann: Annotation) -> Annotation:
    out = Annotation()
    for gene in ann.genes:
        coding = [tx for tx in gene.transcripts if tx.is_coding]
        if not coding:
            continue
        best = coding[0]
        for tx in coding[1:]:
            if tx.cds_length > best.cds_length:
                best = tx
        gene.transcripts = [best]
        update_spans(gene)
        out.genes.append(gene)
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff3", required=True, help="input GFF3")
    parser.add_argument("--out", required=True, help="output GFF3 with one transcript per gene")
    args = parser.parse_args(argv)
    try:
        ann = read_gff3(args.gff3)
    except GFF3Error as exc:
        print(f"longest_isoform.py: {args.gff3}: {exc}", file=sys.stderr)
        return 1
    n_tx = sum(len(g.transcripts) for g in ann.genes)
    out = longest_isoforms(ann)
    with open(args.out, "w", encoding="utf-8") as handle:
        write_gff3(out, handle)
    print(f"longest_isoform.py: {len(ann.genes)} genes, {n_tx} transcripts in; {len(out.genes)} transcripts out",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
