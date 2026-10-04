#!/usr/bin/env python3
"""Merge intron hint files as if bam2hints had run on the merged BAM.

The inputs are the outputs of BAM2HINTS (bam2hints --intronsonly and
filterIntronsFindStrand.pl --score) or the hints.gff of pyVARUS runs, which
are made the same way. The multiplicities of an intron (seqname, start, end,
strand) are summed. The score column is the multiplicity, the attributes are
mult=N;pri=4;src=E (no mult= when N is 1), source, feature and frame are
taken from the input. Seqnames come in the order of their first appearance,
introns of a seqname are sorted by start, end and strand. Writes to stdout.
"""
import argparse
import sys


def parse_attrs(attr):
    """(multiplicity, the other attributes) of a hint attribute column."""
    mult = 1
    rest = []
    for part in attr.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if part.startswith("mult="):
            mult = int(part[5:])
        else:
            rest.append(part)
    return mult, rest


def merge(paths):
    seqnames = {}     # seqname -> {(start, end, strand): [mult, source, feature, frame, rest]}
    for path in paths:
        with open(path, encoding="utf-8") as handle:
            for number, line in enumerate(handle, 1):
                line = line.rstrip("\n")
                if not line.strip() or line.startswith("#"):
                    continue
                cols = line.split("\t")
                if len(cols) != 9:
                    sys.exit(f"{path}:{number}: expected 9 tab-separated columns, got {len(cols)}")
                mult, rest = parse_attrs(cols[8])
                key = (int(cols[3]), int(cols[4]), cols[6])
                introns = seqnames.setdefault(cols[0], {})
                if key in introns:
                    introns[key][0] += mult
                else:
                    introns[key] = [mult, cols[1], cols[2], cols[7], rest]
    return seqnames


def write(seqnames, out):
    for seqname, introns in seqnames.items():
        for (start, end, strand) in sorted(introns):
            mult, source, feature, frame, rest = introns[(start, end, strand)]
            attr = ";".join(([f"mult={mult}"] if mult > 1 else []) + rest)
            out.write(f"{seqname}\t{source}\t{feature}\t{start}\t{end}\t{mult}\t{strand}\t{frame}\t{attr}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("hints", nargs="+", help="intron hint files (GFF)")
    args = ap.parse_args()
    write(merge(args.hints), sys.stdout)


if __name__ == "__main__":
    main()
