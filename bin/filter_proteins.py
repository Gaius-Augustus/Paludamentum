#!/usr/bin/env python3
"""
Copy a protein FASTA for BUSCO/compleasm/OMArk, dropping records that
hmmsearch refuses and, optionally, cleaning stop-codon characters.

  - records with >= --max-length residues (default 100000) are dropped:
    hmmsearch (p7_pipeline.c) refuses targets of 100,000 aa or more, and one
    such record stops the whole BUSCO/compleasm run (BRAKER4
    scripts/compleasm_wrapper.py drops them the same way). A trailing '*' is
    not counted as a residue.
  - with --strip-stop, a trailing '*' is removed and every internal '*' or
    '.' is replaced by 'X' (gffread -y normally drops the terminal stop codon
    already; internal stops or a '*' from other tools would otherwise reach
    hmmsearch).

The number and IDs of the dropped records go to stderr. Sequences are written
on one line each; headers are kept unchanged.

Usage:
    filter_proteins.py --in proteins.fa --out filtered.fa [--max-length 100000] [--strip-stop]
"""
from __future__ import annotations

import argparse
import sys
from typing import Iterator, TextIO

MAX_PROTEIN_LEN = 100000


def read_records(handle: TextIO) -> Iterator[tuple[str, str]]:
    """Yield (header line without '>', sequence) of a FASTA stream."""
    header = None
    seq: list[str] = []
    for line in handle:
        line = line.rstrip("\r\n")
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq)
            header = line[1:]
            seq = []
        elif header is not None:
            seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq)


def clean_stops(seq: str) -> str:
    """Remove a trailing '*' and turn internal '*' and '.' into 'X'."""
    seq = seq.rstrip("*")
    return seq.replace("*", "X").replace(".", "X")


def filter_proteins(fin: TextIO, fout: TextIO, max_length: int = MAX_PROTEIN_LEN,
                    strip_stop: bool = False) -> tuple[int, list[str]]:
    """Copy fin to fout; return (records written, IDs of dropped records)."""
    written = 0
    dropped: list[str] = []
    for header, seq in read_records(fin):
        if len(seq.rstrip("*")) >= max_length:
            dropped.append(header.split()[0] if header.split() else "")
            continue
        if strip_stop:
            seq = clean_stops(seq)
        fout.write(f">{header}\n{seq}\n")
        written += 1
    return written, dropped


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--in", dest="infile", required=True, help="protein FASTA")
    parser.add_argument("--out", required=True, help="filtered protein FASTA")
    parser.add_argument("--max-length", type=int, default=MAX_PROTEIN_LEN,
                        help="drop records with at least this many residues (default %(default)s)")
    parser.add_argument("--strip-stop", action="store_true",
                        help="remove a trailing '*', replace internal '*' and '.' by 'X'")
    args = parser.parse_args(argv)
    if args.max_length < 1:
        parser.error("--max-length must be positive")

    with open(args.infile) as fin, open(args.out, "w") as fout:
        written, dropped = filter_proteins(fin, fout, args.max_length, args.strip_stop)
    if dropped:
        print(f"filter_proteins: dropped {len(dropped)} record(s) with >= {args.max_length} residues: "
              + ", ".join(dropped), file=sys.stderr)
    else:
        print(f"filter_proteins: no record with >= {args.max_length} residues", file=sys.stderr)
    print(f"filter_proteins: wrote {written} record(s) to {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
