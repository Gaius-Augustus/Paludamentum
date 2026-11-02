#!/usr/bin/env python3
"""Warn when many proteins of a gene set have an internal stop codon.

PROTEIN_FROM_GFF writes internal stops as * (gffread -S), which DIAMOND accepts,
so a gene set with wrong CDS phases no longer stops the run at DIAMOND. Its
proteins are mistranslated and lose their DIAMOND support, without any error.

Measured on the proteins of whole gene sets (October 2026): A. thaliana Tiberius
with wrong minus-strand phases 20.4%, the same set with repaired phases 0%. Sets
without such an error reached 2.9% (A. thaliana Drusilla), H. sapiens Tiberius
1.05%. MAX_PERCENT lies between them.

The warning goes to stdout, which the process forwards to the log of the run; no
other output is written, and the exit status stays 0.

Usage: check_internal_stops.py proteins.fa <name of the gene set>
"""

import sys

MAX_PERCENT = 5.0


def count(path):
    """Number of proteins, and of those with a stop before their last residue."""
    n = with_stop = 0
    seq = []

    def finish():
        nonlocal n, with_stop
        if seq:
            n += 1
            body = "".join(seq).rstrip("*.")
            if "*" in body or "." in body:
                with_stop += 1

    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                finish()
                seq = []
            else:
                seq.append(line.strip())
    finish()
    return n, with_stop


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    path, name = sys.argv[1:]
    n, with_stop = count(path)
    percent = 100 * with_stop / n if n else 0.0
    if percent > MAX_PERCENT:
        print(f"WARN: {percent:.1f}% of the {name} proteins ({with_stop} of {n}) have an internal stop codon. "
              f"Check the CDS phases of a reused {name}.result (bin/fix_cds_phases.py repairs them): "
              "with wrong phases the proteins are mistranslated and lose their DIAMOND support.")


if __name__ == "__main__":
    main()
