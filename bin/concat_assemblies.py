#!/usr/bin/env python3
"""One transcript set from a short-read and a long-read StringTie assembly.

The combined flow calls ORFs on ONE assembly. When the short reads and the
Iso-Seq reads of a mixed run do not both come as alignments (so that
`stringtie --mix` cannot assemble them together), their two assemblies are
kept side by side: this script writes the features of both files into one GTF
and prefixes the gene and transcript ids with `sr.` (short reads) and `lr.`
(long reads), because both StringTie runs number their genes from STRG.1.
Everything else is kept as it is, in particular the cov, FPKM and TPM
attributes that the StringTie pre-filter of Drusilla reads; `stringtie --merge`
would drop them (CONCAT_ASSEMBLIES in modules/assembly.nf).

A transcript that both assemblies contain is in the output twice. The ORF
callers then find the same CDS twice, and combine_genesets.py merges
identical CDS structures.

Usage: concat_assemblies.py --short short.gtf --long long.gtf > stringtie.gtf
"""
from __future__ import annotations

import argparse
import re
import sys

ID_ATTRIBUTE = re.compile(r'\b(gene_id|transcript_id) "([^"]*)"')


def prefixed(path: str, tag: str):
    """The feature lines of a GTF with tag + '.' in front of every gene_id and transcript_id."""
    features = 0
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                sys.exit(f"{path}: not a GTF line with 9 columns: {line.rstrip()[:120]}")
            if "transcript_id" not in fields[8]:
                sys.exit(f"{path}: a line without transcript_id; a StringTie assembly (GTF) is expected: "
                         f"{line.rstrip()[:120]}")
            fields[8] = ID_ATTRIBUTE.sub(lambda m: f'{m.group(1)} "{tag}.{m.group(2)}"', fields[8])
            features += 1
            yield "\t".join(fields) + "\n"
    if not features:
        sys.exit(f"{path}: no feature lines; a StringTie assembly (GTF) is expected.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--short", required=True, help="StringTie assembly of the short reads (GTF)")
    parser.add_argument("--long", required=True, help="StringTie assembly of the long reads (GTF)")
    args = parser.parse_args()
    sys.stdout.write("# short-read and long-read StringTie assemblies side by side (concat_assemblies.py)\n")
    for path, tag in ((args.short, "sr"), (args.long, "lr")):
        sys.stdout.writelines(prefixed(path, tag))


if __name__ == "__main__":
    main()
