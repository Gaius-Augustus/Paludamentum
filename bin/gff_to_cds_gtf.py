#!/usr/bin/env python3
"""
Write the transcript and CDS lines of a GFF3 or GTF annotation as GTF lines
with gene_id and transcript_id. The LightGBM filter of the Drusilla flow reads
transcripts from the CDS lines of a GTF, the hint rescue from the transcript
lines.

Usage: gff_to_cds_gtf.py annotation.gff3 > annotation.cds.gtf
"""

import sys


def parse_attributes(text):
    attrs = {}
    for part in text.strip().strip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        else:
            key, _, value = part.partition(" ")
        attrs[key.strip()] = value.strip().strip('"')
    return attrs


def main(path):
    tx_gene = {}
    tx_lines = {}
    cds = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            attrs = parse_attributes(cols[8])
            if cols[2] in ("mRNA", "transcript") and "ID" in attrs:
                tx_gene[attrs["ID"]] = attrs.get("Parent", attrs["ID"])
                tx_lines[attrs["ID"]] = cols
            elif cols[2] in ("mRNA", "transcript") and "transcript_id" in attrs:
                tx_lines[attrs["transcript_id"]] = cols
            elif cols[2] == "CDS":
                cds.append((cols, attrs))
    out = sys.stdout
    written = set()
    for cols, attrs in cds:
        tx = attrs.get("transcript_id") or attrs.get("Parent", "").split(",")[0]
        if not tx:
            continue
        gene = attrs.get("gene_id") or tx_gene.get(tx, tx)
        gtf_attrs = f'gene_id "{gene}"; transcript_id "{tx}";'
        if tx not in written:
            written.add(tx)
            tx_cols = tx_lines.get(tx)
            if tx_cols is not None:
                out.write("\t".join(tx_cols[:2] + ["transcript"] + tx_cols[3:8] + [gtf_attrs]) + "\n")
        out.write("\t".join(cols[:8] + [gtf_attrs]) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
