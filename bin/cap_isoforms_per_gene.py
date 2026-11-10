#!/usr/bin/env python3
"""Cap the number of transcripts per gene in a GFF3 or GTF.

Reads the annotation (file or stdin), keeps at most --max transcripts per
gene (preferring the ones with the largest total CDS length), and writes the
reduced annotation to stdout in the format it came in. Non-feature lines
(`##`, blank) and gene features are preserved. Used before `gffread -y` in
`PROTEIN_FROM_GFF` to work around a gffread malloc() crash on genes with many
isoforms. The proteins go to DIAMOND: the species ranking of
PREPROCESS_PROTEINDB, and the four tables that combine_genesets.py reads (an
isoform beyond the cap has no protein there and so no DIAMOND hit).

GFF3 (Tiberius, Vipsania, TD2): the mRNA / transcript features name
themselves (ID) and their gene (Parent), their children name their
transcript (Parent). GTF (the Drusilla ORFs): every feature names its
transcript (transcript_id) and gene (gene_id).
"""
import argparse
import sys
from collections import defaultdict

TRANSCRIPT_TYPES = {"mRNA", "transcript"}


def parse_attrs(attr_str):
    """GFF3 (key=value) or GTF (key "value") attributes as a dict."""
    attrs = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k.strip()] = v.strip()
        else:
            toks = part.split(None, 1)
            attrs[toks[0]] = toks[1].strip().strip('"') if len(toks) > 1 else ""
    return attrs


def transcript_of(ftype, attrs):
    """(transcript, gene) of a feature. The gene is None where the feature does
    not name it (a child feature in GFF3); both are None for a gene feature
    or a feature without a transcript."""
    if ftype == "gene":
        return None, None
    if "transcript_id" in attrs:                      # GTF
        return attrs["transcript_id"], attrs.get("gene_id", "")
    if ftype in TRANSCRIPT_TYPES:                     # GFF3 transcript
        return attrs.get("ID"), attrs.get("Parent", "").split(",")[0]
    return attrs.get("Parent", "").split(",")[0] or None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max", type=int, default=10)
    ap.add_argument("gff", nargs="?", default="-")
    args = ap.parse_args()

    src = sys.stdin if args.gff == "-" else open(args.gff, "r", encoding="utf-8")
    all_lines = src.readlines()
    if src is not sys.stdin:
        src.close()

    tx_parent = {}
    tx_cds_len = defaultdict(int)
    for line in all_lines:
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) != 9:
            continue
        ftype = cols[2]
        tid, gid = transcript_of(ftype, parse_attrs(cols[8]))
        if tid is None:
            continue
        if gid is not None:
            tx_parent.setdefault(tid, gid)
        tx_cds_len.setdefault(tid, 0)
        if ftype == "CDS":
            tx_cds_len[tid] += int(cols[4]) - int(cols[3]) + 1

    by_gene = defaultdict(list)
    for tid, gid in tx_parent.items():
        by_gene[gid].append(tid)

    keep = set()
    for gid, tids in by_gene.items():
        tids_sorted = sorted(tids, key=lambda t: (-tx_cds_len[t], t))
        keep.update(tids_sorted[: args.max])

    for line in all_lines:
        if not line.strip() or line.startswith("#"):
            sys.stdout.write(line)
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) != 9:
            sys.stdout.write(line)
            continue
        ftype = cols[2]
        tid, _gid = transcript_of(ftype, parse_attrs(cols[8]))
        if ftype == "gene" or (tid is not None and tid in keep):
            sys.stdout.write(line)


if __name__ == "__main__":
    main()
