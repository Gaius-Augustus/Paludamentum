#!/usr/bin/env python3
"""
Rewrite a GFF3 so that it satisfies the GFF3 contract of the published files
(docs/postprocessing.md, items 1 to 5) without changing any coordinate:

  - ##gff-version 3, nine columns, attribute values percent-encoded
  - gene -> transcript -> exon/CDS/UTR; a coding `transcript` becomes `mRNA`;
    start_codon, stop_codon and intron features are dropped; a coding
    transcript without exons gets exons (the union of its CDS and UTRs)
  - IDs: exons <tx>.exon<N>, CDS <tx>.cds (one ID for all segments), UTRs
    <tx>.utr5p<N> / <tx>.utr3p<N>, numbered in genome order; UTR types
    five_prime_UTR / three_prime_UTR by their side of the CDS
  - gene_biotype=protein_coding on genes with a coding transcript
  - CDS phases from the phase of the first segment and the segment lengths
  - gene and transcript spans recomputed, features sorted by seqid, gene
    start, transcript, start

Fails (exit 1) with a message naming the feature on a missing Parent, a
duplicated ID, a CDS or UTR outside the exons of its transcript, or
overlapping exons. Coding transcripts whose first CDS segment does not start
at phase 0 (5' partial) are counted on stderr; with --genome, also those
whose first codon is not ATG.

Usage:
    normalize_gff3.py --gff3 in.gff3 --out out.gff3 [--genome genome.fa]
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gff3_lib import (UTR3_TYPES, UTR5_TYPES, UTR_TYPES, Annotation, Feature, GFF3Error,  # noqa: E402
                      Transcript, contained, merge_intervals, phase_of, read_fasta, read_gff3, set_phases,
                      sort_annotation, spliced_cds, update_spans, write_gff3)


def utr_type(tx: Transcript, utr: Feature) -> str:
    """five_prime_UTR or three_prime_UTR by the side of the CDS the UTR lies on."""
    cds = tx.cds
    if not cds:
        if utr.type in UTR5_TYPES:
            return "five_prime_UTR"
        if utr.type in UTR3_TYPES:
            return "three_prime_UTR"
        return utr.type
    cds_start, cds_end = min(c.start for c in cds), max(c.end for c in cds)
    upstream = utr.end < cds_start if tx.feature.strand != "-" else utr.start > cds_end
    downstream = utr.start > cds_end if tx.feature.strand != "-" else utr.end < cds_start
    if upstream:
        return "five_prime_UTR"
    if downstream:
        return "three_prime_UTR"
    raise GFF3Error(f"{utr.describe()} of {tx.id} overlaps the CDS span {cds_start}-{cds_end}")


def normalize_transcript(tx: Transcript) -> None:
    feat = tx.feature
    for child in tx.children:
        if child.seqid != feat.seqid:
            raise GFF3Error(f"{child.describe()} is on another sequence than its transcript {tx.id}")
        if child.strand != feat.strand:
            raise GFF3Error(f"{child.describe()} is on another strand than its transcript {tx.id}")
    coding = tx.is_coding
    if coding and feat.type == "transcript":
        feat.type = "mRNA"

    for utr in [c for c in tx.children if c.type in UTR_TYPES]:
        utr.type = utr_type(tx, utr)
        utr.phase = "."

    exons = tx.exons
    if not exons and coding:
        for start, end in merge_intervals((c.start, c.end) for c in tx.children):
            tx.children.append(Feature(feat.seqid, feat.source, "exon", start, end, ".", feat.strand, ".", {}))
        exons = tx.exons
    for prev, nxt in zip(exons, exons[1:]):
        if nxt.start <= prev.end:
            raise GFF3Error(f"exons of {tx.id} overlap: {prev.describe()} and {nxt.describe()}")
    intervals = [(e.start, e.end) for e in exons]
    for child in tx.children:
        if child.type == "CDS" or child.type in UTR_TYPES:
            if not contained(child, intervals):
                raise GFF3Error(f"{child.describe()} of {tx.id} lies outside the exons of its transcript")

    # IDs of the contract
    for n, exon in enumerate(tx.exons, 1):
        exon.attrs = {"ID": [f"{tx.id}.exon{n}"], "Parent": [tx.id],
                      **{k: v for k, v in exon.attrs.items() if k not in ("ID", "Parent")}}
        exon.phase = "."
    for cds in tx.cds:
        cds.attrs = {"ID": [f"{tx.id}.cds"], "Parent": [tx.id],
                     **{k: v for k, v in cds.attrs.items() if k not in ("ID", "Parent")}}
    for kind, tag in (("five_prime_UTR", "utr5p"), ("three_prime_UTR", "utr3p")):
        for n, utr in enumerate(tx.of_type(kind), 1):
            utr.attrs = {"ID": [f"{tx.id}.{tag}{n}"], "Parent": [tx.id],
                         **{k: v for k, v in utr.attrs.items() if k not in ("ID", "Parent")}}
    set_phases(tx)


def normalize(ann: Annotation) -> dict:
    """Normalise in place; returns counts for the log."""
    counts = {"genes": 0, "transcripts": 0, "coding": 0, "first_phase_not_0": 0}
    for gene in ann.genes:
        for tx in gene.transcripts:
            normalize_transcript(tx)
            counts["transcripts"] += 1
            if tx.is_coding:
                counts["coding"] += 1
                if phase_of(tx.cds_in_transcription_order()[0]) != 0:
                    counts["first_phase_not_0"] += 1
        if any(tx.is_coding for tx in gene.transcripts):
            gene.feature.attrs["gene_biotype"] = ["protein_coding"]
        update_spans(gene)
        counts["genes"] += 1
    # genes without transcripts carry nothing of the annotation
    ann.genes = [g for g in ann.genes if g.transcripts]
    sort_annotation(ann)
    return counts


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff3", required=True, help="input GFF3")
    parser.add_argument("--out", required=True, help="output GFF3")
    parser.add_argument("--genome", help="genome FASTA: count coding transcripts whose first codon is not ATG")
    args = parser.parse_args(argv)
    try:
        ann = read_gff3(args.gff3)
        counts = normalize(ann)
    except GFF3Error as exc:
        print(f"normalize_gff3.py: {args.gff3}: {exc}", file=sys.stderr)
        return 1
    with open(args.out, "w", encoding="utf-8") as out:
        write_gff3(ann, out)
    print(f"normalize_gff3.py: {counts['genes']} genes, {counts['transcripts']} transcripts "
          f"({counts['coding']} coding)", file=sys.stderr)
    if counts["first_phase_not_0"]:
        print(f"normalize_gff3.py: {counts['first_phase_not_0']} coding transcripts do not start at phase 0 "
              "(5' partial)", file=sys.stderr)
    if args.genome:
        genome = read_fasta(args.genome, {g.feature.seqid for g in ann.genes})
        no_atg = sum(1 for _g, tx in ann.transcripts()
                     if tx.is_coding and phase_of(tx.cds_in_transcription_order()[0]) == 0
                     and spliced_cds(tx, genome)[:3] != "ATG")
        if no_atg:
            print(f"normalize_gff3.py: {no_atg} coding transcripts start with another codon than ATG",
                  file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
