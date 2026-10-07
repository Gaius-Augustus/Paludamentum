#!/usr/bin/env python3
"""
Sanity filter of the final protein-coding annotation.

Per transcript, in this order:
  g. a transcript without CDS (`transcript`, or `mRNA` without CDS): removed, no_cds
  a. sum of the CDS lengths not a multiple of 3: removed, cds_length_mod3
  b. CDS segments that overlap, or lie on different strands or sequences:
     removed, cds_structure
  c. a CDS segment outside the exons of its transcript: removed, cds_outside_exon
  d. a stop codon before the last codon of the translated CDS: removed, internal_stop
  e. the CDS ends without a stop codon, and the next three bases are one:
     the last CDS segment (and its exon) is extended by 3 bp, kept, extended_stop;
     neither: kept, no_stop (a 3' partial gene). The three bases must lie in
     the exon of the CDS end when another exon follows it: a codon that reaches
     into the intron is not added (it would be split, and the GT donor site
     fakes TA|G).
  f. the first codon is not ATG: kept, non_atg_start (a note only)

A gene whose transcripts are all removed is removed. The result is
normalised to the GFF3 contract (normalize_gff3.py), so gene and transcript
spans follow the changes. The report is a TSV `gene_id transcript_id action
reason` (one row per transcript and reason) after comment lines with the counts.

--keep reports what would be removed or extended (action `flagged`) and
writes the normalised input without changes.

Usage:
    sanity_filter_gff3.py --gff3 in.gff3 --genome genome.fa --out out.gff3 --report sanity_filter.tsv [--keep]
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import Counter
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gff3_lib import (STOP_CODONS, UTR_TYPES, Annotation, GFF3Error, Transcript, contained,  # noqa: E402
                      phase_of, read_fasta, read_gff3, spliced_cds, subseq, translate, write_gff3)
from normalize_gff3 import normalize  # noqa: E402

REMOVE_REASONS = ("no_cds", "cds_length_mod3", "cds_structure", "cds_outside_exon", "internal_stop")
NOTE_REASONS = ("extended_stop", "no_stop", "non_atg_start")


def structure_problem(tx: Transcript) -> Optional[str]:
    """Reason of checks g, a, b, c, or None."""
    cds = tx.cds
    if not cds:
        return "no_cds"
    if sum(c.length for c in cds) % 3:
        return "cds_length_mod3"
    if len({c.seqid for c in cds}) > 1 or len({c.strand for c in cds}) > 1 or cds[0].strand not in "+-" \
            or cds[0].seqid != tx.feature.seqid or cds[0].strand != tx.feature.strand:
        return "cds_structure"
    if any(b.start <= a.end for a, b in zip(cds, cds[1:])):
        return "cds_structure"
    exons = [(e.start, e.end) for e in tx.exons]
    if exons and not all(contained(c, exons) for c in cds):
        return "cds_outside_exon"
    return None


def stop_extension(tx: Transcript, genome: Dict[str, str]) -> Optional[Tuple[int, int]]:
    """
    Genomic interval of the stop codon that follows a CDS without one, or None.
    When another exon follows the one that holds the CDS end, all three bases
    must lie inside that exon: a codon reaching into the intron would be split
    by it, and the intron's GT donor completes TA|G to a stop that is not there.
    A CDS whose exon is the last one may extend past the exon end (the exon
    grows with it).
    """
    cds = tx.cds
    strand = tx.feature.strand
    seq = genome[tx.feature.seqid]
    exons = tx.exons
    if strand == "+":
        last = cds[-1]
        start, end = last.end + 1, last.end + 3
        if end > len(seq):
            return None
        holder = [e for e in exons if e.start <= last.end <= e.end]
        followed = any(e.start > last.end for e in exons)
        split = followed and any(e.end < end for e in holder)
    else:
        last = cds[0]
        start, end = last.start - 3, last.start - 1
        if start < 1:
            return None
        holder = [e for e in exons if e.start <= last.start <= e.end]
        followed = any(e.end < last.start for e in exons)
        split = followed and any(e.start > start for e in holder)
    if split:
        return None
    return (start, end) if subseq(genome, tx.feature.seqid, start, end, strand) in STOP_CODONS else None


def extend_stop(tx: Transcript, interval: Tuple[int, int]) -> None:
    """Add the stop codon to the last CDS segment and its exon; UTRs give way."""
    start, end = interval
    strand = tx.feature.strand
    cds = tx.cds
    last = cds[-1] if strand == "+" else cds[0]
    old_start, old_end = last.start, last.end
    if strand == "+":
        last.end = end
    else:
        last.start = start
    for exon in tx.exons:
        if strand == "+" and exon.start <= old_end <= exon.end and exon.end < end:
            exon.end = end
        if strand == "-" and exon.start <= old_start <= exon.end and exon.start > start:
            exon.start = start
    kept = []
    for child in tx.children:
        if child.type in UTR_TYPES and child.start <= end and child.end >= start:
            # a 3' UTR that began right after the old CDS end
            if strand == "+":
                child.start = end + 1
            else:
                child.end = start - 1
            if child.start > child.end:
                continue
        kept.append(child)
    tx.children = kept


def check(tx: Transcript, genome: Dict[str, str], keep: bool) -> Tuple[Optional[str], List[str]]:
    """(removal reason or None, notes) of one transcript; extends the stop codon unless keep."""
    reason = structure_problem(tx)
    if reason:
        return reason, []
    if tx.feature.seqid not in genome:
        raise GFF3Error(f"sequence {tx.feature.seqid} of transcript {tx.id} is not in the genome")
    seq = spliced_cds(tx, genome)
    protein = translate(seq)
    if "*" in protein[:-1]:
        return "internal_stop", []
    notes = []
    if not protein.endswith("*"):
        interval = stop_extension(tx, genome)
        if interval:
            notes.append("extended_stop")
            if not keep:
                extend_stop(tx, interval)
        else:
            notes.append("no_stop")
    if phase_of(tx.cds_in_transcription_order()[0]) == 0 and seq[:3] != "ATG":
        notes.append("non_atg_start")
    return None, notes


def sanity_filter(ann: Annotation, genome: Dict[str, str], keep: bool = False) -> List[Tuple[str, str, str, str]]:
    """Filter in place; returns the report rows (gene_id, transcript_id, action, reason)."""
    rows = []
    for gene in ann.genes:
        kept = []
        for tx in gene.transcripts:
            reason, notes = check(tx, genome, keep)
            if reason:
                rows.append((gene.id, tx.id, "flagged" if keep else "removed", reason))
                if keep:
                    kept.append(tx)
            else:
                kept.append(tx)
            for note in notes:
                action = "flagged" if keep and note == "extended_stop" else "kept"
                rows.append((gene.id, tx.id, action, note))
        gene.transcripts = kept
    return rows


def write_report(path: str, rows, n_genes: int, n_tx: int, ann: Annotation, keep: bool) -> None:
    removed_tx = {r[1] for r in rows if r[3] in REMOVE_REASONS}
    left = {g.id for g in ann.genes if g.transcripts}
    by_reason = Counter(r[3] for r in rows)
    with open(path, "w", encoding="utf-8") as out:
        out.write(f"# sanity_filter_gff3.py{' --keep (nothing removed or extended)' if keep else ''}\n")
        out.write(f"# input: {n_genes} genes, {n_tx} transcripts\n")
        out.write(f"# transcripts {'flagged for removal' if keep else 'removed'}: {len(removed_tx)}; "
                  f"genes left: {len(left)}\n")
        out.write("# removal reasons: " + ", ".join(f"{r} {by_reason[r]}" for r in REMOVE_REASONS) + "\n")
        out.write("# notes on kept transcripts: " + ", ".join(f"{r} {by_reason[r]}" for r in NOTE_REASONS) + "\n")
        out.write("gene_id\ttranscript_id\taction\treason\n")
        for row in rows:
            out.write("\t".join(row) + "\n")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff3", required=True, help="input GFF3 (merge_annotations.py output)")
    parser.add_argument("--genome", required=True, help="genome FASTA")
    parser.add_argument("--out", required=True, help="filtered and normalised GFF3")
    parser.add_argument("--report", required=True, help="TSV of the removed and noted transcripts")
    parser.add_argument("--keep", action="store_true", help="report only; remove and extend nothing")
    args = parser.parse_args(argv)
    try:
        ann = read_gff3(args.gff3)
        n_genes = len(ann.genes)
        n_tx = sum(len(g.transcripts) for g in ann.genes)
        genome = read_fasta(args.genome, {g.feature.seqid for g in ann.genes})
        rows = sanity_filter(ann, genome, args.keep)
        write_report(args.report, rows, n_genes, n_tx, ann, args.keep)
        normalize(ann)
    except GFF3Error as exc:
        hint = " (the sanity filter is off: postprocess.sanity_filter = false)" if args.keep else ""
        print(f"sanity_filter_gff3.py: {args.gff3}: {exc}{hint}", file=sys.stderr)
        return 1
    with open(args.out, "w", encoding="utf-8") as out:
        write_gff3(ann, out)
    removed = len({r[1] for r in rows if r[3] in REMOVE_REASONS})
    print(f"sanity_filter_gff3.py: {n_tx} transcripts, {removed} {'flagged' if args.keep else 'removed'}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
