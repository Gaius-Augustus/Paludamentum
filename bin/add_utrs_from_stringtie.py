#!/usr/bin/env python3
# Port of BRAKER4 scripts/stringtie2utr.py (commit 3535ed3), see the docstring.
# Copyright (c) 2025 Katharina Hoff. MIT License, see LICENSE-BRAKER4.
"""
Add UTRs from StringTie transcript assemblies to the coding transcripts of a
GFF3 that have none.

Port of BRAKER4 scripts/stringtie2utr.py (commit 3535ed3) to GFF3; changes:
GFF3 in and out instead of GTF; only coding transcripts without UTR features
are candidates (the original decorates every braker.gtf transcript, which
never has UTRs); all StringTie transcripts per intron are kept (the original
create_introns_hash keeps the last one only); standard library (dicts,
bisect, bins) instead of intervaltree; the CDS includes the stop codon, so the
3-bp three_prime_UTR artefact of the original cannot occur; several
--stringtie and --longread files; StringTie transcripts are told apart by
transcript_id, sequence and strand (IDs from per-scaffold runs may collide);
spliced length (sum of the exons) instead of exon plus intron length picks
the longest match; a StringTie exon that starts inside the CDS of a
single-exon gene no longer becomes a UTR that overlaps the CDS.

Rules (as in the original):
  - multi-exon transcript (CDS with introns): matches every StringTie
    transcript on the same sequence and strand whose introns include all its
    CDS introns
  - single-exon transcript: matches every StringTie transcript on the same
    strand that overlaps it, unless a StringTie intron overlaps its CDS
  - of several matches, the longest StringTie transcript (spliced length;
    ties: first seen); a match in a --longread file replaces a short-read
    match. StringTie transcripts with strand "." are ignored.
  - StringTie exons are clipped to CDS start - max_utr_extension and CDS end +
    max_utr_extension (0: no UTRs), and at the nearest same-strand transcript
    of another gene that has a match itself (ab initio neighbours without a
    match do not stop a UTR); for a single-exon transcript only StringTie
    exons that overlap its CDS are used; a StringTie exon that overlaps a CDS
    segment but is shorter than it is not used
  - the parts of the remaining StringTie exons outside the CDS span become
    five_prime_UTR (upstream in transcription direction) or three_prime_UTR
    (source stringtie2utr); the exons of the transcript become the union of
    CDS and UTRs. The CDS never changes. Gene and transcript spans follow.

Transcripts that already have UTRs, non-coding transcripts and transcripts
without a match are written unchanged. With neither --stringtie nor
--longread, the input is written unchanged.

The report has one row per coding transcript of the input: transcript_id,
matched_by (short, long, none, or has_utr for a transcript that already had
UTRs and was no candidate), and the bp of 5' and 3' UTR the transcript has in
the output; comment lines with the counts precede the header.

Usage:
    add_utrs_from_stringtie.py --gff3 in.gff3 [--stringtie a.gtf [b.gtf ...]] [--longread iso.gtf ...]
        --out out.gff3 [--max-utr-extension 5000] --report utr_report.tsv

Original: Katharina J. Hoff, University of Greifswald, 2023
(katharina.hoff@uni-greifswald.de). This program is free software; you can
redistribute it and/or modify it under the terms of the Artistic License.
"""
from __future__ import annotations

import argparse
import bisect
import os
import re
import sys
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gff3_lib import (UTR_TYPES, Annotation, Feature, GFF3Error, Transcript, merge_intervals,  # noqa: E402
                      open_text, read_gff3, sort_annotation, update_spans, write_gff3)

SOURCE = "stringtie2utr"
BIN = 100_000          # bin size of the overlap index of StringTie transcripts
_TID = re.compile(r'transcript_id "([^"]+)"')

Interval = Tuple[int, int]


# ---------------------------------------------------------------- StringTie

class Assemblies:
    """
    The StringTie transcripts of all input files, numbered in the order they
    are first seen (file by file), so that a lower index means seen first.
    """

    def __init__(self):
        self.exons: List[List[Interval]] = []
        self.locus: List[Tuple[str, str]] = []
        self.long: List[bool] = []
        self.length: List[int] = []

    def read(self, path: str, long: bool) -> None:
        index: Dict[Tuple[str, str, str], int] = {}
        with open_text(path) as handle:
            for number, raw in enumerate(handle, 1):
                if not raw.strip() or raw.startswith("#"):
                    continue
                cols = raw.rstrip("\n\r").split("\t")
                if len(cols) < 9:
                    raise GFF3Error(f"{path} line {number}: {len(cols)} tab-separated columns instead of 9")
                if cols[2] != "exon" or cols[6] not in ("+", "-"):
                    continue
                match = _TID.search(cols[8])
                if not match:
                    continue
                key = (match.group(1), cols[0], cols[6])
                i = index.get(key)
                if i is None:
                    i = index[key] = len(self.exons)
                    self.exons.append([])
                    self.locus.append((cols[0], cols[6]))
                    self.long.append(long)
                    self.length.append(0)
                start, end = int(cols[3]), int(cols[4])
                self.exons[i].append((start, end))
                self.length[i] += end - start + 1
        for i in index.values():
            self.exons[i].sort()

    def introns(self, i: int) -> List[Interval]:
        exons = self.exons[i]
        return [(a[1] + 1, b[0] - 1) for a, b in zip(exons, exons[1:]) if b[0] > a[1] + 1]

    def __len__(self) -> int:
        return len(self.exons)


# ---------------------------------------------------------------- candidates

class Candidate:
    __slots__ = ("tx", "gene_id", "seqid", "strand", "cds", "cds_start", "cds_end", "introns",
                 "span_start", "span_end")

    def __init__(self, tx: Transcript, gene_id: str):
        self.tx = tx
        self.gene_id = gene_id
        self.seqid = tx.feature.seqid
        self.strand = tx.feature.strand
        self.cds: List[Interval] = [(c.start, c.end) for c in tx.cds]
        self.cds_start = self.cds[0][0]
        self.cds_end = max(e for _s, e in self.cds)
        self.introns = [(a[1] + 1, b[0] - 1) for a, b in zip(self.cds, self.cds[1:]) if b[0] > a[1] + 1]
        self.span_start = min(c.start for c in tx.children)
        self.span_end = max(c.end for c in tx.children)


def overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start <= b_end and a_end >= b_start


def better(st: Assemblies, current: Optional[int], new: int) -> bool:
    """A long-read match replaces a short-read match; else the longer one, ties: first seen."""
    if current is None:
        return True
    if st.long[new] != st.long[current]:
        return st.long[new]
    return st.length[new] > st.length[current] or (st.length[new] == st.length[current] and new < current)


def find_matches(cands: List[Candidate], st: Assemblies) -> Dict[int, int]:
    """Candidate index -> index of the selected StringTie transcript."""
    multi = [k for k, c in enumerate(cands) if c.introns]
    single = [k for k, c in enumerate(cands) if not c.introns]
    best: Dict[int, int] = {}

    # multi-exon: intron chain subset; index only the introns some candidate has
    wanted = {(c.seqid, c.strand) + intron for c in (cands[k] for k in multi) for intron in c.introns}
    by_intron: Dict[Tuple[str, str, int, int], List[int]] = defaultdict(list)
    for i in range(len(st)):
        seqid, strand = st.locus[i]
        for intron in st.introns(i):
            key = (seqid, strand) + intron
            if key in wanted:
                by_intron[key].append(i)
    for k in multi:
        c = cands[k]
        lists = [by_intron.get((c.seqid, c.strand) + intron) for intron in c.introns]
        if not all(lists):
            continue
        lists.sort(key=len)
        common = set(lists[0])
        for other in lists[1:]:
            common.intersection_update(other)
            if not common:
                break
        for i in common:
            if better(st, best.get(k), i):
                best[k] = i

    # single-exon: overlap without a StringTie intron overlapping the CDS
    if single:
        bins: Dict[Tuple[str, str, int], List[int]] = defaultdict(list)
        for i in range(len(st)):
            seqid, strand = st.locus[i]
            exons = st.exons[i]
            for b in range(exons[0][0] // BIN, exons[-1][1] // BIN + 1):
                bins[(seqid, strand, b)].append(i)
        for k in single:
            c = cands[k]
            seen = set()
            for b in range(c.span_start // BIN, c.span_end // BIN + 1):
                for i in bins.get((c.seqid, c.strand, b), ()):
                    if i in seen:
                        continue
                    seen.add(i)
                    exons = st.exons[i]
                    if not overlaps(exons[0][0], exons[-1][1], c.span_start, c.span_end):
                        continue
                    if any(overlaps(s, e, c.cds_start, c.cds_end) for s, e in st.introns(i)):
                        continue
                    if better(st, best.get(k), i):
                        best[k] = i
    return best


def neighbour_bounds(cands: List[Candidate], matched: List[int]) -> Dict[int, Tuple[Optional[int], Optional[int]]]:
    """
    Region the UTRs of each matched candidate may occupy without running into
    a matched transcript of another gene on the same sequence and strand
    (neighbour_bounds of the original): (lo, hi), None where no such
    neighbour lies on that side. Overlapping transcripts are ignored.
    """
    by_locus: Dict[Tuple[str, str], List[Tuple[int, int, str]]] = defaultdict(list)
    for k in matched:
        c = cands[k]
        by_locus[(c.seqid, c.strand)].append((c.span_start, c.span_end, c.gene_id))
    index = {}
    for key, items in by_locus.items():
        by_end = sorted(items, key=lambda x: x[1])
        by_start = sorted(items, key=lambda x: x[0])
        index[key] = (by_end, [x[1] for x in by_end], by_start, [x[0] for x in by_start])

    bounds = {}
    for k in matched:
        c = cands[k]
        by_end, ends, by_start, starts = index[(c.seqid, c.strand)]
        lo = hi = None
        i = bisect.bisect_left(ends, c.span_start) - 1
        while i >= 0:
            if by_end[i][2] != c.gene_id:
                lo = by_end[i][1] + 1
                break
            i -= 1
        i = bisect.bisect_right(starts, c.span_end)
        while i < len(by_start):
            if by_start[i][2] != c.gene_id:
                hi = by_start[i][0] - 1
                break
            i += 1
        if lo is not None or hi is not None:
            bounds[k] = (lo, hi)
    return bounds


def clip(exons: List[Interval], lo: Optional[int], hi: Optional[int]) -> List[Interval]:
    out = []
    for s, e in exons:
        s = s if lo is None else max(s, lo)
        e = e if hi is None else min(e, hi)
        if s <= e:
            out.append((s, e))
    return out


def utr_pieces(c: Candidate, exons: List[Interval], max_ext: int,
               bound: Tuple[Optional[int], Optional[int]]) -> Tuple[List[Interval], List[Interval]]:
    """(5' UTR intervals, 3' UTR intervals) of a candidate from the exons of its StringTie match."""
    exons = clip(exons, max(1, c.cds_start - max_ext), c.cds_end + max_ext)
    exons = clip(exons, *bound)
    if not c.introns:
        exons = [x for x in exons if overlaps(x[0], x[1], c.cds_start, c.cds_end)]
    # a StringTie exon shorter than a CDS segment it overlaps is not used
    exons = [(s, e) for s, e in exons
             if not any(overlaps(s, e, cs, ce) and e - s < ce - cs for cs, ce in c.cds)]
    left: List[Interval] = []
    right: List[Interval] = []
    for s, e in exons:
        if s < c.cds_start:
            left.append((s, min(e, c.cds_start - 1)))
        if e > c.cds_end:
            right.append((max(s, c.cds_end + 1), e))
    return (left, right) if c.strand == "+" else (right, left)


def apply_utrs(c: Candidate, utr5: List[Interval], utr3: List[Interval]) -> None:
    """Replace the exons of the transcript by CDS + UTRs and add the UTR features."""
    tx = c.tx
    feat = tx.feature
    old_exons = tx.exons
    exon_source = old_exons[0].source if old_exons else feat.source
    tx.children = [ch for ch in tx.children if ch.type != "exon"]
    for kind, tag, pieces in (("five_prime_UTR", "utr5p", utr5), ("three_prime_UTR", "utr3p", utr3)):
        for n, (s, e) in enumerate(sorted(pieces), 1):
            tx.children.append(Feature(feat.seqid, SOURCE, kind, s, e, ".", feat.strand, ".",
                                       {"ID": [f"{tx.id}.{tag}{n}"], "Parent": [tx.id]}))
    for n, (s, e) in enumerate(merge_intervals(c.cds + utr5 + utr3), 1):
        tx.children.append(Feature(feat.seqid, exon_source, "exon", s, e, ".", feat.strand, ".",
                                   {"ID": [f"{tx.id}.exon{n}"], "Parent": [tx.id]}))


def cds_lines(ann: Annotation) -> Counter:
    """All CDS features with every column, to check that no CDS changes."""
    return Counter((ch.seqid, ch.source, ch.start, ch.end, ch.score, ch.strand, ch.phase,
                    tuple((k, tuple(v)) for k, v in ch.attrs.items()))
                   for _g, tx in ann.transcripts() for ch in tx.children if ch.type == "CDS")


def utr_bp(tx: Transcript) -> Tuple[int, int]:
    five = sum(u.length for u in tx.utrs if u.type in ("five_prime_UTR", "5UTR", "5'UTR"))
    three = sum(u.length for u in tx.utrs if u.type in ("three_prime_UTR", "3UTR", "3'UTR"))
    return five, three


def add_utrs(ann: Annotation, st: Assemblies, max_ext: int) -> Tuple[Dict[str, str], Dict[str, int]]:
    """Add UTRs in place; returns (transcript ID -> matched_by, counts)."""
    matched_by: Dict[str, str] = {}
    cands: List[Candidate] = []
    gene_of: Dict[int, object] = {}
    for gene, tx in ann.transcripts():
        if not tx.is_coding:
            continue
        if any(ch.type in UTR_TYPES for ch in tx.children):
            matched_by[tx.id] = "has_utr"
            continue
        matched_by[tx.id] = "none"
        if tx.feature.strand in ("+", "-"):
            gene_of[len(cands)] = gene
            cands.append(Candidate(tx, gene.id))

    best = find_matches(cands, st) if len(st) else {}
    matched = sorted(best)
    bounds = neighbour_bounds(cands, matched)
    counts = {"stopped": 0, "utr_added": 0}
    touched = {}
    for k in matched:
        c, i = cands[k], best[k]
        matched_by[c.tx.id] = "long" if st.long[i] else "short"
        lo, hi = bounds.get(k, (None, None))
        if any((lo is not None and s < lo) or (hi is not None and e > hi) for s, e in st.exons[i]):
            counts["stopped"] += 1
        utr5, utr3 = utr_pieces(c, st.exons[i], max_ext, (lo, hi))
        if utr5 or utr3:
            apply_utrs(c, utr5, utr3)
            counts["utr_added"] += 1
            touched[id(gene_of[k])] = gene_of[k]
    for gene in touched.values():
        update_spans(gene)
    return matched_by, counts


def write_report(path: str, ann: Annotation, matched_by: Dict[str, str], counts: Dict[str, int]) -> None:
    tally = Counter(matched_by.values())
    with open(path, "w", encoding="utf-8") as out:
        out.write(f"# coding transcripts: {len(matched_by)}\n")
        out.write(f"# with UTR before: {tally['has_utr']}\n")
        out.write(f"# matched short: {tally['short']}\n")
        out.write(f"# matched long: {tally['long']}\n")
        out.write(f"# none: {tally['none']}\n")
        out.write(f"# UTR added: {counts['utr_added']}\n")
        out.write(f"# UTR stopped at a matched neighbour gene: {counts['stopped']}\n")
        out.write("transcript_id\tmatched_by\tfive_prime_utr_bp\tthree_prime_utr_bp\n")
        for _gene, tx in ann.transcripts():
            if tx.id in matched_by:
                five, three = utr_bp(tx)
                out.write(f"{tx.id}\t{matched_by[tx.id]}\t{five}\t{three}\n")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff3", required=True, help="input GFF3")
    parser.add_argument("--stringtie", nargs="+", default=[], help="StringTie GTF(s) of short reads")
    parser.add_argument("--longread", nargs="+", default=[],
                        help="StringTie GTF(s) of long reads; a match here replaces a short-read match")
    parser.add_argument("--out", required=True, help="output GFF3")
    parser.add_argument("--max-utr-extension", type=int, default=5000,
                        help="bp a UTR may reach beyond the CDS on either side (default 5000; 0: no UTRs)")
    parser.add_argument("--report", required=True, help="report TSV")
    args = parser.parse_args(argv)
    if args.max_utr_extension < 0:
        parser.error("--max-utr-extension must not be negative")

    try:
        ann = read_gff3(args.gff3)
        st = Assemblies()
        for path in args.stringtie:
            st.read(path, long=False)
        for path in args.longread:
            st.read(path, long=True)
    except GFF3Error as exc:
        print(f"add_utrs_from_stringtie.py: {exc}", file=sys.stderr)
        return 1

    before = cds_lines(ann)
    matched_by, counts = add_utrs(ann, st, args.max_utr_extension)
    if cds_lines(ann) != before:
        raise RuntimeError("add_utrs_from_stringtie.py: a CDS changed; this is a bug")

    write_report(args.report, ann, matched_by, counts)
    sort_annotation(ann)
    with open(args.out, "w", encoding="utf-8") as out:
        write_gff3(ann, out)
    tally = Counter(matched_by.values())
    print(f"add_utrs_from_stringtie.py: {len(matched_by)} coding transcripts, {tally['has_utr']} with UTR before, "
          f"{tally['short']} matched short, {tally['long']} matched long, {tally['none']} none; "
          f"UTRs added to {counts['utr_added']} ({counts['stopped']} stopped at a matched neighbour gene); "
          f"{len(st)} StringTie transcripts", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
