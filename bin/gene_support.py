#!/usr/bin/env python3
"""
Evidence support of every coding transcript of a GFF3 by the hints of
hintsfile.gff: how many CDS introns and CDS segments are supported by RNA-Seq
hints, by protein hints and by any hint.

Support rules:
  - intron: the intron between two CDS segments of the transcript matches an
    `intron` hint exactly (seqid, strand, start, end)
  - CDS segment ("exon" in the column names): it overlaps a `CDSpart`,
    `exonpart`, `ep`, `CDS` or `exon` hint on the same seqid by at least one
    base, on the same strand or a hint with strand '.'
  - an intron hint with strand '.' matches both strands

One TSV row per coding transcript (mRNA with CDS) in file order; transcripts
without CDS are skipped. Single-exon transcripts have 0 introns and NA in the
intron percentage columns. A summary of all transcripts is written as '#'
lines above the column header; a one-line summary goes to stderr. An empty
hints file gives zero support everywhere.

Hint sources (the src= attribute) of Paludamentum's hintsfile.gff, which
CONCAT_HINTS (modules/util.nf) concatenates from the protein hints, the
RNA-Seq intron hints and the Iso-Seq intron hints:

  src  class    written by                                       features
  ---  -------  -----------------------------------------------  ---------------
  E    RNA-Seq  bam2hints --intronsonly + filterIntronsFindStrand  intron
                (BAM2HINTS, modules/rnaseq.nf) on the merged
                short-read BAM (HISAT2) and, as BAM2HINTS_ISO, on
                the merged Iso-Seq BAM (minimap2); pyVARUS
                hints.gff (made the same way); several sources
                summed by merge_intron_hints.py (src=E kept)
  P    protein  aln2hints.pl --prg=miniprot --priority=4          CDSpart, intron
                (ALN2HINTS, modules/proteins.nf) on miniprot.gtf
                of miniprothint (scorer2gtf.py of all scored
                miniprot alignments); no --genome_file, so no
                start/stop hints
  W, b2h  RNA-Seq  not written by Paludamentum; accepted as in BRAKER4
  PH      protein  not written by Paludamentum; accepted as in BRAKER4
  C       both     not written by Paludamentum (ETP combined hints);
                   counted as RNA-Seq and protein, as in BRAKER4
  other   -        ignored (counted on stderr)

"RNA-Seq" therefore means transcript evidence: short reads and Iso-Seq give
the same src=E and cannot be told apart in hintsfile.gff. miniprothint's
high-confidence hints (miniprot/hc.gff: intron, start_codon, stop_codon with
al_score/splice_sites/prots/CDS_overlap attributes and no src=) are not part
of hintsfile.gff; they are used by the Drusilla flow only. An intron or CDS
segment supported by both an RNA-Seq and a protein hint (or by a src=C hint)
counts in the RNA-Seq, the protein and the any columns; the summary lines
also count it as "both".

Port of BRAKER4 scripts/gene_support_summary.py (commit 3535ed3) to GFF3;
changes: GFF3 input read with gff3_lib (gene ID = Parent of the mRNA); only
transcripts with CDS, introns and "exons" taken from the CDS segments
(BRAKER4: exon features, CDS if there are none); hint matching strand-aware
(BRAKER4 ignores the strand); sorted lists + bisect instead of intervaltree
(standard library only); summary lines "Supported by both" added; one-line
summary on stderr; options --gff3/--hints/--out.

Usage:
    gene_support.py --gff3 final.gff3 --hints hintsfile.gff --out gene_support.tsv
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from bisect import bisect_right
from collections import defaultdict
from typing import Dict, List, Set, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gff3_lib import GFF3Error, open_text, read_gff3  # noqa: E402

RNASEQ_SOURCES = {"E", "W", "b2h"}
PROTEIN_SOURCES = {"P", "PH"}
COMBINED_SOURCES = {"C"}
EXON_HINT_TYPES = {"exonpart", "CDSpart", "ep", "exon", "CDS"}
CLASSES = ("rnaseq", "protein", "combined")

COLUMNS = ["gene_id", "transcript_id", "chrom", "strand",
           "num_introns",
           "introns_sup_rnaseq", "pct_introns_sup_rnaseq",
           "introns_sup_protein", "pct_introns_sup_protein",
           "introns_sup_any", "pct_introns_sup_any",
           "num_exons",
           "exons_sup_rnaseq", "pct_exons_sup_rnaseq",
           "exons_sup_protein", "pct_exons_sup_protein",
           "exons_sup_any", "pct_exons_sup_any"]

_SRC = re.compile(r"(?:^|;)\s*src=([^;]+)")


def source_class(attributes: str) -> str:
    match = _SRC.search(attributes)
    src = match.group(1).strip() if match else ""
    if src in COMBINED_SOURCES:
        return "combined"
    if src in RNASEQ_SOURCES:
        return "rnaseq"
    if src in PROTEIN_SOURCES:
        return "protein"
    return "other"


class OverlapIndex:
    """Intervals sorted by start with the running maximum of their ends: 'is there any
    interval overlapping start..end' by one bisect."""

    def __init__(self, intervals: List[Tuple[int, int]]):
        intervals = sorted(intervals)
        self.starts = [s for s, _e in intervals]
        self.max_end: List[int] = []
        top = 0
        for _s, e in intervals:
            top = max(top, e)
            self.max_end.append(top)

    def overlaps(self, start: int, end: int) -> bool:
        idx = bisect_right(self.starts, end)          # intervals starting at or before `end`
        return idx > 0 and self.max_end[idx - 1] >= start


class Hints:
    def __init__(self):
        # (seqid, strand, start, end) -> source classes
        self.introns: Dict[Tuple[str, str, int, int], Set[str]] = defaultdict(set)
        # (seqid, strand, class) -> OverlapIndex of the exon-like hints
        self.exon_index: Dict[Tuple[str, str, str], OverlapIndex] = {}
        self.n_intron = 0
        self.n_exon = 0
        self.n_other = 0

    def intron_classes(self, seqid: str, strand: str, start: int, end: int) -> Set[str]:
        return self.introns.get((seqid, strand, start, end), set()) | \
            self.introns.get((seqid, ".", start, end), set())

    def exon_classes(self, seqid: str, strand: str, start: int, end: int) -> Set[str]:
        found = set()
        for cls in CLASSES:
            for hint_strand in (strand, "."):
                index = self.exon_index.get((seqid, hint_strand, cls))
                if index is not None and index.overlaps(start, end):
                    found.add(cls)
                    break
        return found


def parse_hints(path: str) -> Hints:
    hints = Hints()
    exon_hints: Dict[Tuple[str, str, str], List[Tuple[int, int]]] = defaultdict(list)
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n\r").split("\t")
            if len(cols) < 9:
                continue
            feature = cols[2]
            if feature != "intron" and feature not in EXON_HINT_TYPES:
                continue
            seqid, strand = cols[0], cols[6].strip() or "."
            start, end = int(cols[3]), int(cols[4])
            if start > end:
                start, end = end, start
            cls = source_class(cols[8])
            if cls == "other":
                hints.n_other += 1
                continue
            if feature == "intron":
                hints.introns[(seqid, strand, start, end)].add(cls)
                hints.n_intron += 1
            else:
                exon_hints[(seqid, strand, cls)].append((start, end))
                hints.n_exon += 1
    hints.exon_index = {key: OverlapIndex(ivs) for key, ivs in exon_hints.items()}
    return hints


def classify(classes: Set[str]) -> Tuple[bool, bool]:
    """(RNA-Seq support, protein support) of a set of source classes."""
    combined = "combined" in classes
    return ("rnaseq" in classes or combined), ("protein" in classes or combined)


def cds_introns(segments: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    introns = []
    for (_s1, e1), (s2, _e2) in zip(segments, segments[1:]):
        if e1 + 1 <= s2 - 1:
            introns.append((e1 + 1, s2 - 1))
    return introns


def pct_cell(num: int, denom: int) -> str:
    return f"{100 * num / denom:.1f}" if denom > 0 else "NA"


def pct(num: int, denom: int) -> str:
    return f"{100 * num / denom:.1f}%" if denom > 0 else "NA"


def support(ann, hints: Hints):
    """Rows of the TSV and the totals."""
    rows = []
    tot = defaultdict(int)
    for gene, tx in ann.transcripts():
        if not tx.is_coding:
            continue
        seqid, strand = tx.feature.seqid, tx.feature.strand
        segments = [(c.start, c.end) for c in tx.cds]
        introns = cds_introns(segments)
        counts = defaultdict(int)
        for kind, intervals, lookup in (("i", introns, hints.intron_classes),
                                        ("e", segments, hints.exon_classes)):
            for start, end in intervals:
                rna, prot = classify(lookup(seqid, strand, start, end))
                counts[kind + "_rna"] += rna
                counts[kind + "_prot"] += prot
                counts[kind + "_any"] += rna or prot
                counts[kind + "_both"] += rna and prot
        n_i, n_e = len(introns), len(segments)
        rows.append([gene.id, tx.id, seqid, strand,
                     n_i, counts["i_rna"], pct_cell(counts["i_rna"], n_i),
                     counts["i_prot"], pct_cell(counts["i_prot"], n_i),
                     counts["i_any"], pct_cell(counts["i_any"], n_i),
                     n_e, counts["e_rna"], pct_cell(counts["e_rna"], n_e),
                     counts["e_prot"], pct_cell(counts["e_prot"], n_e),
                     counts["e_any"], pct_cell(counts["e_any"], n_e)])
        tot["introns"] += n_i
        tot["exons"] += n_e
        for key, value in counts.items():
            tot[key] += value
        tot["tx"] += 1
        if n_i == 0:
            tot["single_exon"] += 1
        elif counts["i_any"] == n_i:
            tot["fully"] += 1
        elif counts["i_any"] > 0:
            tot["some"] += 1
        else:
            tot["none"] += 1
    return rows, tot


def write_tsv(path: str, rows, tot) -> None:
    n_tx, n_i, n_e = tot["tx"], tot["introns"], tot["exons"]
    with open(path, "w", encoding="utf-8") as out:
        out.write("# Gene support summary\n")
        out.write(f"# Transcripts: {n_tx}\n")
        out.write(f"#   Single-exon (no introns to evaluate): {tot['single_exon']} ({pct(tot['single_exon'], n_tx)})\n")
        out.write(f"#   All introns supported: {tot['fully']} ({pct(tot['fully'], n_tx)})\n")
        out.write(f"#   Some introns supported: {tot['some']} ({pct(tot['some'], n_tx)})\n")
        out.write(f"#   No intron support: {tot['none']} ({pct(tot['none'], n_tx)})\n")
        out.write(f"# Introns: {n_i}\n")
        out.write(f"#   Supported by RNA-Seq: {tot['i_rna']} ({pct(tot['i_rna'], n_i)})\n")
        out.write(f"#   Supported by protein: {tot['i_prot']} ({pct(tot['i_prot'], n_i)})\n")
        out.write(f"#   Supported by both: {tot['i_both']} ({pct(tot['i_both'], n_i)})\n")
        out.write(f"#   Supported by any: {tot['i_any']} ({pct(tot['i_any'], n_i)})\n")
        out.write(f"# Exons: {n_e}\n")
        out.write(f"#   Supported by RNA-Seq: {tot['e_rna']} ({pct(tot['e_rna'], n_e)})\n")
        out.write(f"#   Supported by protein: {tot['e_prot']} ({pct(tot['e_prot'], n_e)})\n")
        out.write(f"#   Supported by both: {tot['e_both']} ({pct(tot['e_both'], n_e)})\n")
        out.write(f"#   Supported by any: {tot['e_any']} ({pct(tot['e_any'], n_e)})\n")
        out.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            out.write("\t".join(str(x) for x in row) + "\n")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff3", required=True, help="gene set (GFF3)")
    parser.add_argument("--hints", required=True, help="hints file (hintsfile.gff), may be empty")
    parser.add_argument("--out", required=True, help="output TSV")
    args = parser.parse_args(argv)
    try:
        ann = read_gff3(args.gff3)
    except GFF3Error as exc:
        print(f"gene_support.py: {args.gff3}: {exc}", file=sys.stderr)
        return 1
    try:
        hints = parse_hints(args.hints)
    except ValueError as exc:
        print(f"gene_support.py: {args.hints}: {exc}", file=sys.stderr)
        return 1
    rows, tot = support(ann, hints)
    write_tsv(args.out, rows, tot)
    other = f"; {hints.n_other} hints of other sources ignored" if hints.n_other else ""
    print(f"gene_support.py: {tot['tx']} coding transcripts; introns supported {tot['i_any']}/{tot['introns']} "
          f"({pct(tot['i_any'], tot['introns'])}; RNA-Seq {tot['i_rna']}, protein {tot['i_prot']}, "
          f"both {tot['i_both']}); CDS segments supported {tot['e_any']}/{tot['exons']} "
          f"({pct(tot['e_any'], tot['exons'])}){other}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
