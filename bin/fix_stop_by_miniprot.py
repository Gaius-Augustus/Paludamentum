#!/usr/bin/env python3
# From tiberius_orf_finder/scripts/fix_stop_by_miniprot.py (Lars Gabriel). Changed in Paludamentum:
# a corrected CDS is only written if it has no in-frame stop codon before its end,
# every intron it adds has canonical splice sites, and a start fix begins with the ATG.
# Copyright (c) 2026 Lars Gabriel. Artistic License 1.0, see LICENSE.
"""Fix or recover stop codons in predicted ORFs using miniprot protein-to-genome alignments.

Two cases are handled:

  1. Partial ORFs (no stop codon predicted — transcript truncated at 3' end):
     Extend the CDS to a miniprot-supported stop codon on the genome.

  2. Complete ORFs with an early spurious stop codon:
     If a miniprot alignment extends past the predicted stop AND a valid
     downstream stop codon is found in the genome, replace the spurious
     stop with the protein-supported position.

For each qualifying miniprot alignment the stop codon is located by:
  a) Checking the triplet immediately after the last aligned CDS block.
  b) If not found there, scanning forward (+ strand) or backward (- strand)
     in-frame up to --max-stop-scan nt.  This recovers cases where miniprot
     did not mark StopCodon=1 (e.g., alignment ends at a contig edge).

The stop codon triplet is always verified in the genome FASTA before any
corrected output is written, and a corrected CDS is only written if its
translation has no stop codon before the last one (an extension that runs
through an in-frame stop is dropped, and the ORF is kept as it was, or
omitted if it was partial). Every intron a correction adds to the ORF must
be at least MIN_INTRON nt long with GT/GC...AG splice sites; otherwise the
next candidate alignment or hint is tried.

Usage
-----
python scripts/fix_stop_by_miniprot.py \\
    --orfs     results/orfs.gtf \\
    --partial  results/orfs.partial.gtf \\
    --miniprot results/miniprot.gff \\
    --genome   data/genome.fa \\
    --out      results/orfs.fixed.gtf

Output
------
One GTF file containing:
  - Complete ORFs corrected to a protein-supported stop (source unchanged).
  - Partial ORFs recovered with a protein-supported stop (source unchanged).
  - Complete ORFs for which no correction was needed, passed through as-is.
  - Partial ORFs with no protein support are silently omitted.

A summary is printed to stderr.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

from pyfaidx import Fasta

STOP_CODONS: frozenset[str] = frozenset({"TAA", "TAG", "TGA"})
DONORS: frozenset[str] = frozenset({"GT", "GC"})
MIN_INTRON = 20
# candidate corrections rejected by splice_sites_ok(), reported at the end
n_bad_splice = 0
_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _rev_comp(seq: str) -> str:
    return seq.translate(_RC)[::-1]


def _parse_attr(col: str, key: str) -> str | None:
    """Extract key from a GTF/GFF3 attribute column."""
    m = re.search(rf'{key}\s+"([^"]+)"', col)
    if m:
        return m.group(1)
    m = re.search(rf'{key}=([^;"\s]+)', col)
    if m:
        return m.group(1)
    return None


# ──────────────────────────────────────────────────────────────────────────────
# Data classes
# ──────────────────────────────────────────────────────────────────────────────

class ORF:
    """Predicted CDS parsed from a GTF file."""
    __slots__ = ("tid", "contig", "strand", "segments", "is_partial", "source",
                 "lorf_class")

    def __init__(self, tid: str, contig: str, strand: str,
                 is_partial: bool, source: str) -> None:
        self.tid = tid
        self.contig = contig
        self.strand = strand
        self.segments: list[tuple[int, int]] = []  # 0-based half-open, sorted ASC
        self.is_partial = is_partial
        self.source = source
        self.lorf_class: str | None = None

    @property
    def span_start(self) -> int:
        return self.segments[0][0]

    @property
    def span_end(self) -> int:
        return self.segments[-1][1]

    @property
    def coding_len(self) -> int:
        return sum(e - s for s, e in self.segments)


class MpAlignment:
    """Miniprot protein-to-genome CDS alignment."""
    __slots__ = ("mid", "contig", "strand", "cds_segments", "identity", "has_stop")

    def __init__(self, mid: str, contig: str, strand: str,
                 identity: float, has_stop: bool) -> None:
        self.mid = mid
        self.contig = contig
        self.strand = strand
        self.cds_segments: list[tuple[int, int]] = []  # 0-based half-open, sorted ASC
        self.identity = identity
        self.has_stop = has_stop  # whether miniprot marked StopCodon=1


# ──────────────────────────────────────────────────────────────────────────────
# Parsing
# ──────────────────────────────────────────────────────────────────────────────

def parse_orf_gtf(path: Path, is_partial: bool) -> dict[str, ORF]:
    """Parse a CDS-only GTF produced by annotate.py."""
    orfs: dict[str, ORF] = {}
    for row in csv.reader(open(path), delimiter="\t"):
        if not row or row[0].startswith("#"):
            continue
        if len(row) < 9 or row[2] != "CDS":
            continue
        tid = _parse_attr(row[8], "transcript_id")
        if tid is None:
            continue
        contig, strand, source = row[0], row[6], row[1]
        start = int(row[3]) - 1   # GTF 1-based inclusive → 0-based
        end = int(row[4])          # GTF inclusive → exclusive
        if tid not in orfs:
            orfs[tid] = ORF(tid, contig, strand, is_partial, source)
        orfs[tid].segments.append((start, end))
    for orf in orfs.values():
        orf.segments.sort()
        orf.lorf_class = None  # populated below during GTF re-scan
    # Second pass: pick up lorf_class (one attribute per transcript, any CDS line)
    for row in csv.reader(open(path), delimiter="\t"):
        if not row or row[0].startswith("#") or len(row) < 9 or row[2] != "CDS":
            continue
        tid = _parse_attr(row[8], "transcript_id")
        if tid and tid in orfs and orfs[tid].lorf_class is None:
            lc = _parse_attr(row[8], "lorf_class")
            if lc:
                orfs[tid].lorf_class = lc
    return orfs


def parse_miniprot_gff(path: Path) -> list[MpAlignment]:
    """Parse a miniprot GFF3, keeping all alignments that have CDS features.

    Unlike the previous version, StopCodon=1 is NOT required.  The flag is
    stored in MpAlignment.has_stop for informational purposes, but stop-codon
    verification is done independently by scan_for_stop().
    """
    pending: dict[str, MpAlignment] = {}
    for row in csv.reader(open(path), delimiter="\t"):
        if not row or row[0].startswith("#"):
            continue
        if len(row) < 9:
            continue
        feat = row[2]
        if feat == "mRNA":
            mid = _parse_attr(row[8], "ID")
            if mid is None:
                continue
            ident_s = _parse_attr(row[8], "Identity")
            identity = float(ident_s) if ident_s else 0.0
            has_stop = _parse_attr(row[8], "StopCodon") == "1"
            pending[mid] = MpAlignment(mid, row[0], row[6], identity, has_stop)
        elif feat == "CDS":
            parent = _parse_attr(row[8], "Parent")
            if parent not in pending:
                continue
            start = int(row[3]) - 1
            end = int(row[4])
            pending[parent].cds_segments.append((start, end))

    result = []
    for aln in pending.values():
        if aln.cds_segments:
            aln.cds_segments.sort()
            result.append(aln)
    return result


# ──────────────────────────────────────────────────────────────────────────────
# Overlap index and lookup
# ──────────────────────────────────────────────────────────────────────────────

MpIndex = dict[tuple[str, str], list[MpAlignment]]
HintIndex = dict[tuple[str, str], list[tuple[int, int]]]
StartIndex = dict[tuple[str, str], list[tuple[int, int]]]  # (contig,strand)->[(atg_pos,score)]


def parse_miniprothint_hints(path: Path) -> HintIndex:
    """Parse miniprothint hc.gff and return (contig, strand) → sorted intron list.

    Each intron is stored as (start, end) in 0-based half-open coordinates:
      start = left exon end  (first intronic base)
      end   = right exon start (first post-intronic base)
    Only 'intron' features are used; stop/start_codon hints are ignored.
    """
    hints: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    for row in csv.reader(open(path), delimiter="\t"):
        if not row or row[0].startswith("#") or len(row) < 9:
            continue
        if row[2] != "intron":
            continue
        s = int(row[3]) - 1   # GFF3 1-based → 0-based half-open start
        e = int(row[4])        # GFF3 inclusive end = 0-based exclusive end
        hints[(row[0], row[6])].append((s, e))
    return {k: sorted(v) for k, v in hints.items()}


def parse_miniprothint_starts(path: Path) -> StartIndex:
    """Parse start_codon features from a miniprothint GFF (hc.gff or miniprothint.gff).

    Returns {(contig, strand): sorted [(atg_pos, score), ...]} where atg_pos is the
    0-based genomic start of the ATG and score is the integer protein-support count
    (GFF score column).
    """
    starts: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    for row in csv.reader(open(path), delimiter="\t"):
        if not row or row[0].startswith("#") or len(row) < 9:
            continue
        if row[2] != "start_codon":
            continue
        atg_pos = int(row[3]) - 1   # 1-based → 0-based
        score = int(float(row[5])) if row[5] not in (".", "") else 0
        starts[(row[0], row[6])].append((atg_pos, score))
    return {k: sorted(v) for k, v in starts.items()}


def build_mp_index(alns: list[MpAlignment]) -> MpIndex:
    idx: MpIndex = defaultdict(list)
    for a in alns:
        idx[(a.contig, a.strand)].append(a)
    return dict(idx)


def find_overlapping(
    orf: ORF,
    mp_index: MpIndex,
    min_overlap_frac: float,
) -> list[MpAlignment]:
    """Return miniprot alignments that overlap orf by ≥ min_overlap_frac of its CDS."""
    candidates = mp_index.get((orf.contig, orf.strand), [])
    if not candidates:
        return []
    orf_cds_len = orf.coding_len
    result = []
    for mp in candidates:
        # Fast span check
        if mp.cds_segments[-1][1] <= orf.span_start or \
           mp.cds_segments[0][0] >= orf.span_end:
            continue
        # CDS-level overlap
        overlap = 0
        for os, oe in orf.segments:
            for ms, me in mp.cds_segments:
                lo = max(os, ms)
                hi = min(oe, me)
                if hi > lo:
                    overlap += hi - lo
        if overlap >= min_overlap_frac * orf_cds_len:
            result.append(mp)
    return result


# ──────────────────────────────────────────────────────────────────────────────
# Stop codon scanning
# ──────────────────────────────────────────────────────────────────────────────

def scan_for_stop(
    contig: str,
    start: int,
    strand: str,
    genome: Fasta,
    contig_lens: dict[str, int],
    max_scan: int,
) -> int | None:
    """Return 0-based position of the first in-frame stop codon within max_scan nt of start.

    + strand: checks start, start+3, start+6, ...
    - strand: checks start, start-3, start-6, ...
              (start is the 0-based start of the highest-coord triplet to check first)

    max_scan=0 checks only the single triplet at start.
    Returns None if no stop is found within the window.
    """
    clen = contig_lens.get(contig, 0)
    if strand == "+":
        for step in range(0, max_scan + 3, 3):
            pos = start + step
            if pos + 3 > clen:
                break
            if str(genome[contig][pos:pos + 3]).upper() in STOP_CODONS:
                return pos
    else:
        for step in range(0, max_scan + 3, 3):
            pos = start - step
            if pos < 0 or pos + 3 > clen:
                break
            if _rev_comp(str(genome[contig][pos:pos + 3]).upper()) in STOP_CODONS:
                return pos
    return None


def scan_for_start(
    contig: str,
    start: int,
    strand: str,
    genome: Fasta,
    contig_lens: dict[str, int],
    max_scan: int,
) -> int | None:
    """Return 0-based position of the first in-frame ATG within max_scan nt of start.

    + strand: checks start, start-3, start-6, ... (upstream = decreasing pos)
              start should be mp.cds_segments[0][0] (5'-most in-frame boundary).
    - strand: checks start, start+3, start+6, ... (upstream = increasing pos)
              start should be mp.cds_segments[-1][1] - 3 (5'-most codon).

    Returns None if no ATG is found within the scan window.
    """
    clen = contig_lens.get(contig, 0)
    if strand == "+":
        for step in range(0, max_scan + 3, 3):
            pos = start - step
            if pos < 0:
                break
            if str(genome[contig][pos:pos + 3]).upper() == "ATG":
                return pos
    else:
        for step in range(0, max_scan + 3, 3):
            pos = start + step
            if pos + 3 > clen:
                break
            if _rev_comp(str(genome[contig][pos:pos + 3]).upper()) == "ATG":
                return pos
    return None


# ──────────────────────────────────────────────────────────────────────────────
# Extension logic
# ──────────────────────────────────────────────────────────────────────────────

def cds_sequence(
    segs: list[tuple[int, int]], contig: str, strand: str, genome: Fasta,
) -> str:
    """Coding sequence of 0-based half-open segments in reading direction."""
    seq = "".join(str(genome[contig][s:e]) for s, e in sorted(segs)).upper()
    return _rev_comp(seq) if strand == "-" else seq


def has_internal_stop(
    segs: list[tuple[int, int]], contig: str, strand: str, genome: Fasta,
) -> bool:
    """True if the CDS has a stop codon in frame before its last codon.

    The extension builders check the length modulo 3 and the terminal stop,
    but an extension through miniprot CDS blocks or hint introns can cross a
    stop codon in frame. Such a CDS would translate to a truncated protein
    and must not be written.
    """
    seq = cds_sequence(segs, contig, strand, genome)
    return any(seq[i:i + 3] in STOP_CODONS for i in range(0, len(seq) - 3, 3))


def _introns(segs: list[tuple[int, int]]) -> set[tuple[int, int]]:
    segs = sorted(segs)
    return {(a[1], b[0]) for a, b in zip(segs, segs[1:])}


def splice_sites_ok(segs: list[tuple[int, int]], orf: ORF, genome: Fasta) -> bool:
    """True if every intron of segs that the ORF does not have is canonical.

    The extensions join the ORF to miniprot CDS blocks or hint introns, so a
    new intron can start where the transcript happened to end, or be one of
    the 1-2 nt gaps miniprot leaves at a frameshift. Such an intron must be
    at least MIN_INTRON nt long and read GT...AG or GC...AG on the ORF strand.
    """
    global n_bad_splice
    for s, e in _introns(segs) - _introns(orf.segments):
        seq = str(genome[orf.contig][s:e]).upper()
        if orf.strand == "-":
            seq = _rev_comp(seq)
        if e - s < MIN_INTRON or seq[:2] not in DONORS or seq[-2:] != "AG":
            n_bad_splice += 1
            return False
    return True


def _upstream_ext(
    mp_segs: list[tuple[int, int]], strand: str, atg_pos: int, orf_first: int, orf_last: int,
) -> list[tuple[int, int]] | None:
    """CDS segments from the ATG at atg_pos to the ORF, along the miniprot exons.

    An ATG upstream of the alignment's 5' end is joined to its first exon
    without an intron; an ATG that is not inside a miniprot exon (it lies in
    one of its introns) gives None.
    """
    atg_end = atg_pos + 3
    if strand == "+":
        ext = [(max(s, atg_pos), min(e, orf_first))
               for s, e in mp_segs if s < orf_first and e > atg_pos]
        if atg_pos < mp_segs[0][0]:
            ext.append((atg_pos, mp_segs[0][0]))
    else:
        ext = [(max(s, orf_last), min(e, atg_end))
               for s, e in mp_segs if e > orf_last and s < atg_end]
        if atg_end > mp_segs[-1][1]:
            ext.append((mp_segs[-1][1], atg_end))
    ext = _merge_adjacent([(s, e) for s, e in ext if s < e])
    if not any(s <= atg_pos and atg_end <= e for s, e in ext):
        return None
    return ext


def _merge_adjacent(segs: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not segs:
        return []
    segs = sorted(segs)
    merged = [list(segs[0])]
    for s, e in segs[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def _build_extended(
    orf_segs: list[tuple[int, int]],
    mp_segs: list[tuple[int, int]],
    is_complete: bool,
    strand: str,
    max_extension: int,
    stop_pos: int,
) -> list[tuple[int, int]] | None:
    """Return merged CDS segments for the extended ORF, or None on failure.

    stop_pos is the 0-based genomic start of the verified stop codon (3 nt).

    For + strand:
      ext_boundary = last CDS end (partial) or last CDS end − 3 (complete).
                  If it lies in a miniprot intron whose donor is in the ORF's
                  last exon, the ORF read into that intron (a retained intron
                  or a transcript that ends in it): ext_boundary moves back to
                  the donor, so the CDS splices with the protein.
      ORF base: orf_segs clipped to ext_boundary.
      Extension:  mp_segs starting at or past ext_boundary.
      Gap fill:   codons between mp_last_end and stop_pos (when stop is further
                  out than the miniprot CDS end — in-frame by scan_for_stop).
      Stop codon: (stop_pos, stop_pos + 3).

    For − strand (mirrored):
      ext_boundary = first CDS start (partial) or first CDS start + 3 (complete),
                  moved up to the miniprot donor in the same way.
      ORF base: orf_segs at or above ext_boundary.
      Extension:  mp_segs ending at or before ext_boundary.
      Gap fill:   codons between stop_pos+3 and mp_first_start.
      Stop codon: (stop_pos, stop_pos + 3).

    Returns None if the extension exceeds max_extension or if the resulting
    total coding length is not divisible by 3.
    """
    if strand == "+":
        orf_last = max(e for _, e in orf_segs)
        ext_boundary = orf_last - 3 if is_complete else orf_last
        mp_last_end = max(e for _, e in mp_segs)
        stop_end = stop_pos + 3

        if stop_end - orf_last > max_extension:
            return None

        last_exon_start = max(orf_segs)[0]
        for (_, donor), (acceptor, _) in zip(mp_segs, mp_segs[1:]):
            if donor < ext_boundary < acceptor and donor > last_exon_start:
                ext_boundary = donor
                break

        base = [(s, min(e, ext_boundary)) for s, e in orf_segs if s < ext_boundary]
        base = [(s, e) for s, e in base if s < e]
        ext = [(max(s, ext_boundary), e) for s, e in mp_segs if e > ext_boundary]
        ext = [(s, e) for s, e in ext if s < e]
        # Codons between miniprot CDS end and the scanned stop (in-frame by construction)
        gap = [(mp_last_end, stop_pos)] if stop_pos > mp_last_end else []
        stop_seg = (stop_pos, stop_end)

    else:  # "-"
        orf_first = min(s for s, _ in orf_segs)
        ext_boundary = orf_first + 3 if is_complete else orf_first
        mp_first_start = min(s for s, _ in mp_segs)
        stop_end = stop_pos + 3  # stop occupies [stop_pos, stop_pos+3]

        if orf_first - stop_pos > max_extension:
            return None

        first_exon_end = min(orf_segs)[1]
        for (_, acceptor), (donor, _) in zip(mp_segs, mp_segs[1:]):
            if acceptor < ext_boundary < donor and donor < first_exon_end:
                ext_boundary = donor
                break

        base = [(max(s, ext_boundary), e) for s, e in orf_segs if e > ext_boundary]
        base = [(s, e) for s, e in base if s < e]
        ext = [(s, min(e, ext_boundary)) for s, e in mp_segs if s < ext_boundary]
        ext = [(s, e) for s, e in ext if s < e]
        # Codons between the scanned stop end and miniprot first CDS start
        gap = [(stop_end, mp_first_start)] if stop_end < mp_first_start else []
        stop_seg = (stop_pos, stop_end)

    new_segs = _merge_adjacent(base + ext + gap + [stop_seg])
    if not new_segs:
        return None
    total = sum(e - s for s, e in new_segs)
    if total % 3 != 0:
        return None
    return new_segs


# ──────────────────────────────────────────────────────────────────────────────
# Hint-based extension (miniprothint intron hints)
# ──────────────────────────────────────────────────────────────────────────────

def _build_from_hints(
    orf_segs: list[tuple[int, int]],
    contig: str,
    strand: str,
    trunc_pos: int,
    sorted_hints: list[tuple[int, int]],
    genome: Fasta,
    contig_lens: dict[str, int],
    max_extension: int,
    max_stop_scan: int,
) -> list[tuple[int, int]] | None:
    """Build a corrected CDS using miniprothint intron positions.

    + strand: sorted_hints is ascending by intron_start; trunc_pos = orf_end.
    - strand: sorted_hints is descending by intron_end; trunc_pos = orf_first.

    Extension exons are placed between consecutive introns; a stop codon is
    located by in-frame scanning after the last intron.  Returns None if no
    stop is found, the extension exceeds max_extension, or the total coding
    length is not divisible by 3.
    """
    if strand == "+":
        pos = trunc_pos
        ext_exons: list[tuple[int, int]] = []

        for hint_start, hint_end in sorted_hints:
            if hint_start <= pos:
                continue
            if hint_start - trunc_pos > max_extension:
                break
            ext_exons.append((pos, hint_start))   # exon before this intron
            pos = hint_end                          # jump past intron

        sp = scan_for_stop(contig, pos, "+", genome, contig_lens, max_stop_scan)
        if sp is None or sp + 3 - trunc_pos > max_extension:
            return None
        ext_exons.append((pos, sp + 3))            # last exon including stop

        base = [(s, min(e, trunc_pos)) for s, e in orf_segs if s < trunc_pos]
        base = [(s, e) for s, e in base if s < e]
        new_segs = _merge_adjacent(base + ext_exons)

    else:  # "-"
        pos = trunc_pos
        ext_exons = []

        for hint_start, hint_end in sorted_hints:  # descending by hint_end
            if hint_end >= pos:
                continue
            if trunc_pos - hint_end > max_extension:
                break
            ext_exons.append((hint_end, pos))      # exon before this intron
            pos = hint_start                        # jump past intron (going left)

        sp = scan_for_stop(contig, pos - 3, "-", genome, contig_lens, max_stop_scan)
        if sp is None or trunc_pos - sp > max_extension:
            return None
        ext_exons.append((sp, pos))                # last exon including stop

        base = [(max(s, trunc_pos), e) for s, e in orf_segs if e > trunc_pos]
        base = [(s, e) for s, e in base if s < e]
        new_segs = _merge_adjacent(ext_exons + base)

    if not new_segs:
        return None
    total = sum(e - s for s, e in new_segs)
    if total % 3 != 0:
        return None
    return new_segs


def _try_fix_with_hints(
    orf: ORF,
    mp_alns: list[MpAlignment],
    hint_index: HintIndex,
    genome: Fasta,
    contig_lens: dict[str, int],
    max_extension: int,
    max_stop_scan: int,
) -> list[tuple[int, int]] | None:
    """Attempt stop recovery via miniprothint intron hints (partial ORFs only).

    Requires at least one miniprot alignment that extends past the truncation
    to confirm protein coverage in the extension region.  Hints are filtered
    to the footprint of the best-extending alignment.  Returns None if no
    relevant hints exist or the stop cannot be located.
    """
    contig, strand = orf.contig, orf.strand
    all_hints = hint_index.get((contig, strand), [])
    if not all_hints:
        return None

    if strand == "+":
        trunc_pos = max(e for _, e in orf.segments)
        ext_ends = [m.cds_segments[-1][1] for m in mp_alns
                    if m.cds_segments[-1][1] > trunc_pos]
        if not ext_ends:
            return None
        search_limit = max(ext_ends) + max_stop_scan + 100
        hints = [(s, e) for s, e in all_hints
                 if s >= trunc_pos and e <= search_limit]
        if not hints:
            return None
        return _build_from_hints(orf.segments, contig, strand, trunc_pos,
                                  sorted(hints),
                                  genome, contig_lens, max_extension, max_stop_scan)

    else:
        trunc_pos = min(s for s, _ in orf.segments)
        ext_starts = [m.cds_segments[0][0] for m in mp_alns
                      if m.cds_segments[0][0] < trunc_pos]
        if not ext_starts:
            return None
        search_limit = min(ext_starts) - max_stop_scan - 100
        hints = sorted(
            [(s, e) for s, e in all_hints
             if e <= trunc_pos and s >= search_limit],
            key=lambda h: h[1], reverse=True,   # closest intron_end first
        )
        if not hints:
            return None
        return _build_from_hints(orf.segments, contig, strand, trunc_pos,
                                  hints,
                                  genome, contig_lens, max_extension, max_stop_scan)


# ──────────────────────────────────────────────────────────────────────────────
# Per-ORF fix attempt
# ──────────────────────────────────────────────────────────────────────────────

def try_fix(
    orf: ORF,
    mp_alns: list[MpAlignment],
    genome: Fasta,
    contig_lens: dict[str, int],
    max_extension: int,
    max_stop_scan: int,
    hint_index: HintIndex | None = None,
) -> tuple[list[tuple[int, int]] | None, bool]:
    """Try to locate a protein-supported stop codon for the ORF.

    Returns (new_segs, via_hints) where via_hints is True when miniprothint
    intron hints were used to build the extension (INTRON case).

    Strategy:
      1. If hint_index is provided and orf is partial: try hint-based extension
         first.  This uses miniprothint intron positions to build correct splice
         sites across any intron(s) between the truncation point and the stop.
      2. Fall back to raw miniprot CDS: scan in-frame from the end of the
         furthest-reaching miniprot CDS block.
    """
    orf_end = max(e for _, e in orf.segments)
    orf_first = min(s for s, _ in orf.segments)

    # ── (1) hint-based path (partial ORFs only) ───────────────────────────────
    if hint_index is not None and orf.is_partial:
        new_segs = _try_fix_with_hints(
            orf, mp_alns, hint_index, genome, contig_lens,
            max_extension, max_stop_scan,
        )
        if new_segs is not None and splice_sites_ok(new_segs, orf, genome):
            return new_segs, True

    # ── (2) raw miniprot CDS path ─────────────────────────────────────────────
    if orf.strand == "+":
        ordered = sorted(mp_alns, key=lambda m: m.cds_segments[-1][1], reverse=True)
    else:
        ordered = sorted(mp_alns, key=lambda m: m.cds_segments[0][0])

    for mp in ordered:
        if orf.strand == "+":
            mp_cds_end = mp.cds_segments[-1][1]
            if mp_cds_end <= orf_end:
                continue
            sp = scan_for_stop(
                orf.contig, mp_cds_end, "+", genome, contig_lens, max_stop_scan,
            )
        else:
            mp_cds_start = mp.cds_segments[0][0]
            if mp_cds_start >= orf_first:
                continue
            sp = scan_for_stop(
                orf.contig, mp_cds_start - 3, "-", genome, contig_lens, max_stop_scan,
            )

        if sp is None:
            continue

        new_segs = _build_extended(
            orf.segments, mp.cds_segments,
            not orf.is_partial, orf.strand, max_extension, sp,
        )
        if new_segs is not None and splice_sites_ok(new_segs, orf, genome):
            return new_segs, False

    return None, False


def try_fix_5prime(
    orf: ORF,
    mp_alns: list[MpAlignment],
    genome: Fasta,
    contig_lens: dict[str, int],
    max_extension: int,
    max_start_scan: int,
) -> list[tuple[int, int]] | None:
    """Extend a 5'-partial ORF upstream to a protein-supported ATG.

    Finds miniprot alignments that begin upstream of the ORF's first CDS base
    (+ strand) or end downstream of the ORF's last CDS base (- strand), then
    scans for an ATG at or near the miniprot 5' boundary.  Uses the miniprot
    CDS exon structure to handle any introns in the extension region.

    Returns merged CDS segments (starting with ATG) or None.
    """
    orf_first = min(s for s, _ in orf.segments)
    orf_last = max(e for _, e in orf.segments)

    if orf.strand == "+":
        ordered = sorted(
            [m for m in mp_alns if m.cds_segments[0][0] < orf_first],
            key=lambda m: m.cds_segments[0][0],
        )
        for mp in ordered:
            mp_cds_start = mp.cds_segments[0][0]
            if orf_first - mp_cds_start > max_extension:
                continue
            atg_pos = scan_for_start(
                orf.contig, mp_cds_start, "+", genome, contig_lens, max_start_scan,
            )
            if atg_pos is None or atg_pos >= orf_first:
                continue
            # Build extension using miniprot exon structure from atg_pos to orf_first
            ext = _upstream_ext(mp.cds_segments, "+", atg_pos, orf_first, orf_last)
            if ext is None:
                continue
            new_segs = _merge_adjacent(ext + list(orf.segments))
            if sum(e - s for s, e in new_segs) % 3 != 0:
                continue
            if not splice_sites_ok(new_segs, orf, genome):
                continue
            return new_segs

    else:  # "-"
        # On - strand, transcript 5' = high genomic coord; ATG is at high coord
        ordered = sorted(
            [m for m in mp_alns if m.cds_segments[-1][1] > orf_last],
            key=lambda m: m.cds_segments[-1][1], reverse=True,
        )
        for mp in ordered:
            mp_cds_end = mp.cds_segments[-1][1]
            if mp_cds_end - orf_last > max_extension:
                continue
            # Scan rightward (= upstream on - strand) from the mp 5'-most codon
            atg_pos = scan_for_start(
                orf.contig, mp_cds_end - 3, "-", genome, contig_lens, max_start_scan,
            )
            if atg_pos is None or atg_pos + 3 <= orf_last:
                continue
            ext = _upstream_ext(mp.cds_segments, "-", atg_pos, orf_first, orf_last)
            if ext is None:
                continue
            new_segs = _merge_adjacent(list(orf.segments) + ext)
            if sum(e - s for s, e in new_segs) % 3 != 0:
                continue
            if not splice_sites_ok(new_segs, orf, genome):
                continue
            return new_segs

    return None


def try_fix_start_with_hints(
    orf: ORF,
    overlapping: list[MpAlignment],
    start_index: StartIndex,
    genome: Fasta,
    contig_lens: dict[str, int],
    max_extension: int,
    max_start_scan: int,
) -> list[tuple[int, int]] | None:
    """Extend an ORF 5' using a miniprothint start_codon hint anchored to an
    overlapping miniprot alignment.

    For each overlapping alignment whose 5' boundary reaches upstream of the
    ORF, look for a start_codon hint within max_start_scan nt of that boundary.
    If found, use the hint position as the new ATG and bridge to orf_first via
    the alignment's own CDS exon structure (same splice logic as try_fix_5prime).

    This mirrors GeneMark-ETP's ProtHint usage: a start hint is only applied
    when it belongs to a protein that already overlaps the predicted gene body.
    Hint score is used to rank candidates; ties broken by proximity to the
    alignment's 5' boundary.
    """
    contig, strand = orf.contig, orf.strand
    all_starts = start_index.get((contig, strand), [])
    if not all_starts:
        return None

    orf_first = min(s for s, _ in orf.segments)
    orf_last = max(e for _, e in orf.segments)
    if strand == "+":
        candidates_mp = sorted(
            [m for m in overlapping if m.cds_segments[0][0] < orf_first],
            key=lambda m: m.cds_segments[0][0],
        )
        for mp in candidates_mp:
            mp_cds_start = mp.cds_segments[0][0]
            if orf_first - mp_cds_start > max_extension:
                continue
            hint_candidates = [
                (p, sc) for p, sc in all_starts
                if p < orf_first and abs(p - mp_cds_start) <= max_start_scan
            ]
            if not hint_candidates:
                continue
            # Best score first; ties: prefer closest to mp 5' boundary
            hint_candidates.sort(key=lambda x: (-x[1], abs(x[0] - mp_cds_start)))
            for atg_pos, _ in hint_candidates:
                ext = _upstream_ext(mp.cds_segments, "+", atg_pos, orf_first, orf_last)
                if ext is None:
                    continue
                new_segs = _merge_adjacent(ext + list(orf.segments))
                if sum(e - s for s, e in new_segs) % 3 != 0:
                    continue
                if not splice_sites_ok(new_segs, orf, genome):
                    continue
                return new_segs

    else:  # "-"
        candidates_mp = sorted(
            [m for m in overlapping if m.cds_segments[-1][1] > orf_last],
            key=lambda m: m.cds_segments[-1][1], reverse=True,
        )
        for mp in candidates_mp:
            mp_cds_end = mp.cds_segments[-1][1]
            if mp_cds_end - orf_last > max_extension:
                continue
            mp_5prime_codon = mp_cds_end - 3  # 5'-most codon on − strand
            hint_candidates = [
                (p, sc) for p, sc in all_starts
                if p + 3 > orf_last and abs(p - mp_5prime_codon) <= max_start_scan
            ]
            if not hint_candidates:
                continue
            hint_candidates.sort(key=lambda x: (-x[1], abs(x[0] - mp_5prime_codon)))
            for atg_pos, _ in hint_candidates:
                ext = _upstream_ext(mp.cds_segments, "-", atg_pos, orf_first, orf_last)
                if ext is None:
                    continue
                new_segs = _merge_adjacent(list(orf.segments) + ext)
                if sum(e - s for s, e in new_segs) % 3 != 0:
                    continue
                if not splice_sites_ok(new_segs, orf, genome):
                    continue
                return new_segs

    return None


# ──────────────────────────────────────────────────────────────────────────────
# GTF writing
# ──────────────────────────────────────────────────────────────────────────────

def _gtf_lines(
    tid: str,
    segments: list[tuple[int, int]],
    contig: str,
    strand: str,
    source: str,
) -> list[str]:
    """Emit one CDS GTF line per segment with correct GTF phase."""
    reading_order = sorted(segments) if strand == "+" else sorted(segments, reverse=True)
    phase_map: dict[tuple[int, int], int] = {}
    coding_so_far = 0
    for seg in reading_order:
        phase_map[seg] = (3 - coding_so_far % 3) % 3
        coding_so_far += seg[1] - seg[0]

    lines = []
    for s, e in sorted(segments):
        ph = phase_map[(s, e)]
        lines.append(
            f"{contig}\t{source}\tCDS\t{s + 1}\t{e}\t.\t{strand}\t{ph}"
            f'\ttranscript_id "{tid}"; gene_id "{tid}";'
        )
    return lines


def _passthrough_lines(orf: ORF) -> list[str]:
    return _gtf_lines(orf.tid, orf.segments, orf.contig, orf.strand, orf.source)


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--orfs", type=Path, required=True,
                    help="Complete predicted ORF GTF (from annotate.py).")
    ap.add_argument("--partial", type=Path, default=None,
                    help="Partial ORF GTF (from annotate.py --partial-out). Optional.")
    ap.add_argument("--partial5", type=Path, default=None,
                    help="5'-partial ORF GTF (from annotate.py --partial5-out). Optional.")
    ap.add_argument("--miniprot", type=Path, required=True,
                    help="Miniprot GFF3 output (protein-to-genome alignments).")
    ap.add_argument("--genome", type=Path, required=True,
                    help="Genome FASTA (pyfaidx-readable).")
    ap.add_argument("--out", type=Path, required=True,
                    help="Output corrected GTF.")
    ap.add_argument("--max-extension", type=int, default=5000,
                    help="Max bp a stop codon may be extended downstream of "
                         "the current ORF end (default 5000).")
    ap.add_argument("--min-overlap-frac", type=float, default=0.3,
                    help="Min fraction of ORF CDS that must overlap the "
                         "miniprot alignment CDS to consider it a match "
                         "(default 0.3).")
    ap.add_argument("--max-stop-scan", type=int, default=30,
                    help="Max nt to scan in-frame beyond the miniprot CDS end "
                         "when searching for a stop codon (default 30, i.e. "
                         "10 codons). Set to 0 to check only the triplet "
                         "immediately after the last aligned CDS block.")
    ap.add_argument("--hints", type=Path, default=None,
                    help="miniprothint hc.gff (intron hints from boundary-scored "
                         "miniprot alignments).  When provided, partial ORFs are "
                         "first extended using hint-derived intron positions before "
                         "falling back to the raw miniprot CDS approach.")
    ap.add_argument("--fix-starts", action="store_true",
                    help="Also attempt to extend complete ORFs upstream to a "
                         "protein-supported ATG using try_fix_5prime(). "
                         "Controlled by --fix-starts-classes.")
    ap.add_argument("--fix-starts-classes", nargs="+",
                    default=["LORF_NOUPSTOP", "upLORF"],
                    metavar="CLASS",
                    help="LORF classes eligible for start-fixing (default: "
                         "LORF_NOUPSTOP upLORF). Ignored unless --fix-starts "
                         "is set. Pass 'ALL' to attempt on every complete ORF.")
    ap.add_argument("--max-start-scan", type=int, default=30,
                    help="Max nt to scan in-frame upstream from the miniprot "
                         "5' boundary when searching for an ATG (default 30). "
                         "Used with --fix-starts and for --partial5 ORFs.")
    return ap.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)

    print(f"Loading complete ORFs:   {args.orfs}", file=sys.stderr)
    complete_orfs = parse_orf_gtf(args.orfs, is_partial=False)
    print(f"  {len(complete_orfs)} transcripts", file=sys.stderr)

    partial_orfs: dict[str, ORF] = {}
    if args.partial is not None:
        print(f"Loading partial ORFs:    {args.partial}", file=sys.stderr)
        partial_orfs = parse_orf_gtf(args.partial, is_partial=True)
        print(f"  {len(partial_orfs)} transcripts", file=sys.stderr)

    partial5_orfs: dict[str, ORF] = {}
    if args.partial5 is not None:
        print(f"Loading 5'-partial ORFs: {args.partial5}", file=sys.stderr)
        partial5_orfs = parse_orf_gtf(args.partial5, is_partial=True)
        print(f"  {len(partial5_orfs)} transcripts", file=sys.stderr)

    print(f"Loading miniprot GFF3:   {args.miniprot}", file=sys.stderr)
    mp_alns = parse_miniprot_gff(args.miniprot)
    n_with_stop = sum(1 for a in mp_alns if a.has_stop)
    print(f"  {len(mp_alns)} alignments total, {n_with_stop} with StopCodon=1",
          file=sys.stderr)
    mp_index = build_mp_index(mp_alns)

    hint_index: HintIndex | None = None
    start_index: StartIndex | None = None
    if args.hints is not None:
        print(f"Loading miniprothint hints: {args.hints}", file=sys.stderr)
        hint_index = parse_miniprothint_hints(args.hints)
        n_hint_introns = sum(len(v) for v in hint_index.values())
        print(f"  {n_hint_introns} intron hints on {len(hint_index)} contig/strand pairs",
              file=sys.stderr)
        if args.fix_starts:
            start_index = parse_miniprothint_starts(args.hints)
            n_start_hints = sum(len(v) for v in start_index.values())
            print(f"  {n_start_hints} start_codon hints for start-fixing",
                  file=sys.stderr)

    print(f"Loading genome FASTA:    {args.genome}", file=sys.stderr)
    genome = Fasta(str(args.genome), as_raw=True, sequence_always_upper=True)
    contig_lens = {k: len(genome[k]) for k in genome.keys()}

    all_orfs: list[ORF] = list(complete_orfs.values()) + list(partial_orfs.values())

    n_complete_fixed = 0
    n_complete_unchanged = 0
    n_complete_start_fixed = 0
    n_complete_start_hint = 0   # of start-fixed, how many used miniprothint hints
    n_partial_recovered = 0
    n_partial_hint = 0      # recovered via miniprothint intron hints
    n_partial_dropped = 0
    n_partial5_recovered = 0
    n_partial5_dropped = 0
    n_internal_stop = 0     # corrections dropped for an in-frame stop codon

    def clean(segs, orf):
        """segs, or None if the corrected CDS has an in-frame stop codon."""
        nonlocal n_internal_stop
        if segs is not None and has_internal_stop(segs, orf.contig, orf.strand, genome):
            n_internal_stop += 1
            return None
        return segs

    fix_starts_classes: set[str] = set()
    if args.fix_starts:
        if "ALL" in args.fix_starts_classes:
            fix_starts_classes = {"ALL"}
        else:
            fix_starts_classes = set(args.fix_starts_classes)

    with open(args.out, "w") as fh:
        for orf in sorted(all_orfs, key=lambda o: o.tid):
            overlapping = find_overlapping(orf, mp_index, args.min_overlap_frac)
            new_segs, via_hints = try_fix(
                orf, overlapping, genome, contig_lens,
                args.max_extension, args.max_stop_scan,
                hint_index=hint_index,
            )
            new_segs = clean(new_segs, orf)

            if new_segs is not None:
                lines = _gtf_lines(
                    orf.tid, new_segs, orf.contig, orf.strand, orf.source,
                )
                if orf.is_partial:
                    n_partial_recovered += 1
                    if via_hints:
                        n_partial_hint += 1
                else:
                    n_complete_fixed += 1
            elif not orf.is_partial:
                # Attempt upstream start-fix for eligible complete ORFs.
                if fix_starts_classes and (
                    "ALL" in fix_starts_classes
                    or orf.lorf_class in fix_starts_classes
                ):
                    start_segs = None
                    via_hint_start = False
                    if start_index is not None:
                        # miniprothint hints available: only fix when a hint is
                        # anchored to an overlapping alignment (GeneMark-ETP style).
                        # No fallback to raw miniprot — unanchored ATG picks are
                        # too noisy and consistently hurt accuracy.
                        start_segs = try_fix_start_with_hints(
                            orf, overlapping, start_index, genome, contig_lens,
                            args.max_extension, args.max_start_scan,
                        )
                        if start_segs is not None:
                            via_hint_start = True
                    else:
                        # No hints provided: use raw miniprot 5' boundary scan.
                        start_segs = try_fix_5prime(
                            orf, overlapping, genome, contig_lens,
                            args.max_extension, args.max_start_scan,
                        )
                    start_segs = clean(start_segs, orf)
                    if start_segs is not None:
                        lines = _gtf_lines(
                            orf.tid, start_segs, orf.contig, orf.strand, orf.source,
                        )
                        n_complete_start_fixed += 1
                        if via_hint_start:
                            n_complete_start_hint += 1
                        for line in lines:
                            fh.write(line + "\n")
                        continue
                lines = _passthrough_lines(orf)
                n_complete_unchanged += 1
            else:
                n_partial_dropped += 1
                continue

            for line in lines:
                fh.write(line + "\n")

        # Process 5'-partial ORFs with try_fix_5prime().
        for orf in sorted(partial5_orfs.values(), key=lambda o: o.tid):
            overlapping = find_overlapping(orf, mp_index, args.min_overlap_frac)
            new_segs = clean(try_fix_5prime(
                orf, overlapping, genome, contig_lens,
                args.max_extension, args.max_start_scan,
            ), orf)
            if new_segs is not None:
                lines = _gtf_lines(
                    orf.tid, new_segs, orf.contig, orf.strand, orf.source,
                )
                n_partial5_recovered += 1
                for line in lines:
                    fh.write(line + "\n")
            else:
                n_partial5_dropped += 1

    hint_note = f" (of which {n_partial_hint} via miniprothint intron hints)" \
                if hint_index is not None else ""
    partial5_note = (
        f"\n5'-partial ORFs: {n_partial5_recovered} recovered, "
        f"{n_partial5_dropped} dropped (no protein support)."
        if partial5_orfs else ""
    )
    start_note = (
        f"\nStart-fixed    : {n_complete_start_fixed} complete ORFs extended upstream "
        f"({n_complete_start_hint} via miniprothint hints, "
        f"{n_complete_start_fixed - n_complete_start_hint} via raw miniprot fallback)."
        if fix_starts_classes else ""
    )
    print(
        f"Complete ORFs : {n_complete_fixed} stop-fixed, "
        f"{n_complete_unchanged} unchanged.{start_note}\n"
        f"Partial ORFs  : {n_partial_recovered} recovered{hint_note}, "
        f"{n_partial_dropped} dropped (no protein support).{partial5_note}\n"
        f"Corrections dropped for an in-frame stop codon: {n_internal_stop}.\n"
        f"Candidate corrections rejected for a non-canonical new intron: {n_bad_splice}.",
        file=sys.stderr,
    )
    print(f"Output: {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
