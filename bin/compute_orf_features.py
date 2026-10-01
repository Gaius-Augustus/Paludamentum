#!/usr/bin/env python3
# From tiberius_orf_finder/scripts/compute_orf_features.py (Lars Gabriel), used unchanged in the Drusilla flow.
# Copyright (c) 2026 Lars Gabriel. Artistic License 1.0, see LICENSE.
"""Compute per-ORF feature table for scoring and analysis.

The table holds the protein, LORF and hint features of every transcript; it
is the input of apply_lgb_model_gtf.py (the LightGBM filter of the Drusilla
flow).

Inputs
------
--orfs-gtf       CDS-only GTF from annotate.py (with lorf_class attribute)
--miniprot-gff   Miniprot protein-to-genome GFF3 (miniprot_scored.gff)
--proteins-fasta Protein FASTA used for miniprot alignment (enables best_protein_coverage)
--hints-gff      miniprothint hc.gff (intron/start_codon/stop_codon hints)
--genome         Genome FASTA (pyfaidx-readable)
--ref-tmap       Optional gffcompare .tmap output for match labels
--upstream-scan  Max nt to scan upstream for stop codon (default 1500)
--out            Output TSV path

Output columns (one row per transcript)
----------------------------------------
transcript_id contig strand n_exons cds_length_nt lorf_class
dist_upstream_stop_nt n_upstream_atgs
has_protein_support n_overlapping_alignments
best_identity best_score best_norm_bitscore best_target_coverage best_protein_coverage
protein_extends_5prime_codons protein_extends_3prime_codons
n_introns_supported frac_introns_supported
has_start_hint has_stop_hint support_level
has_conflict conflict_identity_delta
cds_length_pct has_upstream_partner has_downstream_partner
[gffcompare_class if --ref-tmap provided]

Notes on approximations
-----------------------
dist_upstream_stop_nt: scanned in GENOMIC coordinates in coding frame.
For genes with intron-containing 5' UTRs this underestimates the true
transcript-space distance; lorf_class (computed in transcript space by
annotate.py) is the authoritative categorical version.

partialness_score (GeneMark-ETP): full formula needs qstart_norm
(protein N-terminal coverage fraction) which requires the protein database.
protein_extends_5prime_codons covers the targetStartDiff component.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import math
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from pyfaidx import Fasta

_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")
STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})


def _rev_comp(seq: str) -> str:
    return seq.translate(_RC)[::-1]


def _attr(col: str, key: str) -> str | None:
    m = re.search(rf'{key}\s+"([^"]+)"', col)
    if m:
        return m.group(1)
    m = re.search(rf'{key}=([^;"\s]+)', col)
    if m:
        return m.group(1)
    return None


def _parse_target(col: str) -> tuple[str, int, int]:
    """Parse GFF3 Target=proteinID start end attribute.  Returns ("", 0, 0) when absent."""
    m = re.search(r'Target=([^;\s]+)\s+(\d+)\s+(\d+)', col)
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3))
    return "", 0, 0


# ── Data structures ────────────────────────────────────────────────────────────

@dataclass
class ORFRec:
    tid: str
    contig: str
    strand: str
    segments: list[tuple[int, int]] = field(default_factory=list)
    lorf_class: str | None = None

    @property
    def n_exons(self) -> int:
        return len(self.segments)

    @property
    def cds_length(self) -> int:
        return sum(e - s for s, e in self.segments)

    @property
    def orf_first(self) -> int:
        return self.segments[0][0]

    @property
    def orf_last(self) -> int:
        return self.segments[-1][1]

    @property
    def atg_pos(self) -> int:
        return self.orf_first if self.strand == "+" else self.orf_last - 3

    @property
    def stop_pos(self) -> int:
        return self.orf_last - 3 if self.strand == "+" else self.orf_first

    @property
    def introns(self) -> list[tuple[int, int]]:
        return [(self.segments[i][1], self.segments[i + 1][0])
                for i in range(len(self.segments) - 1)]


@dataclass
class MpAln:
    mid: str
    contig: str
    strand: str
    segments: list[tuple[int, int]] = field(default_factory=list)
    identity: float = 0.0
    score: int = 0
    target_id: str = ""      # protein sequence ID from Target= attribute
    target_aln_aa: int = 0   # aligned amino acids on the target protein side

    @property
    def aln_first(self) -> int:
        return self.segments[0][0]

    @property
    def aln_last(self) -> int:
        return self.segments[-1][1]

    @property
    def aligned_aa(self) -> int:
        return max(1, sum(e - s for s, e in self.segments) // 3)


# ── Parsing ────────────────────────────────────────────────────────────────────

def parse_orfs(path: Path) -> list[ORFRec]:
    orfs: dict[str, ORFRec] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            tid = _attr(f[8], "transcript_id")
            if not tid:
                continue
            if tid not in orfs:
                orfs[tid] = ORFRec(tid=tid, contig=f[0], strand=f[6])
            orfs[tid].segments.append((int(f[3]) - 1, int(f[4])))
            if orfs[tid].lorf_class is None:
                lc = _attr(f[8], "lorf_class")
                if lc:
                    orfs[tid].lorf_class = lc
    for o in orfs.values():
        o.segments.sort()
    return list(orfs.values())


def parse_miniprot(path: Path) -> list[MpAln]:
    pending: dict[str, MpAln] = {}
    for row in csv.reader(open(path), delimiter="\t"):
        if not row or row[0].startswith("#") or len(row) < 9:
            continue
        if row[2] == "mRNA":
            mid = _attr(row[8], "ID")
            if not mid:
                continue
            ident = float(_attr(row[8], "Identity") or "0")
            # Score is in GFF column 6 (0-based index 5), not an attribute
            sc = int(row[5]) if row[5] not in (".", "", "*") else 0
            tid, tstart, tend = _parse_target(row[8])
            target_aln_aa = max(0, tend - tstart)
            pending[mid] = MpAln(mid=mid, contig=row[0], strand=row[6],
                                 identity=ident, score=sc,
                                 target_id=tid, target_aln_aa=target_aln_aa)
        elif row[2] == "CDS":
            parent = _attr(row[8], "Parent")
            if parent in pending:
                pending[parent].segments.append((int(row[3]) - 1, int(row[4])))
    result = [a for a in pending.values() if a.segments]
    for a in result:
        a.segments.sort()
    return result


def parse_hints(path: Path) -> tuple[
    dict[tuple[str, str], list[tuple[int, int]]],
    dict[tuple[str, str], list[int]],
    dict[tuple[str, str], list[int]],
]:
    """Return (intron_index, start_codon_index, stop_codon_index)."""
    introns: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    starts: dict[tuple[str, str], list[int]] = defaultdict(list)
    stops: dict[tuple[str, str], list[int]] = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            if len(f) < 9:
                continue
            key = (f[0], f[6])
            s0 = int(f[3]) - 1  # 0-based
            if f[2] == "intron":
                introns[key].append((s0, int(f[4])))
            elif f[2] == "start_codon":
                starts[key].append(s0)
            elif f[2] == "stop_codon":
                stops[key].append(s0)
    return (
        {k: sorted(v) for k, v in introns.items()},
        {k: sorted(v) for k, v in starts.items()},
        {k: sorted(v) for k, v in stops.items()},
    )


def parse_tmap(path: Path) -> dict[str, str]:
    """Return {qry_id: class_code} from a gffcompare .tmap file."""
    result: dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("ref_gene_id") or not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) >= 4:
                result[cols[4].strip()] = cols[2].strip()
    return result


# ── Overlap index ──────────────────────────────────────────────────────────────

MpIndex = dict[tuple[str, str], list[MpAln]]


def build_mp_index(alns: list[MpAln]) -> MpIndex:
    idx: MpIndex = defaultdict(list)
    for a in alns:
        idx[(a.contig, a.strand)].append(a)
    return dict(idx)


def overlapping_mp(orf: ORFRec, mp_index: MpIndex,
                   min_frac: float) -> list[MpAln]:
    candidates = mp_index.get((orf.contig, orf.strand), [])
    orf_len = orf.cds_length
    result = []
    for mp in candidates:
        if mp.aln_last <= orf.orf_first or mp.aln_first >= orf.orf_last:
            continue
        overlap = sum(
            max(0, min(oe, me) - max(os, ms))
            for os, oe in orf.segments
            for ms, me in mp.segments
        )
        if overlap >= min_frac * orf_len:
            result.append(mp)
    return result


# ── Feature computation ────────────────────────────────────────────────────────

def upstream_features(orf: ORFRec, genome: Fasta,
                      contig_lens: dict[str, int],
                      max_scan: int) -> dict:
    """Scan genomic sequence upstream of ATG in coding frame (approximation —
    ignores possible introns in 5' UTR; use lorf_class for transcript-accurate
    categorical classification)."""
    contig = orf.contig
    clen = contig_lens.get(contig, 0)
    n_atgs = 0
    dist = None

    if orf.strand == "+":
        atg = orf.orf_first
        for step in range(3, max_scan + 3, 3):
            pos = atg - step
            if pos < 0:
                break
            codon = str(genome[contig][pos:pos + 3]).upper()
            if codon in STOP_CODONS:
                dist = step
                break
            if codon == "ATG":
                n_atgs += 1
    else:
        atg = orf.orf_last - 3
        for step in range(3, max_scan + 3, 3):
            pos = atg + step
            if pos + 3 > clen:
                break
            codon = _rev_comp(str(genome[contig][pos:pos + 3]).upper())
            if codon in STOP_CODONS:
                dist = step
                break
            if codon == "ATG":
                n_atgs += 1

    return {"dist_upstream_stop_nt": dist, "n_upstream_atgs": n_atgs}


def protein_features(orf: ORFRec, alns: list[MpAln],
                     protein_lens: dict[str, int]) -> dict:
    nan = float("nan")
    if not alns:
        return dict(
            has_protein_support=False,
            n_overlapping_alignments=0,
            best_identity=nan,
            best_score=nan,
            best_norm_bitscore=nan,
            best_target_coverage=nan,
            best_protein_coverage=nan,
            protein_extends_5prime_codons=0,
            protein_extends_3prime_codons=0,
            has_conflict=False,
            conflict_identity_delta=nan,
        )

    orf_len = orf.cds_length

    def _tcov(mp: MpAln) -> float:
        ov = sum(max(0, min(oe, me) - max(os, ms))
                 for os, oe in orf.segments for ms, me in mp.segments)
        return ov / orf_len if orf_len else 0.0

    ranked = sorted(alns, key=lambda m: m.identity * _tcov(m), reverse=True)
    best = ranked[0]
    tc = _tcov(best)

    if orf.strand == "+":
        ext5 = max(0, orf.orf_first - best.aln_first) // 3
        ext3 = max(0, best.aln_last  - orf.orf_last)  // 3
    else:
        ext5 = max(0, best.aln_last  - orf.orf_last)  // 3
        ext3 = max(0, orf.orf_first  - best.aln_first) // 3

    # Conflict: another alignment suggests a different 5' boundary by ≥ 5 codons
    has_conflict = False
    conflict_delta = nan
    for mp in ranked[1:]:
        if orf.strand == "+":
            diff = abs(mp.aln_first - best.aln_first) // 3
        else:
            diff = abs(mp.aln_last - best.aln_last) // 3
        if diff >= 5:
            has_conflict = True
            conflict_delta = mp.identity - best.identity
            break

    prot_len = protein_lens.get(best.target_id, 0)
    best_protein_coverage = best.target_aln_aa / prot_len if prot_len > 0 else nan

    return dict(
        has_protein_support=True,
        n_overlapping_alignments=len(alns),
        best_identity=best.identity,
        best_score=float(best.score),
        best_norm_bitscore=best.score / best.aligned_aa,
        best_target_coverage=tc,
        best_protein_coverage=best_protein_coverage,
        protein_extends_5prime_codons=ext5,
        protein_extends_3prime_codons=ext3,
        has_conflict=has_conflict,
        conflict_identity_delta=conflict_delta,
    )


def hint_features(orf: ORFRec,
                  introns: dict[tuple[str, str], list[tuple[int, int]]],
                  starts: dict[tuple[str, str], list[int]],
                  stops: dict[tuple[str, str], list[int]]) -> dict:
    key = (orf.contig, orf.strand)
    hint_intron_set = set(introns.get(key, []))
    hint_starts     = starts.get(key, [])
    hint_stops      = stops.get(key, [])

    orf_introns = orf.introns
    n_sup = sum(1 for iv in orf_introns if iv in hint_intron_set)
    frac_sup = n_sup / len(orf_introns) if orf_introns else float("nan")

    # Within ±3 nt to tolerate minor GFF coordinate conventions
    has_start = any(abs(p - orf.atg_pos)  <= 3 for p in hint_starts)
    has_stop  = any(abs(p - orf.stop_pos) <= 3 for p in hint_stops)

    if (not orf_introns and has_start and has_stop) or \
       (orf_introns and n_sup == len(orf_introns) and has_start and has_stop):
        level = "fullSupport"
    elif n_sup > 0 or has_start or has_stop:
        level = "anySupport"
    else:
        level = "noSupport"

    return dict(
        n_introns_supported=n_sup,
        frac_introns_supported=frac_sup,
        has_start_hint=has_start,
        has_stop_hint=has_stop,
        support_level=level,
    )


def split_gene_flags(
    orfs: list[ORFRec],
    max_gap: int,
) -> tuple[set[str], set[str]]:
    """Flag ORFs that have a plausible split-gene partner on the same strand.

    An upstream partner: another ORF on the same (contig, strand) whose 3' end
    is within max_gap of our 5' end and does not overlap us.
    A downstream partner: symmetric.
    The frame compatibility check (gap % 3 == 0) filters out incompatible pairs.
    """
    groups: dict[tuple[str, str], list[ORFRec]] = defaultdict(list)
    for orf in orfs:
        groups[(orf.contig, orf.strand)].append(orf)
    for g in groups.values():
        g.sort(key=lambda o: o.orf_first)

    has_up: set[str] = set()
    has_dn: set[str] = set()

    for group in groups.values():
        for i, orf in enumerate(group):
            for j in range(i - 1, -1, -1):
                prev = group[j]
                gap = orf.orf_first - prev.orf_last
                if gap < 0:
                    break  # overlapping
                if gap > max_gap:
                    break
                if gap % 3 != 0:
                    continue  # frame-incompatible
                has_up.add(orf.tid)
                has_dn.add(prev.tid)
                break

    return has_up, has_dn


# ── Output ────────────────────────────────────────────────────────────────────

COLUMNS = [
    "transcript_id", "contig", "strand", "n_exons", "cds_length_nt", "lorf_class",
    "dist_upstream_stop_nt", "n_upstream_atgs",
    "has_protein_support", "n_overlapping_alignments",
    "best_identity", "best_score", "best_norm_bitscore", "best_target_coverage",
    "best_protein_coverage",
    "protein_extends_5prime_codons", "protein_extends_3prime_codons",
    "n_introns_supported", "frac_introns_supported",
    "has_start_hint", "has_stop_hint", "support_level",
    "has_conflict", "conflict_identity_delta",
    "cds_length_pct", "n_overlapping_alignments_pct",
    "has_upstream_partner", "has_downstream_partner",
]


def _fmt(v) -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "NA"
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, float):
        return f"{v:.4f}"
    return str(v)


# ── CLI ───────────────────────────────────────────────────────────────────────

def _parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--orfs-gtf",       required=True, type=Path)
    ap.add_argument("--miniprot-gff",   required=False, default=None, type=Path,
                    help="Miniprot scored GFF (optional; protein features zeroed when absent)")
    ap.add_argument("--proteins-fasta", required=False, default=None, type=Path,
                    help="Protein FASTA used for miniprot alignment; enables best_protein_coverage")
    ap.add_argument("--hints-gff",      required=False, default=None, type=Path,
                    help="miniprothint hc.gff (optional; hint features zeroed when absent)")
    ap.add_argument("--genome",         required=True, type=Path)
    ap.add_argument("--out",            required=True, type=Path)
    ap.add_argument("--ref-tmap",       default=None,  type=Path,
                    help="gffcompare .tmap for match labels (optional)")
    ap.add_argument("--upstream-scan",  type=int, default=1500,
                    help="Max nt upstream genomic scan for stop codon (default 1500)")
    ap.add_argument("--min-overlap-frac", type=float, default=0.3)
    ap.add_argument("--split-gene-max-gap", type=int, default=5000,
                    help="Max intergenic gap to flag as split-gene partner (default 5000)")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = _parse_args(argv)
    cols = list(COLUMNS)
    if args.ref_tmap:
        cols.append("gffcompare_class")

    print(f"Parsing ORF GTF: {args.orfs_gtf}", file=sys.stderr)
    orfs = parse_orfs(args.orfs_gtf)
    print(f"  {len(orfs)} ORFs", file=sys.stderr)

    protein_lens: dict[str, int] = {}
    if args.proteins_fasta is not None:
        print(f"Loading protein lengths: {args.proteins_fasta}", file=sys.stderr)
        prot_fa = Fasta(str(args.proteins_fasta), as_raw=True)
        protein_lens = {k: len(prot_fa[k]) for k in prot_fa.keys()}
        print(f"  {len(protein_lens)} proteins", file=sys.stderr)

    if args.miniprot_gff is not None:
        print(f"Parsing miniprot GFF: {args.miniprot_gff}", file=sys.stderr)
        mp_alns = parse_miniprot(args.miniprot_gff)
        mp_index = build_mp_index(mp_alns)
        print(f"  {len(mp_alns)} alignments", file=sys.stderr)
        if not protein_lens:
            print("  Warning: --proteins-fasta not provided; best_protein_coverage will be NaN",
                  file=sys.stderr)
    else:
        print("No miniprot GFF provided — protein features will be zero", file=sys.stderr)
        mp_alns, mp_index = [], {}

    if args.hints_gff is not None:
        print(f"Parsing miniprothint hints: {args.hints_gff}", file=sys.stderr)
        hint_introns, hint_starts, hint_stops = parse_hints(args.hints_gff)
        n_hi = sum(len(v) for v in hint_introns.values())
        print(f"  {n_hi} intron hints", file=sys.stderr)
    else:
        print("No hints GFF provided — hint features will be zero", file=sys.stderr)
        hint_introns, hint_starts, hint_stops = {}, {}, {}

    print(f"Loading genome: {args.genome}", file=sys.stderr)
    genome = Fasta(str(args.genome), as_raw=True, sequence_always_upper=True)
    contig_lens = {k: len(genome[k]) for k in genome.keys()}

    tmap: dict[str, str] = {}
    if args.ref_tmap:
        print(f"Parsing tmap: {args.ref_tmap}", file=sys.stderr)
        tmap = parse_tmap(args.ref_tmap)

    print("Computing split-gene partner flags ...", file=sys.stderr)
    has_up, has_dn = split_gene_flags(orfs, args.split_gene_max_gap)

    # Pre-compute alignments for all ORFs (needed for percentile ranking).
    print("Pre-computing protein alignments ...", file=sys.stderr)
    orf_alns = [(orf, overlapping_mp(orf, mp_index, args.min_overlap_frac)) for orf in orfs]

    # CDS length percentile across all ORFs
    lengths = sorted(o.cds_length for o in orfs)
    n = len(lengths)

    def _pct(length: int) -> float:
        rank = bisect.bisect_left(lengths, length)
        return rank / n if n else float("nan")

    # n_overlapping_alignments percentile across all ORFs
    aln_counts_sorted = sorted(len(alns) for _, alns in orf_alns)
    n_aln = len(aln_counts_sorted)

    def _aln_pct(count: int) -> float:
        rank = bisect.bisect_left(aln_counts_sorted, count)
        return rank / n_aln if n_aln else float("nan")

    print(f"Writing features to {args.out} ...", file=sys.stderr)
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for orf, alns in orf_alns:
            pf = protein_features(orf, alns, protein_lens)
            row: dict = {
                "transcript_id": orf.tid,
                "contig": orf.contig,
                "strand": orf.strand,
                "n_exons": orf.n_exons,
                "cds_length_nt": orf.cds_length,
                "lorf_class": orf.lorf_class or "NA",
                **upstream_features(orf, genome, contig_lens, args.upstream_scan),
                **pf,
                **hint_features(orf, hint_introns, hint_starts, hint_stops),
                "cds_length_pct": _pct(orf.cds_length),
                "n_overlapping_alignments_pct": _aln_pct(pf["n_overlapping_alignments"]),
                "has_upstream_partner": orf.tid in has_up,
                "has_downstream_partner": orf.tid in has_dn,
            }
            if args.ref_tmap:
                row["gffcompare_class"] = tmap.get(orf.tid, "NA")
            w.writerow([_fmt(row[c]) for c in cols])

    print(f"Done -> {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
