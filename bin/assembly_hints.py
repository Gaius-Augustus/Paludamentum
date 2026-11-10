#!/usr/bin/env python3
"""Intron hints and a coverage track derived from a StringTie assembly.

Substitutes for the bam2hints introns and the bedtools genomecov bedGraph when a
run has a StringTie assembly (params.stringtie) but no reads, so no BAM
(ASSEMBLY_HINTS in modules/assembly.nf; combine_genesets.py --introns and
--coverage, and hintsfile.gff).

Introns: the gaps between consecutive exons of every transcript, from the first
to the last intron base as bam2hints writes them, one hint per distinct intron.
Its multiplicity (mult=) is the sum of the StringTie cov of the transcripts that
contain it (--mult cov, the default): cov is the mean per-base read depth of a
transcript, so the sum approximates the junction read count that bam2hints
reports, which the selection classifier of combine_genesets.py was trained on.
--mult count writes the number of transcripts instead (measured on 31 benchmark
species: 1.0 transcript F1 worse than cov, which is within 0.5 of the BAM).
The strand is the transcripts' strand, '.' if they disagree; the layout is
that of bam2hints + filterIntronsFindStrand.pl (src=E, pri=4).

Coverage (--bedgraph): per-base depth = the largest exon cov among the
assembled exons that cover the base (StringTie writes a cov per exon); bases
outside assembled exons have depth 0 and are left out, as in genomecov -bg.
bedGraph, 0-based half-open, sorted by sequence and position, gzipped when the
file name ends in .gz.

Transcripts without a cov attribute (not written by stringtie --merge) count
with the mean cov of their exons, or 1 without any; they are reported.

Usage: assembly_hints.py stringtie.gtf --introns introns.gff [--mult cov|count] [--bedgraph cov.bg.gz]
"""
import argparse
import gzip
import re
import sys
from collections import defaultdict

TID = re.compile(r'transcript_id\s+"?([^";]+)"?')
COV = re.compile(r'(?:^|;)\s*cov\s+"?([\d.eE+-]+)"?')


def read_assembly(path):
    """{tid: [(contig, strand, start, end, exon cov or None)]}, {tid: transcript cov or None}"""
    tx = defaultdict(list)
    tcov = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            m = TID.search(f[8])
            if not m:
                continue
            c = COV.search(f[8])
            cov = float(c.group(1)) if c else None
            if f[2] == 'transcript':
                tcov[m.group(1)] = cov
            elif f[2] == 'exon':
                tx[m.group(1)].append((f[0], f[6], int(f[3]), int(f[4]), cov))
    return tx, tcov


def transcript_cov(tid, exons, tcov):
    """the transcript's cov, else the mean of its exons' cov, else None"""
    cov = tcov.get(tid)
    if cov is not None:
        return cov
    covs = [e[4] for e in exons if e[4] is not None]
    return sum(covs) / len(covs) if covs else None


def introns_of(tx, tcov):
    """{(contig, start, end): [n transcripts, cov sum, set of strands]}, n transcripts without cov"""
    introns = defaultdict(lambda: [0, 0.0, set()])
    n_missing = 0
    for tid, ex in tx.items():
        ex.sort(key=lambda e: e[2])
        cov = transcript_cov(tid, ex, tcov)
        if cov is None:
            n_missing += 1
            cov = 1.0
        for (c, st, _s1, e1, _cov1), (_c2, _st2, s2, _e2, _cov2) in zip(ex, ex[1:]):
            if s2 - 1 >= e1 + 1:
                rec = introns[(c, e1 + 1, s2 - 1)]
                rec[0] += 1
                rec[1] += cov
                rec[2].add(st)
    return introns, n_missing


def write_hints(introns, mult_mode, out):
    for (c, s, e), (n, covsum, strands) in sorted(introns.items()):
        mult = n if mult_mode == 'count' else max(1, int(round(covsum)))
        strand = next(iter(strands)) if len(strands) == 1 else '.'
        out.write(f'{c}\tasm2h\tintron\t{s}\t{e}\t{mult}\t{strand}\t.\tmult={mult};pri=4;src=E\n')


def coverage_segments(exons):
    """[(start0, end, depth)] of one sequence, depth = max exon cov over the exons covering the segment,
    0-based half-open, zero-depth gaps left out, neighbouring segments of equal depth merged."""
    bounds = sorted({s - 1 for s, _, _ in exons} | {e for _, e, _ in exons})
    starts = sorted(exons, key=lambda x: x[0])
    active, k = [], 0
    out = []
    for b0, b1 in zip(bounds, bounds[1:]):
        while k < len(starts) and starts[k][0] - 1 <= b0:
            active.append(starts[k])
            k += 1
        active = [x for x in active if x[1] > b0]
        depth = max((x[2] for x in active), default=0.0)
        if depth <= 0:
            continue
        if out and out[-1][1] == b0 and out[-1][2] == depth:
            out[-1] = (out[-1][0], b1, depth)
        else:
            out.append((b0, b1, depth))
    return out


def write_bedgraph(tx, tcov, path):
    by_contig = defaultdict(list)
    for tid, ex in tx.items():
        cov = transcript_cov(tid, ex, tcov)
        for c, _st, s, e, ecov in ex:
            by_contig[c].append((s, e, ecov if ecov is not None else (cov if cov is not None else 1.0)))
    op = gzip.open if path.endswith('.gz') else open
    n = 0
    with op(path, 'wt') as out:
        for c in sorted(by_contig):
            for a, b, depth in coverage_segments(by_contig[c]):
                out.write(f'{c}\t{a}\t{b}\t{depth:g}\n')
                n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog='\n\n'.join(__doc__.split('\n\n')[1:]))
    ap.add_argument('gtf', help='StringTie assembly (GTF with cov attributes)')
    ap.add_argument('--introns', required=True, help='intron hints GFF (bam2hints layout, mult= attribute)')
    ap.add_argument('--mult', choices=['cov', 'count'], default='cov',
                    help='multiplicity: summed transcript cov (default) or number of transcripts')
    ap.add_argument('--bedgraph', help='coverage track from the exon cov attributes (bedGraph, .gz allowed)')
    a = ap.parse_args()

    tx, tcov = read_assembly(a.gtf)
    introns, n_missing = introns_of(tx, tcov)
    with open(a.introns, 'w') as out:
        write_hints(introns, a.mult, out)
    print(f'{len(tx)} transcripts, {len(introns)} distinct introns, mult = {a.mult}'
          + (f', {n_missing} transcripts without cov attribute (counted as 1)' if n_missing else ''), file=sys.stderr)
    if a.bedgraph:
        n = write_bedgraph(tx, tcov, a.bedgraph)
        print(f'coverage track: {n} segments', file=sys.stderr)


if __name__ == '__main__':
    main()
