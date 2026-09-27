#!/usr/bin/env python3
# From tiberius_orf_finder/scripts/prepare_hint_rescue_loci.py (Lars Gabriel), used unchanged in the Drusilla flow.
"""
Prepare a multi-FASTA + combined hints GFF for Tiberius hint-guided rescue.

Logic per locus:
  1. Only rescue partial transcripts that have NO overlapping tib_correct
     prediction on the SAME strand.
  2. Use the single best protein chain (highest sum of al_score) per merged
     locus — one FASTA entry, one hint set, one Tiberius prediction.
  3. Skip rescue if any existing ORF transcript already contains ALL intron
     positions from the top chain (ORF agrees with the protein → no rescue
     needed).

Loci are merged per-strand so evidence from opposite strands stays separate.
Start and stop codon hints from the chain are included alongside intron hints.

Usage:
    python prepare_hint_rescue_loci.py \\
        --partial_gtf  tiberius_lgb_partial.gtf \\
        --correct_gtf  tiberius_lgb_correct.gtf \\
        --chained_hints chained_hints.gff \\
        --orfs_gtf     orfs.gtf [orfs.partial.gtf ...] \\
        --genome       genome.fa \\
        --outdir       hint_rescue/ \\
        [--flank 25000]
"""
import argparse
import bisect
import subprocess
import sys
from collections import defaultdict


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--partial_gtf',   required=True, help='tiberius_lgb_partial.gtf')
    p.add_argument('--correct_gtf',   required=True, help='tiberius_lgb_correct.gtf')
    p.add_argument('--chained_hints', required=True,
                   help='chain_id-tagged hints (output of chainedHints.py)')
    p.add_argument('--orfs_gtf', required=True, nargs='+',
                   help='ORF GTF file(s) for intron-agreement check (orfs.gtf etc.)')
    p.add_argument('--genome',  required=True, help='genome.fa (must have .fai)')
    p.add_argument('--outdir',  required=True, help='output directory')
    p.add_argument('--flank', type=int, default=25000,
                   help='bp flanking each side of a locus (default 25000)')
    return p.parse_args()


# ── GFF/GTF attribute helpers ─────────────────────────────────────────────────

def gff_attr(col, key):
    """Return value of 'key=' from a GFF9 attribute column, or None."""
    for field in col.split(';'):
        field = field.strip()
        if field.startswith(key + '='):
            return field[len(key)+1:].strip()
    return None

def gtf_attr(col, key):
    """Return value of 'key "value"' from a GTF attribute column, or None."""
    import re
    m = re.search(r'%s "([^"]+)"' % re.escape(key), col)
    return m.group(1) if m else None


# ── FAI / chrom sizes ─────────────────────────────────────────────────────────

def load_chrom_sizes(fai_path):
    sizes = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.split('\t')
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes


# ── Transcript loaders ────────────────────────────────────────────────────────

def load_transcripts(gtf_path):
    """(chrom, start0, end0, strand) for each transcript feature."""
    txs = []
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9 or parts[2] != 'transcript':
                continue
            try:
                txs.append((parts[0], int(parts[3]) - 1, int(parts[4]), parts[6]))
            except ValueError:
                continue
    return txs


# ── Strand-aware interval index ───────────────────────────────────────────────

def build_strand_index(transcripts):
    """O(log n) overlap query index keyed by (chrom, strand)."""
    raw = defaultdict(list)
    for chrom, start, end, strand in transcripts:
        raw[(chrom, strand)].append((start, end))
    index = {}
    for key, ivs in raw.items():
        ivs.sort()
        starts  = [iv[0] for iv in ivs]
        max_end = []
        cur = 0
        for iv in ivs:
            cur = max(cur, iv[1])
            max_end.append(cur)
        index[key] = (starts, max_end)
    return index

def any_overlap(index, chrom, strand, start0, end0):
    key = (chrom, strand)
    if key not in index:
        return False
    starts, max_end = index[key]
    right = bisect.bisect_left(starts, end0)
    return right > 0 and max_end[right - 1] > start0


# ── ORF intron index ──────────────────────────────────────────────────────────

def load_orf_introns(gtf_paths):
    """
    Build a per-(chrom, strand) list of (tx_start, tx_end, frozenset_of_introns).
    Introns are derived from consecutive CDS features (1-based inclusive coords).
    Single-exon transcripts are excluded (no introns to compare).
    """
    cds_by_tx = defaultdict(list)   # (chrom, strand, tx_id) -> [(start, end)]
    for path in gtf_paths:
        with open(path) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 9 or parts[2] != 'CDS':
                    continue
                tx_id = gtf_attr(parts[8], 'transcript_id')
                if not tx_id:
                    continue
                try:
                    txs_key = (parts[0], parts[6], tx_id)
                    cds_by_tx[txs_key].append((int(parts[3]), int(parts[4])))
                except ValueError:
                    continue

    raw = defaultdict(list)
    for (chrom, strand, _tx_id), cds in cds_by_tx.items():
        if len(cds) < 2:
            continue
        cds = sorted(cds)
        tx_start = cds[0][0]
        tx_end   = cds[-1][1]
        introns  = frozenset(
            (cds[i][1] + 1, cds[i+1][0] - 1)
            for i in range(len(cds) - 1)
        )
        raw[(chrom, strand)].append((tx_start, tx_end, introns))

    # Sort each chromosome list and build a max-end prefix for fast lookup
    orf_index = {}
    for key, txs in raw.items():
        txs.sort()
        starts  = [tx[0] for tx in txs]
        max_end = []
        cur = 0
        for tx in txs:
            cur = max(cur, tx[1])
            max_end.append(cur)
        orf_index[key] = (starts, max_end, txs)
    return orf_index

def orf_agrees(orf_index, chrom, strand, chain_introns, locus_start0, locus_end0):
    """
    True if any ORF transcript overlapping the locus contains all chain introns.
    chain_introns : frozenset of (start, end) tuples in 1-based genome coords.
    """
    if not chain_introns:
        return False
    key = (chrom, strand)
    if key not in orf_index:
        return False
    starts, max_end, txs = orf_index[key]
    right = bisect.bisect_left(starts, locus_end0)
    if right == 0 or max_end[right - 1] <= locus_start0:
        return False
    for i in range(right):
        tx_start, tx_end, orf_introns = txs[i]
        if tx_end <= locus_start0:
            continue
        if chain_introns.issubset(orf_introns):
            return True
    return False


# ── Chained hint index ────────────────────────────────────────────────────────

def load_chained_hints(gff_path):
    """
    Parse chain_id-tagged hints (genome coordinates).

    Returns
    -------
    intron_index   : chrom → (sorted_starts, max_end_prefix, hint_tuples)
                     hint_tuples: (start0, end0, chain_id, strand)
    hints_by_chain : (chrom, chain_id) → [raw_gff_lines]
    chain_scores   : (chrom, chain_id) → sum of al_score
    """
    introns_by_chrom = defaultdict(list)
    hints_by_chain   = defaultdict(list)
    chain_scores     = defaultdict(float)

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom   = parts[0]
            feature = parts[2].lower()
            strand  = parts[6] if len(parts) > 6 else '.'
            chain_id = gff_attr(parts[8], 'chain_id')
            if not chain_id:
                continue
            try:
                start0 = int(parts[3]) - 1   # 0-based
                end0   = int(parts[4])         # 0-based exclusive
            except ValueError:
                continue
            if end0 <= start0:
                continue

            al_score_str = gff_attr(parts[8], 'al_score')
            try:
                al_score = float(al_score_str) if al_score_str else 0.0
            except ValueError:
                al_score = 0.0

            hints_by_chain[(chrom, chain_id)].append(line.rstrip('\n'))
            chain_scores[(chrom, chain_id)] += al_score

            if feature == 'intron':
                introns_by_chrom[chrom].append((start0, end0, chain_id, strand))

    intron_index = {}
    for chrom, ivs in introns_by_chrom.items():
        ivs.sort()
        starts  = [iv[0] for iv in ivs]
        max_end = []
        cur = 0
        for iv in ivs:
            cur = max(cur, iv[1])
            max_end.append(cur)
        intron_index[chrom] = (starts, max_end, ivs)

    return intron_index, hints_by_chain, chain_scores


def best_chain_at_locus(intron_index, chain_scores, chrom, strand, start0, end0):
    """Return (chain_id, score) of the highest-scoring chain with intron hints at locus."""
    if chrom not in intron_index:
        return None, 0.0
    starts, max_end, ivs = intron_index[chrom]
    right = bisect.bisect_left(starts, end0)
    if right == 0 or max_end[right - 1] <= start0:
        return None, 0.0

    local_scores = defaultdict(float)
    for i in range(right):
        iv_start0, iv_end0, chain_id, iv_strand = ivs[i]
        if iv_end0 > start0 and iv_strand == strand:
            local_scores[chain_id] = chain_scores.get((chrom, chain_id), 0.0)

    if not local_scores:
        return None, 0.0
    best = max(local_scores, key=local_scores.get)
    return best, local_scores[best]


def get_chain_introns_genome(hints_lines):
    """frozenset of (start, end) 1-based intron positions from chain's hint lines."""
    introns = set()
    for line in hints_lines:
        parts = line.split('\t')
        if len(parts) < 9 or parts[2].lower() != 'intron':
            continue
        try:
            introns.add((int(parts[3]), int(parts[4])))
        except ValueError:
            pass
    return frozenset(introns)


# ── Interval merge (per strand) ───────────────────────────────────────────────

def merge_by_strand(ivs):
    """Merge overlapping (chrom, strand, start, end) tuples per strand."""
    if not ivs:
        return []
    ivs = sorted(ivs)
    merged = [list(ivs[0])]
    for chrom, strand, start, end in ivs[1:]:
        last = merged[-1]
        if chrom == last[0] and strand == last[1] and start <= last[3]:
            last[3] = max(last[3], end)
        else:
            merged.append([chrom, strand, start, end])
    return merged


# ── FASTA extraction ──────────────────────────────────────────────────────────

def fetch_sequence(genome, chrom, start0, end0):
    """Return bare sequence string for [start0, end0) via samtools faidx."""
    result = subprocess.run(
        ['samtools', 'faidx', genome, f'{chrom}:{start0 + 1}-{end0}'],
        capture_output=True, text=True, check=True,
    )
    lines = result.stdout.split('\n')
    return ''.join(lines[1:]).replace('\n', '')


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    import os
    os.makedirs(args.outdir, exist_ok=True)

    chrom_sizes = load_chrom_sizes(args.genome + '.fai')
    print(f'[index]  {len(chrom_sizes)} chromosomes', file=sys.stderr)

    print('[load]   chained hints ...', file=sys.stderr)
    intron_index, hints_by_chain, chain_scores = load_chained_hints(args.chained_hints)
    print(f'[load]   {len(hints_by_chain)} (chrom, chain) pairs', file=sys.stderr)

    correct_txs = load_transcripts(args.correct_gtf)
    correct_idx = build_strand_index(correct_txs)
    print(f'[load]   {len(correct_txs)} tib_correct transcripts', file=sys.stderr)

    orf_index = load_orf_introns(args.orfs_gtf)
    n_orf = sum(len(v[2]) for v in orf_index.values())
    print(f'[load]   {n_orf} multi-exon ORF transcripts', file=sys.stderr)

    partial_txs = load_transcripts(args.partial_gtf)
    print(f'[load]   {len(partial_txs)} tib_partial transcripts', file=sys.stderr)

    # ── Filter and collect eligible loci ─────────────────────────────────────
    # Keep every partial transcript that is not already covered by a
    # same-strand tib_correct prediction.  Chain-less loci are emitted as
    # hint-free FASTA entries so Tiberius still runs ab-initio on them
    # (option B).  ORF-agrees filter is dropped (option A).
    eligible  = []
    n_correct = 0

    for chrom, start0, end0, strand in partial_txs:
        if any_overlap(correct_idx, chrom, strand, start0, end0):
            n_correct += 1
            continue
        s = max(0, start0 - args.flank)
        e = min(chrom_sizes.get(chrom, end0 + args.flank), end0 + args.flank)
        eligible.append((chrom, strand, s, e))

    print(f'[filter] {n_correct} skipped: tib_correct overlaps same strand', file=sys.stderr)
    print(f'[filter] {len(eligible)} eligible partial loci', file=sys.stderr)

    merged = merge_by_strand(eligible)
    print(f'[merge]  {len(merged)} merged loci (flank={args.flank})', file=sys.stderr)

    # ── Build output ──────────────────────────────────────────────────────────
    out_fa   = os.path.join(args.outdir, 'combined_loci.fa')
    out_gff  = os.path.join(args.outdir, 'combined_hints.gff')
    out_mfst = os.path.join(args.outdir, 'loci_manifest.tsv')

    n_entries       = 0
    n_with_chain    = 0
    n_without_chain = 0

    with open(out_fa,   'w') as fa_fh, \
         open(out_gff,  'w') as gff_fh, \
         open(out_mfst, 'w') as mfst_fh:

        mfst_fh.write('locus_id\tchr\tbed_start\tbed_end\tstrand\tchain_id\n')

        for locus_idx, (chrom, strand, bed_start, bed_end) in enumerate(merged):
            best_chain, _ = best_chain_at_locus(
                intron_index, chain_scores, chrom, strand, bed_start, bed_end)

            chain_hints = hints_by_chain.get((chrom, best_chain), []) if best_chain else []

            seq = fetch_sequence(args.genome, chrom, bed_start, bed_end)
            if not seq:
                continue

            locus_id = f'locus_{n_entries:07d}'
            n_entries += 1
            if best_chain:
                n_with_chain += 1
            else:
                n_without_chain += 1

            # FASTA
            fa_fh.write(f'>{locus_id}\n')
            for i in range(0, len(seq), 60):
                fa_fh.write(seq[i:i+60] + '\n')

            # Hints (local 1-based coords; include intron + start/stop codons).
            # Chain-less loci write nothing here → Tiberius runs ab-initio on them.
            for raw_line in chain_hints:
                parts = raw_line.split('\t')
                try:
                    g_start = int(parts[3])
                    g_end   = int(parts[4])
                except (ValueError, IndexError):
                    continue
                if g_start > bed_end or g_end <= bed_start:
                    continue
                locus_len   = bed_end - bed_start
                local_start = max(1, g_start - bed_start)
                local_end   = min(locus_len, g_end - bed_start)
                if local_end <= 0 or local_start > locus_len:
                    continue
                parts[0] = locus_id
                parts[3] = str(local_start)
                parts[4] = str(local_end)
                gff_fh.write('\t'.join(parts) + '\n')

            chain_id_str = best_chain if best_chain else 'none'
            mfst_fh.write(
                f'{locus_id}\t{chrom}\t{bed_start}\t{bed_end}\t{strand}\t{chain_id_str}\n'
            )

            if (locus_idx + 1) % 500 == 0:
                print(f'[progress] {locus_idx+1}/{len(merged)} loci, '
                      f'{n_entries} entries so far', file=sys.stderr)

    print(f'[out]    {n_entries} rescue loci written '
          f'({n_with_chain} with chain hints, {n_without_chain} ab-initio)',
          file=sys.stderr)
    print(f'[out]    FASTA   : {out_fa}',   file=sys.stderr)
    print(f'[out]    hints   : {out_gff}',  file=sys.stderr)
    print(f'[out]    manifest: {out_mfst}', file=sys.stderr)


if __name__ == '__main__':
    main()
