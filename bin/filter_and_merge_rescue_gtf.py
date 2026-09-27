#!/usr/bin/env python3
# From tiberius_orf_finder/scripts/filter_and_merge_rescue_gtf.py (Lars Gabriel), used unchanged in the Drusilla flow.
"""
Filter and deduplicate Tiberius hint-rescue predictions.

Filters applied to each predicted transcript:
  1. Strand must match the dominant strand of intron hints for that locus entry.
  2. The transcript's genomic span (in local coords) must overlap the region
     covered by the chain's hints — this excludes spurious flanking predictions.

Deduplication:
  Transcripts with identical (chr, strand, frozenset-of-CDS-intervals) are
  collapsed to one representative.  Gene/transcript IDs are renumbered.

Usage:
    python filter_and_merge_rescue_gtf.py \\
        raw_gtf hints_gff manifest_tsv out_gtf
"""
import argparse
import re
import sys
from collections import defaultdict


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('raw_gtf',      help='Tiberius output on locus_id sequences (local coords)')
    p.add_argument('hints_gff',    help='combined_hints.gff (local coords, chain-tagged)')
    p.add_argument('manifest_tsv', help='loci_manifest.tsv (locus_id→chr,bed_start,...)')
    p.add_argument('out_gtf',      help='filtered, deduplicated genome-coordinate GTF')
    return p.parse_args()


# ── Helpers ───────────────────────────────────────────────────────────────────

def gtf_attr(attr_col, key):
    m = re.search(r'%s "([^"]+)"' % re.escape(key), attr_col)
    return m.group(1) if m else None


# ── Load hints (local coordinates) ───────────────────────────────────────────

def load_hints(hints_gff):
    """
    Per locus_id: dominant strand of intron hints + total span of all hints.
    Strand is decided by majority vote among intron features; all hints
    (intron + codon) contribute to the span.
    """
    raw = defaultdict(lambda: {'strand_votes': defaultdict(int), 'spans': []})
    with open(hints_gff) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) < 8:
                continue
            locus_id = parts[0]
            feature  = parts[2].lower()
            strand   = parts[6]
            try:
                start = int(parts[3])
                end   = int(parts[4])
            except ValueError:
                continue
            if feature == 'intron':
                raw[locus_id]['strand_votes'][strand] += 1
            raw[locus_id]['spans'].append((start, end))

    result = {}
    for lid, d in raw.items():
        votes = d['strand_votes']
        if not votes:
            # No intron hints — fall back to first hint's strand from spans
            # (we won't know strand; skip strand filtering for this locus)
            result[lid] = {'strand': None,
                           'span_start': min(s[0] for s in d['spans']),
                           'span_end':   max(s[1] for s in d['spans'])}
            continue
        dominant = max(votes, key=votes.get)
        result[lid] = {
            'strand':     dominant,
            'span_start': min(s[0] for s in d['spans']),
            'span_end':   max(s[1] for s in d['spans']),
        }
    return result


# ── Load manifest ─────────────────────────────────────────────────────────────

def load_manifest(manifest_tsv):
    # columns: locus_id  chr  bed_start  bed_end  strand  chain_id
    offset = {}   # locus_id -> (chr, bed_start, strand)
    with open(manifest_tsv) as fh:
        next(fh)  # header
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue
            locus_id  = parts[0]
            chrom     = parts[1]
            bed_start = int(parts[2])
            strand    = parts[4] if len(parts) > 4 else None
            offset[locus_id] = (chrom, bed_start, strand)
    return offset


# ── Parse raw GTF into per-transcript records ─────────────────────────────────

def load_transcripts(raw_gtf):
    """
    Group GTF lines by (locus_id, transcript_id).  Returns:
      tx_lines : uid -> [raw lines]
      tx_meta  : uid -> dict(locus_id, strand, local_min, local_max, cds)
    """
    tx_lines = defaultdict(list)
    tx_meta  = {}

    with open(raw_gtf) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            locus_id = parts[0]
            feature  = parts[2]
            strand   = parts[6]
            try:
                start = int(parts[3])
                end   = int(parts[4])
            except ValueError:
                continue
            tx_id = gtf_attr(parts[8], 'transcript_id')
            if not tx_id:
                continue

            # Prefix with locus_id to avoid ID collisions across loci
            uid = f'{locus_id}::{tx_id}'
            tx_lines[uid].append(line)

            if uid not in tx_meta:
                tx_meta[uid] = {
                    'locus_id':  locus_id,
                    'strand':    strand,
                    'local_min': start,
                    'local_max': end,
                    'cds':       [],
                }
            else:
                tx_meta[uid]['local_min'] = min(tx_meta[uid]['local_min'], start)
                tx_meta[uid]['local_max'] = max(tx_meta[uid]['local_max'], end)

            if feature == 'CDS':
                tx_meta[uid]['cds'].append((start, end))

    return tx_lines, tx_meta


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    locus_hints = load_hints(args.hints_gff)
    offset      = load_manifest(args.manifest_tsv)
    tx_lines, tx_meta = load_transcripts(args.raw_gtf)

    print(f'[filter] {len(tx_lines)} transcripts loaded from raw GTF', file=sys.stderr)

    n_no_locus     = 0     # locus_id not in manifest
    n_wrong_strand = 0
    n_no_overlap   = 0
    n_chainless    = 0     # accepted from ab-initio (chain-less) loci
    candidates     = []   # (fingerprint, sort_key, uid, genome_lines)

    for uid, lines in tx_lines.items():
        meta     = tx_meta[uid]
        locus_id = meta['locus_id']

        # Manifest lookup is mandatory (we need chrom + bed_start to project coords)
        if locus_id not in offset:
            n_no_locus += 1
            continue
        chrom, bed_start, manifest_strand = offset[locus_id]

        # Manifest-strand cross-check
        if manifest_strand is not None and meta['strand'] != manifest_strand:
            n_wrong_strand += 1
            continue

        # Chain-less loci (absent from combined_hints.gff) skip hint-based filters
        # and are accepted based on the manifest strand check alone.
        if locus_id in locus_hints:
            hint_info = locus_hints[locus_id]

            # 1. Strand check against hint dominant strand
            if hint_info['strand'] is not None and meta['strand'] != hint_info['strand']:
                n_wrong_strand += 1
                continue

            # 2. Transcript span must overlap the hint-covered region (local coords)
            if (meta['local_max'] < hint_info['span_start'] or
                    meta['local_min'] > hint_info['span_end']):
                n_no_overlap += 1
                continue
        else:
            n_chainless += 1

        genome_lines = []
        cds_genome   = []
        g_start_min  = float('inf')
        for line in lines:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 5:
                try:
                    gs = int(parts[3]) + bed_start
                    ge = int(parts[4]) + bed_start
                    parts[0] = chrom
                    parts[3] = str(gs)
                    parts[4] = str(ge)
                    if parts[2] == 'CDS':
                        cds_genome.append((gs, ge))
                    g_start_min = min(g_start_min, gs)
                    genome_lines.append('\t'.join(parts) + '\n')
                except ValueError:
                    genome_lines.append(line)
            else:
                genome_lines.append(line)

        fp       = (chrom, meta['strand'], frozenset(cds_genome))
        sort_key = (chrom, g_start_min)
        candidates.append((fp, sort_key, uid, genome_lines))

    print(f'[filter] dropped: {n_wrong_strand} wrong-strand | '
          f'{n_no_overlap} outside hint region | {n_no_locus} no manifest entry | '
          f'{n_chainless} accepted from chain-less loci',
          file=sys.stderr)

    # Deduplicate by CDS fingerprint (keep first occurrence of each structure)
    candidates.sort(key=lambda x: x[1])   # stable sort by (chr, start)
    seen   = set()
    unique = []
    for fp, _, uid, genome_lines in candidates:
        if fp not in seen:
            seen.add(fp)
            unique.append((uid, genome_lines))

    n_dupes = len(candidates) - len(unique)
    print(f'[dedup]  {n_dupes} duplicate structures removed → '
          f'{len(unique)} unique transcripts', file=sys.stderr)

    # Write with renumbered gene/transcript IDs
    with open(args.out_gtf, 'w') as fout:
        for i, (uid, genome_lines) in enumerate(unique, 1):
            gid = f'rescued_g{i}'
            tid = f'rescued_t{i}'
            for line in genome_lines:
                line = re.sub(r'gene_id "[^"]+"',       f'gene_id "{gid}"',  line)
                line = re.sub(r'transcript_id "[^"]+"', f'transcript_id "{tid}"', line)
                fout.write(line)

    print(f'[out]    wrote {len(unique)} transcripts to {args.out_gtf}', file=sys.stderr)


if __name__ == '__main__':
    main()
