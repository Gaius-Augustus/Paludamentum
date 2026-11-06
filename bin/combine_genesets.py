#!/usr/bin/env python3
r"""
combine_genesets.py - combine the Tiberius, Vipsania, Drusilla and TD2 gene sets into one annotation.

 1. MERGE    transcripts with an identical CDS (contig, strand, CDS segments incl. stop codon) across the
             four inputs become one structure; the set of tools that produced it is recorded.
 2. CLUSTER  structures with overlapping CDS spans on the same strand form a locus.
 3. TIER     from the tool set (built-in table below, or a clade table via --tier-table = 10 x precision):
               10 Tiberius+Vipsania+(Drusilla|TD2)   9 Tiberius+(Drusilla|TD2)   8 Tiberius+Vipsania
                7 Vipsania+Drusilla+TD2               5 Tiberius alone, plausible (>=5 exons | >=3 exons and
                4 Vipsania+(Drusilla|TD2)               >=600 nt | >=1200 nt)       3 Tiberius alone, other
                0 single transcript tool, Vipsania alone, Drusilla+TD2 without an ab initio tool
 4. ADJUST   lone ab initio models (Tiberius alone; with a tier table also Vipsania alone):
               up to 6    all introns RNA-supported / single exon >= --cov-full covered / all introns also
                          predicted by another tool / CDS chain inside one full-length transcript (--chains)
               up to 4.5  CDS identical to a miniprot model (--rescue-lone-miniprot)
               down to 2  best DIAMOND hit covers only part of a protein (length ratio outside 0.8-1.25),
                          unless identical to a miniprot model or all introns have protein hints
               down to 2  introns only partly RNA-supported (contradicting junctions); the model still marks
                          the locus, but another candidate must become primary (--no-drop-lone-rna-conflict)
               down to 2  single exon, for the tools in --drop-lone-single-exon (insects: tiberius)
 5. DECIDE   per locus: keep it if its best model reaches --min-locus (default 4); primary = best-ranked
             model (tier; Drusilla+TD2 agreement ranks 5.5; -2 for an RNA-unsupported intron; +1 for a CDS
             chain inside one full-length transcript); further isoforms if tier >= --iso-min-tier (8), or
             their CDS chain lies in one full-length transcript and they reach --iso-chain-min-tier (4), or
             Drusilla+TD2 agree on them, they share the primary's stop or >= 95% of its length, and their new
             introns are RNA-supported (--iso-dt filtered; --iso-dt-chains: and they lie in one full-length
             transcript).

Full-length transcripts (--chains) are transcripts assembled from long reads (StringTie -L). A model's CDS
chain lies in one of them ('full') if its introns are a run of the transcript's introns and its terminal CDS
segments lie in the flanking exons (a single-exon CDS: inside one exon).
All rules only act when their evidence is given (--introns, --coverage, --chains, DIAMOND tables, miniprot files).
--keep-all keeps every structure and switches all rules off (labelling / calibration runs).

Learned selection (--classifier MODEL.json, the Paludamentum default) replaces steps 3-5 by a probability:
the candidate pool is built as with --keep-all (every structure, all tier rules off), a gradient-boosted tree
model predicts P(structure = reference CDS) per candidate from the report features (tools, structure, locus
context, RNA and protein evidence; no clade input), and per locus the best candidate becomes primary if its
P reaches t_locus, further candidates are isoforms if their P reaches t_iso. The thresholds come from the
model file (--unknown-clade: its more conservative locus threshold for clades outside the training clades;
--t-locus / --t-iso override). --drop-lone-single-exon still applies: those structures are never selected.
Full-length transcripts (--chains) were not available for the training species, so they enter as a correction
of P: the log-odds of a structure whose CDS chain lies in one full-length transcript move by
--chain-logit-full, those of a 'partial' structure by --chain-logit-partial (calibrated on H. sapiens).
--min-locus, --iso-* and the RNA / protein rule switches have no effect in this mode; --tier-table is refused.
The model was trained on the keep-all pools of 15 species with gffcompare labels (see conf/ in Paludamentum).

Usage (recommended):
  combine_genesets.py --tiberius tib.gtf --vipsania vip.gtf --drusilla dru.gtf --td2 td2.gtf \
      --introns introns.gff --coverage cov.bg.gz \
      --diamond-tiberius t.tsv --diamond-vipsania v.tsv --diamond-drusilla d.tsv --diamond-td2 x.tsv \
      --miniprot-train miniprot_trainingGenes.gff --miniprot-rep miniprot_representatives.gff \
      --protein-hints hc.gff [--classifier combiner_classifier_v1.json [--unknown-clade]] \
      [--tier-table clade.tsv --min-locus 3] [--drop-lone-single-exon tiberius] \
      -o combined.gtf --report combined.report.tsv
"""
import argparse
import bisect
import json
import math
import re
import sys
from collections import defaultdict

TOOLS = ['tiberius', 'vipsania', 'drusilla', 'td2']

# ----------------------------------------------------------------------------- parsing

GTF_TID = re.compile(r'transcript_id\s+"?([^";]+)"?')


def parse_attrs_gff3(s):
    d = {}
    for kv in s.strip().strip(';').split(';'):
        if '=' in kv:
            k, v = kv.split('=', 1)
            d[k.strip()] = v.strip()
    return d


def detect_format(path):
    with open(path) as fh:
        for line in fh:
            if line.startswith('##gff-version 3'):
                return 'gff3'
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            if 'transcript_id' in f[8] and '=' not in f[8].split(';')[0]:
                return 'gtf'
            if re.search(r'(^|;)\s*(ID|Parent)=', f[8]):
                return 'gff3'
    return 'gtf'


def read_annotation(path, tool):
    """Return {tid: dict(contig, strand, cds=[(s,e)], stop=[(s,e)], exons=[(s,e)])}."""
    fmt = detect_format(path)
    tx = {}

    def get(tid, contig, strand):
        t = tx.get(tid)
        if t is None:
            t = tx[tid] = {'contig': contig, 'strand': strand, 'cds': [], 'stop': [], 'exons': [],
                           'tool': tool, 'orig_id': tid}
        return t

    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            ftype = f[2]
            if ftype not in ('CDS', 'stop_codon', 'exon'):
                continue
            contig, start, end, strand = f[0], int(f[3]), int(f[4]), f[6]
            if fmt == 'gtf':
                m = GTF_TID.search(f[8])
                if not m:
                    continue
                parents = [m.group(1)]
            else:
                a = parse_attrs_gff3(f[8])
                if 'Parent' not in a:
                    continue
                parents = a['Parent'].split(',')
            key = {'CDS': 'cds', 'stop_codon': 'stop', 'exon': 'exons'}[ftype]
            for p in parents:
                get(p, contig, strand)[key].append((start, end))
    return tx


def merge_intervals(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1] + 1:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def normalise(t):
    """CDS chain with the stop codon always included (AUGUSTUS-style GTFs list it separately)."""
    if not t['cds']:
        return None
    segs = merge_intervals(t['cds'] + t['stop'])
    t['segs'] = tuple(segs)
    t['cds_len'] = sum(e - s + 1 for s, e in segs)
    t['key'] = (t['contig'], t['strand'], t['segs'])
    return t


# ----------------------------------------------------------------------------- RNA-seq hints

MULT = re.compile(r'mult=(\d+)')


def read_intron_hints(paths):
    """bam2hints intron hints -> {(contig, start, end): mult} (strand ignored; summed over lines)."""
    d = defaultdict(int)
    for p in paths:
        with open(p) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f) < 9 or f[2] not in ('intron', 'intronpart'):
                    continue
                m = MULT.search(f[8])
                d[(f[0], int(f[3]), int(f[4]))] += int(m.group(1)) if m else 1
    return d


def read_bedgraph(path, min_depth):
    """bedGraph (e.g. bedtools genomecov -ibam rna.bam -bg -split) -> merged intervals with depth >= min_depth.
    bedGraph is 0-based half-open; converted to 1-based closed."""
    import gzip
    op = gzip.open if path.endswith('.gz') else open
    iv = defaultdict(list)
    with op(path, 'rt') as fh:
        for line in fh:
            if line.startswith(('track', '#', 'browser')):
                continue
            f = line.split()
            if len(f) < 4 or float(f[3]) < min_depth:
                continue
            c, a, b = f[0], int(f[1]) + 1, int(f[2])
            l = iv[c]
            if l and a <= l[-1][1] + 1:          # bedGraph is sorted: merge on the fly
                l[-1][1] = max(l[-1][1], b)
            else:
                l.append([a, b])
    out = {}
    for c, l in iv.items():
        l.sort()
        mg = merge_intervals([tuple(x) for x in l])
        out[c] = ([x for x, _ in mg], [y for _, y in mg])
    return out


def covered_bases(cov, contig, a, b):
    if contig not in cov:
        return 0
    starts, ends = cov[contig]
    i = max(bisect.bisect_right(starts, a) - 1, 0)
    n = 0
    while i < len(starts) and starts[i] <= b:
        lo, hi = max(a, starts[i]), min(b, ends[i])
        if hi >= lo:
            n += hi - lo + 1
        i += 1
    return n


def introns_of(s):
    segs = s['segs']
    return [(s['contig'], segs[i][1] + 1, segs[i + 1][0] - 1) for i in range(len(segs) - 1)]


# ----------------------------------------------------------------------------- full-length transcripts

def read_chains(paths):
    """Transcripts of full-length reads (e.g. the StringTie -L assembly of long reads): their exon chains,
    an index of every intron occurrence and one of every exon. Unstranded transcripts count for both strands."""
    chains = []                   # id -> exons in genomic order
    occ = defaultdict(list)       # (contig, strand, (a, b)) -> [(chain id, intron index)]
    exons = defaultdict(list)     # (contig, strand) -> [(start, end)] of every exon
    for p in paths:
        for t in read_annotation(p, 'chain').values():
            ex = merge_intervals(t['exons'])
            if not ex:
                continue
            cid = len(chains)
            chains.append(ex)
            for sd in ((t['strand'],) if t['strand'] in ('+', '-') else ('+', '-')):
                for i in range(len(ex) - 1):
                    occ[(t['contig'], sd, (ex[i][1] + 1, ex[i + 1][0] - 1))].append((cid, i))
                exons[(t['contig'], sd)].extend(ex)
    index = {}
    for k, l in exons.items():
        l.sort()
        index[k] = ([a for a, _ in l], l, max((b - a for a, b in l), default=0))
    return chains, occ, index


def chain_support(s, chains):
    """'full' if the CDS chain of the model lies in one full-length transcript: its introns are a run of the
    transcript's introns and the terminal CDS segments lie in the flanking exons (a single-exon CDS lies in one
    exon). 'partial' if an intron of the model occurs in such a transcript (a single-exon CDS: overlaps one)
    without a full match, else 'none'. Returns (support, number of transcripts with the full chain)."""
    ch, occ, index = chains
    c, st, segs = s['contig'], s['strand'], s['segs']
    lo, hi = segs[0][0], segs[-1][1]
    ints = [(a, b) for _, a, b in introns_of(s)]
    if ints:
        m, n = len(ints), 0
        for cid, i in occ.get((c, st, ints[0]), []):
            ex = ch[cid]
            if i + m > len(ex) - 1:
                continue
            if all((ex[i + k][1] + 1, ex[i + k + 1][0] - 1) == ints[k] for k in range(m)) \
                    and ex[i][0] <= lo and hi <= ex[i + m][1]:
                n += 1
        if n:
            return 'full', n
        return ('partial', 0) if any((c, st, it) in occ for it in ints) else ('none', 0)
    if (c, st) not in index:
        return 'none', 0
    starts, iv, maxlen = index[(c, st)]
    n = overlap = 0
    i = bisect.bisect_left(starts, lo - maxlen)
    while i < len(iv) and iv[i][0] <= hi:
        a, b = iv[i]
        i += 1
        if b < lo:
            continue
        overlap += 1
        if a <= lo and hi <= b:
            n += 1
    if n:
        return 'full', n
    return ('partial', 0) if overlap else ('none', 0)


# ----------------------------------------------------------------------------- protein evidence (features only)

DMD_TOP = 10   # hits used for median length ratio / terminal offsets


def _median(v):
    v = sorted(v)
    n = len(v)
    return None if n == 0 else (v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2)


DMD_FORMATS = {
    # DIAMOND default (--outfmt 6 without fields)
    'std': 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
    # recommended: appends query/subject lengths to the default. Matches the
    # outfmt of DIAMOND_BLASTP_FULL in modules/diamond.nf.
    'full': 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen',
}


def _is_int(x):
    try:
        int(x)
        return True
    except ValueError:
        return False


def detect_diamond_format(path):
    with open(path) as fh:
        for line in fh:
            x = line.rstrip('\n').split('\t')
            if len(x) < 12:
                continue
            # 14 columns == 'full' (std + qlen slen); 12 == 'std'
            if len(x) >= 14 and all(_is_int(v) for v in x[4:10]) and _is_int(x[12]) and _is_int(x[13]):
                return DMD_FORMATS['full']
            if all(_is_int(v) for v in x[4:10]):
                return DMD_FORMATS['std']
            break
    sys.exit(f'ERROR: cannot recognise the column layout of {path}; give it with --diamond-format')


def read_diamond(path, fmt=None):
    """DIAMOND outfmt 6 in any column order (names as in --outfmt). Hits per query must be in DIAMOND's
    order (best first). Returns {qid: dict(best=row, top=[rows], n_hits, n_taxa)}; rows are dicts."""
    if fmt in DMD_FORMATS:
        fmt = DMD_FORMATS[fmt]
    cols = (fmt or detect_diamond_format(path)).split()
    print(f'DIAMOND {path}: columns = {" ".join(cols)}', file=sys.stderr)
    numeric = {'pident', 'evalue', 'bitscore', 'length', 'mismatch', 'gapopen', 'qlen', 'slen',
               'qstart', 'qend', 'sstart', 'send'}
    out = {}
    cur, top, n, taxa = None, [], 0, set()

    def flush():
        if cur is not None and top:
            out[cur] = {'best': top[0], 'top': top, 'n_hits': n, 'n_taxa': len(taxa)}

    with open(path) as fh:
        for line in fh:
            x = line.rstrip('\n').split('\t')
            if len(x) < len(cols):
                continue
            r = {c: (float(v) if c in numeric else v) for c, v in zip(cols, x)}
            q = r['qseqid']
            if q != cur:
                flush()
                cur, top, n, taxa = q, [], 0, set()
            n += 1
            taxa.add(r['sseqid'].split('_')[0])
            if len(top) < DMD_TOP:
                top.append(r)
        flush()
    return out


def diamond_features(h, qlen_fallback):
    """features from one query's hits; qlen falls back to the CDS length if the table has no qlen."""
    b = h['best']
    qlen = b.get('qlen', qlen_fallback)
    f = {'dmd_bits': b['bitscore'], 'dmd_bits_per_aa': b['bitscore'] / qlen if qlen else '',
         'dmd_pident': b['pident'], 'dmd_n_hits': h['n_hits'], 'dmd_n_taxa': h['n_taxa'],
         'dmd_qcov': (b['qend'] - b['qstart'] + 1) / qlen if ('qstart' in b and qlen) else '',
         'dmd_tcov': (b['send'] - b['sstart'] + 1) / b['slen'] if ('sstart' in b and 'slen' in b) else '',
         'dmd_len_ratio': _median([r.get('qlen', qlen_fallback) / r['slen'] for r in h['top'] if r.get('slen')])
         if 'slen' in b else '',
         # >0: model has extra residues at that end (start too early / stop too late); <0: residues missing
         'dmd_nterm_offset': _median([(r['qstart'] - 1) - (r['sstart'] - 1) for r in h['top']])
         if 'sstart' in b else '',
         'dmd_cterm_offset': _median([(r.get('qlen', qlen_fallback) - r['qend']) - (r['slen'] - r['send'])
                                      for r in h['top']]) if ('slen' in b and 'send' in b) else ''}
    return f


def read_protein_hints(paths):
    """miniprothint hints (hc.gff / miniprothint.gff): intron hints with #proteins in the score column,
    and start/stop hints if present."""
    introns = {}
    stranded = defaultdict(dict)          # (contig, strand) -> {(a, b): n}
    starts, stops = set(), set()
    for p in paths:
        with open(p) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f) < 9:
                    continue
                ft = f[2].lower()
                c, a, b, st = f[0], int(f[3]), int(f[4]), f[6]
                if ft in ('intron', 'intronpart'):
                    try:
                        n = int(float(f[5]))
                    except ValueError:
                        n = 1
                    introns[(c, a, b)] = max(introns.get((c, a, b), 0), n)
                    for sd in ((st,) if st in '+-' else ('+', '-')):
                        stranded[(c, sd)][(a, b)] = max(stranded[(c, sd)].get((a, b), 0), n)
                elif ft in ('start', 'start_codon', 'stop', 'stop_codon'):
                    tgt = starts if ft.startswith('start') else stops
                    for pos in range(a, b + 1):
                        tgt.add((c, pos))
    index = {}
    for k, d in stranded.items():
        iv = sorted((a, b, n) for (a, b), n in d.items())
        index[k] = ([x[0] for x in iv], iv, max((b - a for a, b, _ in iv), default=0))
    return introns, starts, stops, index


def protein_conflicts(s, index):
    """Compare the model with protein intron hints on the same strand within/around its CDS span.
    extra   : hint intron lying completely inside one CDS segment (protein has an intron the model lacks)
    shifted : hint intron overlapping a model intron but with different splice site(s)
    beyond  : hint intron crossing the CDS start/end (protein suggests the gene continues)
    Returns counts and the max protein support (n) per category, plus the number of hint introns in span."""
    key = (s['contig'], s['strand'])
    res = dict(prot_hints_in_span=0, prot_extra=0, prot_extra_max_n=0, prot_shifted=0, prot_shifted_max_n=0,
               prot_beyond=0, prot_beyond_max_n=0)
    if key not in index:
        return res
    starts, iv, maxlen = index[key]
    segs = s['segs']
    lo, hi = segs[0][0], segs[-1][1]
    model_int = set((segs[i][1] + 1, segs[i + 1][0] - 1) for i in range(len(segs) - 1))
    i = bisect.bisect_left(starts, lo - maxlen)
    while i < len(iv) and iv[i][0] <= hi:
        a, b, n = iv[i]
        i += 1
        if b < lo:
            continue
        if a < lo or b > hi:
            res['prot_beyond'] += 1
            res['prot_beyond_max_n'] = max(res['prot_beyond_max_n'], n)
            continue
        res['prot_hints_in_span'] += 1
        if (a, b) in model_int:
            continue
        if any(x <= a and b <= y for x, y in segs):
            res['prot_extra'] += 1
            res['prot_extra_max_n'] = max(res['prot_extra_max_n'], n)
        elif any(a <= y and x <= b for x, y in model_int):
            res['prot_shifted'] += 1
            res['prot_shifted_max_n'] = max(res['prot_shifted_max_n'], n)
        # remaining case: hint intron spanning an exon boundary without overlapping a model intron
        # (e.g. skips a model exon) - counted as shifted too
        else:
            res['prot_shifted'] += 1
            res['prot_shifted_max_n'] = max(res['prot_shifted_max_n'], n)
    return res


def codon_positions(s):
    """first base of start codon and of stop codon (genomic, strand-aware)."""
    lo, hi = s['segs'][0][0], s['segs'][-1][1]
    return (lo, hi - 2) if s['strand'] == '+' else (hi - 2, lo)


def protein_features(S, a):
    """Adds protein evidence columns to every structure. No effect on selection."""
    # miniprot gene models as reference sets (identical CDS)
    mp = {}
    for name, path in (('mp_train', a.miniprot_train), ('mp_rep', a.miniprot_rep)):
        keys = set()
        if path:
            for t in read_annotation(path, name).values():
                t = normalise(t)
                if t is not None:
                    keys.add(t['key'])
        mp[name] = keys if path else None
    for s in S:
        key = (s['contig'], s['strand'], s['segs'])
        for name in ('mp_train', 'mp_rep'):
            s[name + '_match'] = '' if mp[name] is None else int(key in mp[name])

    # miniprothint
    if a.protein_hints:
        pint, pst, psp, pidx = read_protein_hints(a.protein_hints)
        have_codon_hints = bool(pst or psp)
    for s in S:
        if not a.protein_hints:
            s['prot_intron_frac'] = s['prot_intron_min_n'] = s['prot_start'] = s['prot_stop'] = ''
            for c in CONFLICT_COLS:
                s[c] = ''
            continue
        s.update(protein_conflicts(s, pidx))
        ints = introns_of(s)
        ns = [pint.get(i, 0) for i in ints]
        s['prot_intron_frac'] = sum(n > 0 for n in ns) / len(ns) if ns else ''
        s['prot_intron_min_n'] = min(ns) if ns else ''
        st_pos, sp_pos = codon_positions(s)
        s['prot_start'] = int((s['contig'], st_pos) in pst) if have_codon_hints else ''
        s['prot_stop'] = int((s['contig'], sp_pos) in psp) if have_codon_hints else ''

    # DIAMOND per tool, keyed by the tool's transcript id; best member per structure
    dmd = {}
    for tool in TOOLS:
        path = getattr(a, 'diamond_' + tool)
        if path:
            dmd[tool] = read_diamond(path, a.diamond_format)
    cols = ['dmd_bits', 'dmd_bits_per_aa', 'dmd_pident', 'dmd_n_hits', 'dmd_n_taxa', 'dmd_len_ratio',
            'dmd_qcov', 'dmd_tcov', 'dmd_nterm_offset', 'dmd_cterm_offset']
    matched = defaultdict(lambda: [0, 0])
    for s in S:
        best = None
        searched = False
        for m in s['members']:
            if m['tool'] not in dmd:
                continue
            searched = True
            matched[m['tool']][0] += 1
            h = dmd[m['tool']].get(m['orig_id'])
            if h is None:
                h = dmd[m['tool']].get(m['orig_id'] + '.p1')   # common protein-id suffix
            f = diamond_features(h, s['cds_len'] // 3 - 1) if h is not None else None
            if f is not None:
                matched[m['tool']][1] += 1
                if best is None or f['dmd_bits'] > best['dmd_bits']:
                    best = f
        s['dmd_searched'] = int(searched) if dmd else ''
        for c in cols:
            s[c] = '' if best is None else best.get(c, '')
        s['dmd_hit'] = ('' if not dmd else int(best is not None))
    # relative bitscore within cluster
    by_cl = defaultdict(list)
    for s in S:
        by_cl[s['cluster']].append(s)
    for members in by_cl.values():
        mx = max((x['dmd_bits'] for x in members if x['dmd_bits'] != ''), default=None)
        for x in members:
            x['dmd_bits_rel_cluster'] = (x['dmd_bits'] / mx if (mx and x['dmd_bits'] != '') else
                                         (0.0 if mx else ''))
    for tool, (n, k) in matched.items():
        print(f'DIAMOND {tool}: {k}/{n} transcripts have hits ({k / max(n, 1):.1%})', file=sys.stderr)
    if dmd and not any(k for _, k in matched.values()):
        print('WARNING: no DIAMOND query id matched a transcript id - check FASTA headers vs GTF/GFF ids',
              file=sys.stderr)


CONFLICT_COLS = ['prot_hints_in_span', 'prot_extra', 'prot_extra_max_n', 'prot_shifted', 'prot_shifted_max_n',
                 'prot_beyond', 'prot_beyond_max_n']
PROT_COLS = CONFLICT_COLS + ['mp_train_match', 'mp_rep_match', 'prot_intron_frac', 'prot_intron_min_n', 'prot_start',
             'prot_stop', 'dmd_searched', 'dmd_hit', 'dmd_bits', 'dmd_bits_per_aa', 'dmd_bits_rel_cluster',
             'dmd_pident', 'dmd_n_hits', 'dmd_n_taxa', 'dmd_len_ratio', 'dmd_qcov', 'dmd_tcov',
             'dmd_nterm_offset', 'dmd_cterm_offset']


# ----------------------------------------------------------------------------- rules

TIER_TABLE = None          # {(tools, struct_class): tier}, set by --tier-table
RESCUE_TIER = 6


def struct_class(s):
    n = len(s['segs'])
    return 'e1' if n == 1 else ('e2_4' if n <= 4 else 'e5')


def read_tier_table(path):
    """TSV with columns tools, struct_class, tier (tier = 10 x expected exact-CDS precision).
    struct_class '*' = fallback for all classes of that tool combination."""
    t = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        ix = {h: i for i, h in enumerate(hdr)}
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < len(hdr):
                continue
            t[(f[ix['tools']], f[ix['struct_class']])] = float(f[ix['tier']])
    return t


def is_lone(s):
    """single ab initio tool, no transcript tool. Default rules: Tiberius only; with a tier table: either tool."""
    T = s['tools']
    if TIER_TABLE is not None:
        return T in ({'tiberius'}, {'vipsania'})
    return T == {'tiberius'}


def has_abinitio(s):
    return 'tiberius' in s['tools'] or 'vipsania' in s['tools']


def locus_tier(s):
    if TIER_TABLE is not None:
        key = ','.join(sorted(s['tools']))
        v = TIER_TABLE.get((key, struct_class(s)))
        if v is None:
            v = TIER_TABLE.get((key, '*'), 0.0)
        return v
    return locus_tier_default(s)


def locus_tier_default(s):
    T = s['tools']
    tib, vip = 'tiberius' in T, 'vipsania' in T
    tx = 'drusilla' in T or 'td2' in T
    E, L = len(s['segs']), s['cds_len']
    if tib and vip and tx:
        return 10
    if tib and tx:
        return 9
    if tib and vip:
        return 8
    if vip and 'drusilla' in T and 'td2' in T:
        return 7
    if tib:  # Tiberius alone
        return 5 if (E >= 5 or (E >= 3 and L >= 600) or L >= 1200) else 3
    if vip and tx:
        return 4
    return 0


def adjust_tier(s, cfg):
    """RNA / agreement adjustments for lone Tiberius models (tier 3 and 5)."""
    s['tier_base'] = s['tier']
    if not is_lone(s):
        return
    n_int = s['n_introns']
    if cfg.get('chain_rules') and s.get('chain_support') == 'full' and RESCUE_TIER > s['tier']:
        s['tier'], s['adj'] = RESCUE_TIER, 'up_chain_support'
        return
    if cfg['agree_rules'] and n_int > 0 and s['tool_intron_frac'] == 1.0:
        if RESCUE_TIER > s['tier']:
            s['tier'], s['adj'] = RESCUE_TIER, 'up_intron_agreement'
        return
    if not cfg['rna_rules'] or not cfg['have_rna']:
        return
    full = (n_int > 0 and s['rna_frac'] == 1.0) or (n_int == 0 and s['rna_cov'] is not None
                                                   and s['rna_cov'] >= cfg['cov_full'])
    if full and RESCUE_TIER > s['tier']:
        s['tier'], s['adj'] = RESCUE_TIER, 'up_rna_support'


def is_dt(s):
    """Drusilla and TD2 agree, no ab initio tool."""
    T = s['tools']
    return 'drusilla' in T and 'td2' in T and 'tiberius' not in T and 'vipsania' not in T


def rank_score(s, rna_penalty=0.0, chain_bonus=0.0):
    if TIER_TABLE is not None:
        r = s['tier']                     # the table already rates Drusilla+TD2 agreement
    else:
        r = 5.5 if is_dt(s) else s['tier']   # Drusilla+TD2 outranks lone Tiberius (3/5) inside a locus
    # inside a locus, prefer candidates whose introns are all RNA-supported
    if rna_penalty and s.get('n_introns') and s.get('rna_frac') not in ('', None) and s['rna_frac'] < 1:
        r -= rna_penalty
    # and those whose CDS chain lies in one full-length transcript
    if chain_bonus and s.get('chain_support') == 'full':
        r += chain_bonus
    return r + 0.01 * len(s['tools'])


def protein_adjust(s, lo, hi):
    """Drop lone-Tiberius models (base tier 3/5) whose best DIAMOND hit covers only part of a protein
    (length ratio outside [lo, hi]) unless an identical miniprot model or full protein intron support
    confirms them. Overrides RNA/agreement upgrades (measured on Mimulus: F1 75.3 -> 75.6)."""
    if not is_lone(s) or s.get('dmd_hit') != 1:
        return
    lr = s.get('dmd_len_ratio')
    if lr in ('', None) or lo <= lr <= hi:
        return
    if s.get('mp_train_match') == 1 or s.get('mp_rep_match') == 1:
        return
    if s['n_introns'] > 0 and s.get('prot_intron_frac') == 1.0:
        return
    s['tier'], s['adj'] = 2, 'down_partial_protein'



def cluster(structs):
    order = sorted(structs, key=lambda s: (s['contig'], s['strand'], s['segs'][0][0]))
    cid, cur, cend = -1, None, -1
    for s in order:
        a, b = s['segs'][0][0], s['segs'][-1][1]
        if (s['contig'], s['strand']) != cur or a > cend:
            cid += 1
            cur, cend = (s['contig'], s['strand']), b
        else:
            cend = max(cend, b)
        s['cluster'] = cid
    return order


def stop_pos(s):
    return s['segs'][-1][1] if s['strand'] == '+' else s['segs'][0][0]


def decide(members, min_locus, iso_min_tier, iso_dt, iso_rna=False, min_mult=2, iso_chain_min_tier=None,
           iso_dt_chains=False, **_):
    members.sort(key=lambda s: (-s['score'], s['segs']))
    open_locus = max((s.get('open_tier', s['tier']) if has_abinitio(s) else 0)
                     for s in members) >= min_locus
    if open_locus and members[0].get('no_primary'):
        # a demoted model still marks the locus as a gene, but the primary must be another candidate that
        # would itself rank at the locus threshold (e.g. Drusilla+TD2 agreement); otherwise the locus is dropped
        alt = [s for s in members if not s.get('no_primary') and (5.5 if (is_dt(s) and TIER_TABLE is None) else s['tier']) >= min_locus]
        if alt:
            members.remove(alt[0])
            members.insert(0, alt[0])
        else:
            open_locus = False
    prim = members[0]
    prim_introns = set(introns_of(prim))
    for s in members:
        nov = [i for i in introns_of(s) if i not in prim_introns]
        s['novel_introns'] = len(nov)
        s['novel_rna_frac'] = (sum(s['_imult'].get(i, 0) >= min_mult for i in nov) / len(nov)) if nov else ''
    for r, s in enumerate(members):
        s['rank'] = r + 1
        if not open_locus:
            s['keep'], s['why'] = False, 'locus_below_min'
        elif r == 0:
            s['keep'], s['why'] = True, 'primary'
        elif s['tier'] >= iso_min_tier:
            s['keep'], s['why'] = True, 'isoform_tier'
        elif iso_chain_min_tier is not None and s.get('chain_support') == 'full' and (
                5.5 if (is_dt(s) and TIER_TABLE is None) else s['tier']) >= iso_chain_min_tier:
            s['keep'], s['why'] = True, 'isoform_chain'
        elif is_dt(s) and (iso_dt == 'all' or (iso_dt == 'filtered' and (
                stop_pos(s) == stop_pos(prim) or s['cds_len'] >= 0.95 * prim['cds_len']))) and (
                not iso_rna or s['novel_rna_frac'] in ('', 1.0)) and (
                not iso_dt_chains or s.get('chain_support') == 'full'):
            s['keep'], s['why'] = True, 'isoform_drusilla_td2'
        else:
            s['keep'], s['why'] = False, 'isoform_rejected'


# ----------------------------------------------------------------------------- report

REPORT_COLS = ['cluster', 'contig', 'strand', 'cds_start', 'cds_end', 'n_cds_exons', 'cds_len', 'tools', 'locus_tier',
               'rank_score', 'rank_in_cluster', 'kept', 'reason', 'output_id', 'base_tier', 'adjustment', 'n_introns',
               'rna_introns', 'rna_frac', 'rna_min_mult', 'rna_cov', 'tool_intron_frac', 'novel_introns',
               'novel_rna_frac', 'chain_support', 'chain_n'] + PROT_COLS + ['member_ids']


def rnd(x):
    return '' if x is None or x == '' else round(x, 3)


def report_row(s):
    """the --report values of a structure, in REPORT_COLS order (also the input of the classifier features)"""
    return [s['cluster'], s['contig'], s['strand'], s['segs'][0][0], s['segs'][-1][1], len(s['segs']),
            s['cds_len'], ','.join(sorted(s['tools'])), s['tier'], round(s['score'], 2), s['rank'],
            int(s['keep']), s['why'], s.get('out_tid', ''), s['tier_base'], s['adj'] or '.',
            s['n_introns'], s['rna_introns'], rnd(s['rna_frac']), s['rna_min_mult'], rnd(s['rna_cov']),
            rnd(s['tool_intron_frac']), s['novel_introns'], rnd(s['novel_rna_frac']),
            s['chain_support'], s['chain_n'],
            *[rnd(s[c]) if isinstance(s[c], float) else s[c] for c in PROT_COLS],
            ','.join(f"{m['tool']}:{m['orig_id']}" for m in s['members'])]


# ----------------------------------------------------------------------------- learned selection

CLASSIFIER_FORMAT = 'combine_genesets-classifier-1'
NAN = float('nan')
# defaults of --chain-logit-full / --chain-logit-partial. Calibrated on the H. sapiens runs with Iso-Seq (long-read
# StringTie -L transcripts, 2026-11-06): full support is worth +1.5 to +1.7 log-odds beyond structures without
# support (maximum likelihood), 'partial' support is worth nothing beyond 'none'. gffcompare: transcript F1
# +1.2 / +1.7 (runs S / M) over the classifier without chains at 1.5; 2.0 gives +0.2 more at -1.3 precision.
CHAIN_LOGIT_FULL = 1.5
CHAIN_LOGIT_PARTIAL = 0.0


def load_classifier(path):
    """JSON export of the selection model: features (report-derived, see classifier_features), thresholds,
    baseline log-odds and the trees (per tree: feat, thr, left, right, leaf, value, nan_left per node)."""
    with open(path) as fh:
        clf = json.load(fh)
    if clf.get('format') != CLASSIFIER_FORMAT:
        sys.exit(f'{path}: not a combine_genesets classifier (format {clf.get("format")!r}, '
                 f'expected {CLASSIFIER_FORMAT!r})')
    for k in ('features', 'trees', 'baseline', 't_locus', 't_iso'):
        if k not in clf:
            sys.exit(f'{path}: classifier file lacks {k!r}')
    return clf


def _num(v):
    """report value -> float; '' / None / non-numeric -> NaN (pandas.to_numeric(errors='coerce') on the TSV)"""
    if v is None or v == '':
        return NAN
    try:
        return float(v)
    except (TypeError, ValueError):
        return NAN


def _isnan(x):
    return x != x


def _agg(vals, fn):
    """aggregate over the non-NaN values; NaN if there is none (pandas groupby max / min / median)"""
    v = [x for x in vals if not _isnan(x)]
    return fn(v) if v else NAN


CLF_PASS = ['prot_intron_frac', 'prot_intron_min_n', 'prot_start', 'prot_stop', 'prot_hints_in_span', 'prot_extra',
            'prot_shifted', 'prot_beyond', 'dmd_hit', 'dmd_bits_per_aa', 'dmd_bits_rel_cluster', 'dmd_pident',
            'dmd_n_taxa', 'dmd_qcov', 'dmd_tcov', 'dmd_nterm_offset', 'dmd_cterm_offset', 'rna_frac', 'rna_cov',
            'tool_intron_frac', 'novel_introns', 'novel_rna_frac', 'base_tier', 'n_cds_exons', 'cds_len', 'n_introns']


def classifier_features(S, by_cl):
    """One feature dict per structure (S order), computed from the report values exactly as the training code
    (select_classifier.py features()) computes them from a --keep-all report: missing evidence is NaN, locus
    context comes from all structures of the cluster, RNA multiplicity is normalised by the species' median
    over the four-tool structures."""
    F = []
    for s in S:
        r = dict(zip(REPORT_COLS, report_row(s)))
        f = {'by_' + t: int(t in r['tools']) for t in TOOLS}
        f['ntools'] = sum(f['by_' + t] for t in TOOLS)
        f['dt'] = int(f['by_drusilla'] == 1 and f['by_td2'] == 1 and f['by_tiberius'] == 0 and f['by_vipsania'] == 0)
        for c in CLF_PASS:
            f[c] = _num(r[c])
        f['log_cds_len'] = math.log10(max(f['cds_len'], 1.0))
        mm = _num(r['rna_min_mult'])
        f['log_mult'] = NAN if _isnan(mm) else math.log1p(mm)
        f['mp_match'] = _agg([_num(r['mp_train_match']), _num(r['mp_rep_match'])], max)
        nh, nt = _num(r['dmd_n_hits']), f['dmd_n_taxa']
        f['hits_per_taxon'] = NAN if (_isnan(nh) or _isnan(nt)) else (nh / nt if nt else (math.inf if nh > 0 else NAN))
        lr = _num(r['dmd_len_ratio'])
        f['lr_dev'] = NAN if _isnan(lr) else abs(math.log(min(max(lr, 0.01), 100.0)))
        F.append(f)
    ref = _agg([f['log_mult'] for f in F if f['ntools'] == 4], _median)
    ref = 0.0 if _isnan(ref) else ref
    for members in by_cl.values():
        fs = [F[s['_ix']] for s in members]
        maxtier = _agg([f['base_tier'] for f in fs], max)
        max_ntools = max(f['ntools'] for f in fs)
        max_len = _agg([f['cds_len'] for f in fs], max)
        max_rf = _agg([f['rna_frac'] for f in fs], max)
        min_lr = _agg([f['lr_dev'] for f in fs], min)
        max_pif = _agg([f['prot_intron_frac'] for f in fs], max)
        max_hit = _agg([f['dmd_hit'] for f in fs], max)
        for f in fs:
            f['cl_n'] = len(fs)
            f['cl_maxtier'] = maxtier
            f['tier_gap'] = f['base_tier'] - maxtier
            f['cl_max_ntools'] = max_ntools
            f['ntools_gap'] = f['ntools'] - max_ntools
            f['len_rel_cluster'] = f['cds_len'] / max_len
            f['log_mult_rel'] = f['log_mult'] - ref
            f['rna_frac_rel'] = (1.0 if _isnan(f['rna_frac']) else f['rna_frac']) - (1.0 if _isnan(max_rf) else max_rf)
            f['lr_dev_rel'] = f['lr_dev'] - min_lr
            f['pif_rel'] = (0.0 if _isnan(f['prot_intron_frac']) else f['prot_intron_frac']) - \
                           (0.0 if _isnan(max_pif) else max_pif)
            f['cl_has_hit'] = max_hit
    return F


def classifier_predict(clf, F):
    """P per structure: sigmoid(baseline + sum of the leaf values); a node sends x <= thr left, NaN as nan_left."""
    try:
        import numpy as np
    except ImportError:
        sys.exit('--classifier needs numpy')
    names = clf['features']
    missing = [c for c in names if c not in F[0]] if F else []
    if missing:
        sys.exit(f'classifier features not computed by this version of combine_genesets.py: {missing}')
    X = np.array([[f[c] for c in names] for f in F], dtype=float).reshape(len(F), len(names))
    n = len(F)
    raw = np.full(n, float(clf['baseline']))
    rows = np.arange(n)
    for t in clf['trees']:
        feat, left, right = (np.asarray(t[k], dtype=int) for k in ('feat', 'left', 'right'))
        thr, value = (np.asarray(t[k], dtype=float) for k in ('thr', 'value'))
        leaf, nan_left = (np.asarray(t[k], dtype=bool) for k in ('leaf', 'nan_left'))
        node = np.zeros(n, dtype=int)
        active = ~leaf[node]
        while active.any():
            x = X[rows, feat[node]]
            go_left = np.where(np.isnan(x), nan_left[node], x <= thr[node])
            node = np.where(active, np.where(go_left, left[node], right[node]), node)
            active = ~leaf[node]
        raw += value[node]
    return (1.0 / (1.0 + np.exp(-raw))).tolist()


def chain_adjust(p, support, w_full, w_partial):
    """P corrected for full-length transcript support (the model was trained without --chains): the log-odds
    of a structure whose CDS chain lies in one full-length transcript move by w_full, those of a 'partial'
    structure by w_partial. 'none' and '' (no --chains) are unchanged."""
    w = w_full if support == 'full' else (w_partial if support == 'partial' else 0.0)
    if not w:
        return p
    q = min(max(p, 1e-6), 1 - 1e-6)
    return 1.0 / (1.0 + math.exp(-(math.log(q / (1 - q)) + w)))


def classifier_select(S, by_cl, p, t_locus, t_iso, unselectable):
    """Per locus: candidates ordered by P (ties: input order); the first is primary if its P >= t_locus, the
    others are isoforms if P >= t_iso. Structures in `unselectable` (indices into S) are never selected and
    never primary (--drop-lone-single-exon)."""
    for s, pi in zip(S, p):
        s['p_clf'] = float(pi)
    for members in by_cl.values():
        members.sort(key=lambda s: (-(-1.0 if s['_ix'] in unselectable else s['p_clf']), s['_ix']))
        prim = members[0]
        open_locus = prim['_ix'] not in unselectable and prim['p_clf'] >= t_locus
        for r, s in enumerate(members):
            s['rank'] = r + 1
            if s['_ix'] in unselectable:
                s['keep'], s['why'] = False, 'dropped_lone_single_exon'
            elif not open_locus:
                s['keep'], s['why'] = False, 'locus_below_min'
            elif r == 0:
                s['keep'], s['why'] = True, 'primary'
            elif s['p_clf'] >= t_iso:
                s['keep'], s['why'] = True, 'isoform_classifier'
            else:
                s['keep'], s['why'] = False, 'isoform_rejected'


# ----------------------------------------------------------------------------- output

def gtf_line(c, src, ft, s, e, st, frame, gid, tid, extra=''):
    return f'{c}\t{src}\t{ft}\t{s}\t{e}\t.\t{st}\t{frame}\tgene_id "{gid}"; transcript_id "{tid}";{extra}\n'


def gff3_line(c, src, ft, s, e, st, frame, attrs):
    return f'{c}\t{src}\t{ft}\t{s}\t{e}\t.\t{st}\t{frame}\t{attrs}\n'


def _cds_wo_stop(segs, strand):
    """CDS blocks with the trailing 3 nt (stop codon) stripped off, in genomic order."""
    ordered = segs if strand == '+' else segs[::-1]
    stop_left = 3
    out = [list(x) for x in ordered]
    while stop_left > 0 and out:
        last = out[-1]
        ln = last[1] - last[0] + 1
        if ln <= stop_left:
            stop_left -= ln
            out.pop()
        else:
            if strand == '+':
                last[1] -= stop_left
            else:
                last[0] += stop_left
            stop_left = 0
    return out


def write_transcript(out, s, gid, tid, rep):
    c, st = s['contig'], s['strand']
    segs = list(s['segs'])
    exons = merge_intervals(rep['exons']) if rep['exons'] else list(segs)
    exons = merge_intervals(exons + segs)
    lo, hi = exons[0][0], exons[-1][1]
    extra = f' tools "{",".join(sorted(s["tools"]))}"; tier "{s["tier"]}"; source_id "{rep["orig_id"]}";'
    out.write(gtf_line(c, 'combined', 'transcript', lo, hi, st, '.', gid, tid, extra))
    for a, b in exons:
        out.write(gtf_line(c, 'combined', 'exon', a, b, st, '.', gid, tid))
    cds_wo_stop = _cds_wo_stop(segs, st)
    done = 0
    for a, b in cds_wo_stop:
        frame = (3 - done % 3) % 3
        out.write(gtf_line(c, 'combined', 'CDS', a, b, st, frame, gid, tid))
        done += b - a + 1
    if st == '+':
        out.write(gtf_line(c, 'combined', 'start_codon', segs[0][0], segs[0][0] + 2, st, 0, gid, tid))
        out.write(gtf_line(c, 'combined', 'stop_codon', segs[-1][1] - 2, segs[-1][1], st, 0, gid, tid))
    else:
        out.write(gtf_line(c, 'combined', 'start_codon', segs[-1][1] - 2, segs[-1][1], st, 0, gid, tid))
        out.write(gtf_line(c, 'combined', 'stop_codon', segs[0][0], segs[0][0] + 2, st, 0, gid, tid))


def write_gene_gff3(out, c, st, gid, lo, hi):
    out.write(gff3_line(c, 'combined', 'gene', lo, hi, st, '.', f'ID={gid}'))


def write_transcript_gff3(out, s, gid, tid, rep):
    c, st = s['contig'], s['strand']
    segs = list(s['segs'])
    exons = merge_intervals(rep['exons']) if rep['exons'] else list(segs)
    exons = merge_intervals(exons + segs)
    lo, hi = exons[0][0], exons[-1][1]
    attrs = (f'ID={tid};Parent={gid};tools={",".join(sorted(s["tools"]))};'
             f'tier={s["tier"]};source_id={rep["orig_id"]}')
    out.write(gff3_line(c, 'combined', 'mRNA', lo, hi, st, '.', attrs))
    for i, (a, b) in enumerate(exons, 1):
        out.write(gff3_line(c, 'combined', 'exon', a, b, st, '.', f'ID={tid}.exon{i};Parent={tid}'))
    cds_wo_stop = _cds_wo_stop(segs, st)
    done = 0
    for i, (a, b) in enumerate(cds_wo_stop, 1):
        frame = (3 - done % 3) % 3
        out.write(gff3_line(c, 'combined', 'CDS', a, b, st, frame, f'ID={tid}.cds{i};Parent={tid}'))
        done += b - a + 1
    if st == '+':
        out.write(gff3_line(c, 'combined', 'start_codon', segs[0][0], segs[0][0] + 2, st, 0,
                            f'ID={tid}.start;Parent={tid}'))
        out.write(gff3_line(c, 'combined', 'stop_codon', segs[-1][1] - 2, segs[-1][1], st, 0,
                            f'ID={tid}.stop;Parent={tid}'))
    else:
        out.write(gff3_line(c, 'combined', 'start_codon', segs[-1][1] - 2, segs[-1][1], st, 0,
                            f'ID={tid}.start;Parent={tid}'))
        out.write(gff3_line(c, 'combined', 'stop_codon', segs[0][0], segs[0][0] + 2, st, 0,
                            f'ID={tid}.stop;Parent={tid}'))


# ----------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_argument_group('gene sets and output')
    for t in TOOLS:
        g.add_argument(f'--{t}', required=True, help=f'{t} gene set (GTF or GFF3)')
    g.add_argument('-o', '--out', required=True,
                   help='output annotation; GFF3 when the extension is .gff3/.gff, GTF otherwise')
    g.add_argument('--report', help='TSV with every unique CDS structure, its features, tier and decision')
    g.add_argument('--prefer', default='drusilla,td2,tiberius,vipsania',
                   help='tool order for the representative record (UTRs, IDs) of a merged structure')

    g = ap.add_argument_group('selection')
    g.add_argument('--min-locus', type=float, default=4,
                   help='keep a locus if its best model reaches this tier [4]; with --tier-table: 10 x precision')
    g.add_argument('--tier-table', help='clade-specific tier table from calibrate_tiers.py (replaces built-in tiers)')
    g.add_argument('--iso-min-tier', type=float, default=8, help='min tier for additional isoforms [8]')
    g.add_argument('--iso-dt', choices=['none', 'filtered', 'all'], default='filtered',
                   help='Drusilla+TD2-agreed structures as additional isoforms [filtered]')
    g.add_argument('--keep-all', action='store_true',
                   help='keep every unique structure and switch off all rules (for labelling / calibration runs)')

    g = ap.add_argument_group('RNA-seq evidence')
    g.add_argument('--introns', nargs='*', default=[], help='intron hints with mult= (bam2hints)')
    g.add_argument('--coverage', help='bedGraph of read depth (bedtools genomecov -ibam x.bam -bg -split), .gz ok')
    g.add_argument('--min-mult', type=int, default=2, help='min mult for a supported intron [2]')
    g.add_argument('--min-depth', type=float, default=2, help='min read depth counted as covered [2]')
    g.add_argument('--cov-full', type=float, default=0.9,
                   help='CDS coverage for full support of a single-exon model [0.9]')
    g.add_argument('--rank-rna-penalty', type=float, default=2.0,
                   help='rank penalty inside a locus for a model with an RNA-unsupported intron [2]')
    g.add_argument('--no-rna-rules', action='store_true', help='compute RNA features but do not use them')
    g.add_argument('--no-agree-rules', action='store_true', help='no upgrade by intron-level tool agreement')
    g.add_argument('--no-drop-lone-rna-conflict', dest='drop_lone_rna_conflict', action='store_false',
                   help='keep multi-exon models of one ab initio tool whose introns are only partly RNA-supported '
                        '(dropped by default when --introns is given)')

    g = ap.add_argument_group('full-length transcripts')
    g.add_argument('--chains', nargs='*', default=[],
                   help='transcripts of full-length reads with exon lines (GTF/GFF3), e.g. the StringTie -L '
                        'assembly of the long reads')
    g.add_argument('--iso-chain-min-tier', type=float, default=4, metavar='TIER',
                   help='keep an additional isoform whose CDS chain lies in one full-length transcript if its '
                        'tier (Drusilla+TD2: 5.5) reaches this [4]')
    g.add_argument('--iso-dt-chains', action='store_true',
                   help='Drusilla+TD2 isoforms need a full-length transcript holding their CDS chain')
    g.add_argument('--rank-chain-bonus', type=float, default=1.0,
                   help='rank bonus inside a locus for a model whose CDS chain lies in one full-length transcript [1]')
    g.add_argument('--no-chain-rules', action='store_true',
                   help='compute the full-length transcript features but do not use them')

    g = ap.add_argument_group('protein evidence')
    for t in TOOLS:
        g.add_argument(f'--diamond-{t}', help=f'DIAMOND outfmt-6 table for the {t} proteins')
    g.add_argument('--diamond-format', help="DIAMOND column names as given to --outfmt 6, or 'std' / 'full'; "
                                            'auto-detected if omitted')
    g.add_argument('--miniprot-train', help='miniprot_trainingGenes.gff')
    g.add_argument('--miniprot-rep', help='miniprot_representatives.gff')
    g.add_argument('--protein-hints', nargs='*', default=[], help='miniprothint hc.gff')
    g.add_argument('--rescue-lone-miniprot', type=float, default=4.5, metavar='TIER',
                   help='raise lone Tiberius/Vipsania models identical to a miniprot model to this tier [4.5]')
    g.add_argument('--no-rescue-lone-miniprot', dest='rescue_on', action='store_false', help='switch the rescue off')
    g.add_argument('--no-prot-drop-partial', dest='prot_drop_partial', action='store_false',
                   help='keep lone ab initio models with a partial-length DIAMOND hit')
    g.add_argument('--full-len-min', type=float, default=0.8, help='min query/subject length ratio, full hit [0.8]')
    g.add_argument('--full-len-max', type=float, default=1.25, help='max query/subject length ratio, full hit [1.25]')

    g = ap.add_argument_group('clade-specific')
    g.add_argument('--drop-lone-single-exon', metavar='TOOLS',
                   help='comma-separated tools (insects: tiberius): drop single-exon models of that tool alone '
                        'unless identical to a miniprot model')

    g = ap.add_argument_group('learned selection')
    g.add_argument('--classifier', metavar='MODEL.json',
                   help='select by the probability of this model instead of the tier rules (candidate pool as with '
                        '--keep-all; see the module docstring)')
    g.add_argument('--unknown-clade', action='store_true',
                   help='the clade is outside the training clades of the model: use its conservative locus threshold')
    g.add_argument('--t-locus', type=float, help='locus threshold on P, overrides the model file')
    g.add_argument('--t-iso', type=float, help='isoform threshold on P, overrides the model file')
    g.add_argument('--chain-logit-full', type=float, default=CHAIN_LOGIT_FULL, metavar='W',
                   help='with --chains: add W to the log-odds of P for a structure whose CDS chain lies in one '
                        f'full-length transcript [{CHAIN_LOGIT_FULL}]')
    g.add_argument('--chain-logit-partial', type=float, default=CHAIN_LOGIT_PARTIAL, metavar='W',
                   help=f"with --chains: the same for 'partial' chain support [{CHAIN_LOGIT_PARTIAL}]")
    a = ap.parse_args()

    clf = None
    clf_drop_se = set()
    if a.classifier:
        if a.tier_table:
            sys.exit('--classifier and --tier-table exclude each other')
        clf = load_classifier(a.classifier)
        clf_drop_se = set(a.drop_lone_single_exon.split(',')) if a.drop_lone_single_exon else set()
        a.keep_all = True                          # candidate pool and features as in the training runs
    if a.keep_all:
        a.min_locus, a.iso_min_tier, a.iso_dt = 0, 0, 'all'
        a.no_rna_rules = a.no_agree_rules = a.no_chain_rules = True
        a.drop_lone_rna_conflict = a.prot_drop_partial = a.rescue_on = False
        a.drop_lone_single_exon = None
    a.rescue_lone_miniprot = a.rescue_lone_miniprot if a.rescue_on else None
    if a.rescue_lone_miniprot is not None and float(a.rescue_lone_miniprot).is_integer():
        a.rescue_lone_miniprot = int(a.rescue_lone_miniprot)
    global TIER_TABLE
    if a.tier_table:
        TIER_TABLE = read_tier_table(a.tier_table)
        print(f'tier table: {a.tier_table} ({len(TIER_TABLE)} cells)', file=sys.stderr)
    cfg = dict(min_locus=a.min_locus, iso_min_tier=a.iso_min_tier, iso_dt=a.iso_dt)

    prefer = {t: i for i, t in enumerate(a.prefer.split(','))}
    structs = {}
    n_in = defaultdict(int)
    for tool in TOOLS:
        for t in read_annotation(getattr(a, tool), tool).values():
            t = normalise(t)
            if t is None:
                continue
            n_in[tool] += 1
            s = structs.get(t['key'])
            if s is None:
                s = structs[t['key']] = {'contig': t['contig'], 'strand': t['strand'], 'segs': t['segs'],
                                         'cds_len': t['cds_len'], 'tools': set(), 'members': []}
            s['tools'].add(tool)
            s['members'].append(t)

    imult = read_intron_hints(a.introns) if a.introns else {}
    cov = read_bedgraph(a.coverage, a.min_depth) if a.coverage else None
    have_rna = bool(a.introns) or cov is not None
    chains = read_chains(a.chains) if a.chains else None
    chain_rules = chains is not None and not a.no_chain_rules
    rcfg = dict(rna_rules=not a.no_rna_rules, agree_rules=not a.no_agree_rules, have_rna=have_rna,
                cov_full=a.cov_full, chain_rules=chain_rules)
    cfg.update(iso_rna=bool(a.introns) and not a.no_rna_rules, min_mult=a.min_mult,
               iso_chain_min_tier=a.iso_chain_min_tier if chain_rules else None,
               iso_dt_chains=a.iso_dt_chains and chain_rules)
    rna_pen = a.rank_rna_penalty if a.introns else 0.0
    chain_bonus = a.rank_chain_bonus if chain_rules else 0.0

    S = list(structs.values())
    intron_tools = defaultdict(set)
    for s in S:
        for i in introns_of(s):
            intron_tools[i] |= s['tools']
    for s in S:
        ints = introns_of(s)
        s['_imult'] = imult
        s['n_introns'] = len(ints)
        s['rna_introns'] = sum(imult.get(i, 0) >= a.min_mult for i in ints)
        s['rna_frac'] = s['rna_introns'] / len(ints) if ints else ''
        s['rna_min_mult'] = min((imult.get(i, 0) for i in ints), default='') if have_rna else ''
        s['rna_cov'] = (sum(covered_bases(cov, s['contig'], x, y) for x, y in s['segs']) / s['cds_len']
                        if cov is not None else None)
        s['tool_intron_frac'] = (sum(bool(intron_tools[i] - s['tools']) for i in ints) / len(ints)) if ints else ''
        s['chain_support'], s['chain_n'] = chain_support(s, chains) if chains is not None else ('', '')
        s['adj'] = ''
        s['tier'] = locus_tier(s)
        adjust_tier(s, rcfg)
        s['score'] = rank_score(s, rna_pen, chain_bonus)
        s['rep'] = min(s['members'], key=lambda m: (prefer.get(m['tool'], 99), m['orig_id']))

    S = cluster(S)
    for i, s in enumerate(S):
        s['_ix'] = i
    protein_features(S, a)
    if a.prot_drop_partial and any(getattr(a, 'diamond_' + t) for t in TOOLS):
        for s in S:
            protein_adjust(s, a.full_len_min, a.full_len_max)
    drop_se = set(a.drop_lone_single_exon.split(',')) if a.drop_lone_single_exon else set()
    for s in S:
        mp_ok = s.get('mp_train_match') == 1 or s.get('mp_rep_match') == 1
        single_tool = next(iter(s['tools'])) if len(s['tools']) == 1 else None
        rna_conflict = (a.drop_lone_rna_conflict and bool(a.introns) and single_tool in ('tiberius', 'vipsania')
                        and s['n_introns'] > 0 and s['rna_introns'] >= 1 and s['rna_frac'] != '' and s['rna_frac'] < 1)
        if drop_se and single_tool in drop_se and len(s['segs']) == 1 and not mp_ok:
            s['open_tier'], s['no_primary'] = s['tier'], True
            s['tier'], s['adj'] = 2, 'down_lone_single_exon'
        elif rna_conflict and not mp_ok:
            s['open_tier'], s['no_primary'] = s['tier'], True
            s['tier'], s['adj'] = 2, 'down_lone_rna_conflict'
        elif a.rescue_lone_miniprot is not None and single_tool in ('tiberius', 'vipsania') and mp_ok \
                and s['tier'] < a.rescue_lone_miniprot:
            s['tier'], s['adj'] = a.rescue_lone_miniprot, 'up_lone_miniprot'
        s['score'] = rank_score(s, rna_pen, chain_bonus)
    by_cl = defaultdict(list)
    for s in S:
        by_cl[s['cluster']].append(s)
    for members in by_cl.values():
        decide(members, **cfg)

    t_locus = t_iso = None
    if clf is not None:
        t_locus = a.t_locus if a.t_locus is not None else \
            float(clf.get('t_locus_unknown_clade', clf['t_locus']) if a.unknown_clade else clf['t_locus'])
        t_iso = a.t_iso if a.t_iso is not None else float(clf['t_iso'])
        unselectable = {s['_ix'] for s in S
                        if clf_drop_se and len(s['tools']) == 1 and next(iter(s['tools'])) in clf_drop_se
                        and len(s['segs']) == 1 and not (s.get('mp_train_match') == 1 or s.get('mp_rep_match') == 1)}
        p = classifier_predict(clf, classifier_features(S, by_cl))
        if chains is not None:
            p = [chain_adjust(pi, s['chain_support'], a.chain_logit_full, a.chain_logit_partial) for pi, s in zip(p, S)]
        classifier_select(S, by_cl, p, t_locus, t_iso, unselectable)

    def cl_pos(cid):
        m = by_cl[cid]
        return (m[0]['contig'], min(x['segs'][0][0] for x in m))

    # Output format from the extension: .gff3 / .gff -> GFF3; everything else -> GTF.
    out_lower = a.out.lower()
    emit_gff3 = out_lower.endswith('.gff3') or out_lower.endswith('.gff')

    gene_no = n_kept = 0
    with open(a.out, 'w') as out:
        if emit_gff3:
            out.write('##gff-version 3\n')
        for cid in sorted(by_cl, key=cl_pos):
            sel = sorted((s for s in by_cl[cid] if s['keep']), key=lambda s: s['rank'])
            if not sel:
                continue
            gene_no += 1
            gid = f'g{gene_no}'
            c, st = sel[0]['contig'], sel[0]['strand']
            lo = min(min([e[0] for e in s['rep']['exons']] + [s['segs'][0][0]]) for s in sel)
            hi = max(max([e[1] for e in s['rep']['exons']] + [s['segs'][-1][1]]) for s in sel)
            if emit_gff3:
                write_gene_gff3(out, c, st, gid, lo, hi)
            else:
                out.write(f'{c}\tcombined\tgene\t{lo}\t{hi}\t.\t{st}\t.\tgene_id "{gid}";\n')
            for k, s in enumerate(sel, 1):
                s['out_tid'] = f'{gid}.t{k}'
                if emit_gff3:
                    write_transcript_gff3(out, s, gid, s['out_tid'], s['rep'])
                else:
                    write_transcript(out, s, gid, s['out_tid'], s['rep'])
                n_kept += 1

    if a.report:
        with open(a.report, 'w') as rp:
            rp.write('\t'.join(REPORT_COLS + (['p_classifier'] if clf is not None else [])) + '\n')
            for s in S:
                row = report_row(s) + ([round(s['p_clf'], 4)] if clf is not None else [])
                rp.write('\t'.join(map(str, row)) + '\n')

    adj = defaultdict(int)
    for s in S:
        adj[s['adj'] or 'none'] += 1
    why = defaultdict(int)
    for s in S:
        why[s['why']] += 1
    print('settings: ' + ' '.join(f'{k}={v}' for k, v in cfg.items()) +
          f' drop_lone_rna_conflict={a.drop_lone_rna_conflict and bool(a.introns)}'
          f' rescue_lone_miniprot={a.rescue_lone_miniprot} prot_drop_partial={a.prot_drop_partial}'
          f' drop_lone_single_exon={a.drop_lone_single_exon} keep_all={a.keep_all}', file=sys.stderr)
    if clf is not None:
        print(f'classifier: {a.classifier} ({len(clf["trees"])} trees, {len(clf["features"])} features) '
              f't_locus={t_locus} t_iso={t_iso} unknown_clade={a.unknown_clade} '
              f'drop_lone_single_exon={",".join(sorted(clf_drop_se)) or None}'
              + (f' chain_logit_full={a.chain_logit_full} chain_logit_partial={a.chain_logit_partial}'
                 if chains is not None else ''), file=sys.stderr)
    print('input transcripts with CDS: ' + ', '.join(f'{t}={n_in[t]}' for t in TOOLS), file=sys.stderr)
    print(f'unique CDS structures: {len(S)}; clusters: {len(by_cl)}; kept {n_kept} transcripts in '
          f'{gene_no} genes', file=sys.stderr)
    if a.introns:
        hi = [s for s in S if s['tier_base'] >= 9 and s['n_introns'] > 0]
        n_i = sum(s['n_introns'] for s in hi)
        n_hit = sum(sum(imult.get(i, 0) > 0 for i in introns_of(s)) for s in hi)
        n_sup = sum(s['rna_introns'] for s in hi)
        print(f'intron-hint sanity check (tier>=9 models, {n_i} introns): any hint {n_hit / max(n_i, 1):.1%}, '
              f'mult>={a.min_mult} {n_sup / max(n_i, 1):.1%}', file=sys.stderr)
        if n_i and n_hit / n_i < 0.5:
            print('WARNING: <50% of introns of well-supported models match a hint - check contig names and '
                  'intron coordinates (hints must span first..last intron base)', file=sys.stderr)
    if chains is not None:
        cs = defaultdict(int)
        for s in S:
            cs[s['chain_support']] += 1
        hi = [s for s in S if s['tier_base'] >= 9]
        n_full = sum(s['chain_support'] == 'full' for s in hi)
        print(f'full-length transcripts: {len(chains[0])}; models ' +
              ', '.join(f'{k}={v}' for k, v in sorted(cs.items())) +
              f'; tier>=9 models with the chain in one transcript: {n_full / max(len(hi), 1):.1%}', file=sys.stderr)
    print('tier adjustments: ' + ', '.join(f'{k}={v}' for k, v in sorted(adj.items())), file=sys.stderr)
    keys = (['primary', 'isoform_classifier', 'isoform_rejected', 'locus_below_min', 'dropped_lone_single_exon']
            if clf is not None else
            ['primary', 'isoform_tier', 'isoform_chain', 'isoform_drusilla_td2', 'isoform_rejected', 'locus_below_min'])
    for k in keys:
        print(f'  {k}: {why[k]}', file=sys.stderr)


if __name__ == '__main__':
    main()