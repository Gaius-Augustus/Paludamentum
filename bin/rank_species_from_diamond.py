#!/usr/bin/env python3
import re
import sys
from collections import defaultdict

if len(sys.argv) < 2:
    sys.stderr.write(f"Usage: {sys.argv[0]} diamond_hits.tsv [TOP_N]\n")
    sys.exit(1)

in_path = sys.argv[1]
TOP_N = int(sys.argv[2]) if len(sys.argv) > 2 else 8

# thresholds (tune if you want)
MAX_EVALUE = 1e-5
MIN_PIDENT = 30.0
MIN_QCOV   = 0.5   # coverage = aligned_length / qlen

ORTHODB_ID = re.compile(r"^\d+_\d+:")   # OrthoDB v10-v12 ids: <taxid>_<n>:<hex>


def get_species_id(sseqid: str) -> str:
    """
    Species of a subject: the NCBI taxonomy id in front of the first '_' of an
    OrthoDB-style id ('101020_0:000003' -> '101020'). For any other id format
    the whole id is returned, so that every protein counts as its own
    "species" and nothing is grouped by a meaningless prefix.
    """
    head = sseqid.split()[0]
    if ORTHODB_ID.match(head):
        return head.split("_", 1)[0]
    return head

# best hit *per query* (over all species)
# q -> (best_score, best_species)
best_hit_for_query = {}

with open(in_path, "r") as fh:
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 8:
            continue

        qseqid, sseqid = cols[0], cols[1]
        pident = float(cols[2])
        alen   = float(cols[3])
        evalue = float(cols[4])
        bits   = float(cols[5])
        qlen   = float(cols[6])
        slen   = float(cols[7])

        if evalue > MAX_EVALUE:
            continue
        if pident < MIN_PIDENT:
            continue
        if qlen <= 0:
            continue

        qcov = alen / qlen
        if qcov < MIN_QCOV:
            continue

        species = get_species_id(sseqid)
        score   = bits * qcov   # combined strength

        if qseqid not in best_hit_for_query or score > best_hit_for_query[qseqid][0]:
            best_hit_for_query[qseqid] = (score, species)

# aggregate per species from the best hits
species_votes = defaultdict(int)
species_score = defaultdict(float)

for q, (score, sp) in best_hit_for_query.items():
    species_votes[sp] += 1
    species_score[sp] += score

ranking = []
for sp in species_votes:
    ranking.append((sp, species_votes[sp], species_score[sp]))

# sort: 1) #queries, 2) total score
ranking.sort(key=lambda x: (-x[1], -x[2]))

with open("top_species.txt", "w") as out_sp:
    for i, (sp, n_q, tot) in enumerate(ranking):
        mark = "TOP" if i < TOP_N else ""
        sys.stdout.write(f"{i+1}\t{sp}\tqueries={n_q}\ttotal_score={tot:.1f}\t{mark}\n")
        if i < TOP_N:
            out_sp.write(sp + "\n")

if not ranking:
    sys.stderr.write("No DIAMOND hit passed the thresholds: top_species.txt is empty. "
                     "The caller must then keep the whole protein database.\n")
else:
    sys.stderr.write(f"Wrote top {min(TOP_N, len(ranking))} species to top_species.txt\n")
