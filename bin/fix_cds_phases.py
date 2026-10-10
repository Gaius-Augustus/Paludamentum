#!/usr/bin/env python3
"""
Recompute the CDS phases of the transcripts of a GTF that have a start codon.

The Tiberius of the hint rescue image (gaiusaugustus/paludamentum-hint-rescue)
writes the phases of minus-strand transcripts counted from the wrong end of
the CDS: the 5' CDS that holds the ATG gets phase 1 or 2 instead of 0, and the
phases after it follow suit (T. rubripes, 2026-10-07: 1,105 of 1,182 rescued
minus-strand transcripts; plus strand and the Tiberius of the tiberius image
are correct). The structures are right; only the phase column is not, and the
sanity filter reads a 5' phase other than 0 as a 5'-partial CDS whose length
is then no multiple of 3.

A transcript with a start_codon line has a complete 5' end, so its CDS starts
at phase 0 in transcription order and every following phase is determined by
the lengths before it. Those phases are rewritten; transcripts without a
start_codon line, and every line that is not a CDS line, pass unchanged.

Usage: fix_cds_phases.py in.gtf > out.gtf   (stderr: number of transcripts changed)
"""

import re
import sys
from collections import defaultdict

TX_ID = re.compile(r'transcript_id "([^"]+)"')


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__.strip().splitlines()[-1])
    lines = open(sys.argv[1]).read().splitlines(keepends=True)

    cds = defaultdict(list)          # transcript -> [(index of line, start, end, strand)]
    complete = set()
    for i, line in enumerate(lines):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) < 9:
            continue
        m = TX_ID.search(f[8])
        if not m:
            continue
        if f[2] == "start_codon":
            complete.add(m.group(1))
        elif f[2] == "CDS":
            cds[m.group(1)].append((i, int(f[3]), int(f[4]), f[6]))

    changed = 0
    for tx in complete & cds.keys():
        segments = sorted(cds[tx], key=lambda s: s[1], reverse=cds[tx][0][3] == "-")
        done = 0
        before = False
        for i, start, end, _ in segments:
            f = lines[i].rstrip("\n").split("\t")
            phase = str((3 - done % 3) % 3)
            if f[7] != phase:
                f[7] = phase
                lines[i] = "\t".join(f) + "\n"
                before = True
            done += end - start + 1
        changed += before

    sys.stdout.writelines(lines)
    print(f"fix_cds_phases.py: phases rewritten in {changed} of {len(complete & cds.keys())} "
          f"transcripts with a start codon", file=sys.stderr)


if __name__ == "__main__":
    main()
