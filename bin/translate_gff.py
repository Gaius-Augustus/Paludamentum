#!/usr/bin/env python3
"""Proteins of the coding transcripts of a GFF3 or GTF file, for a genome
whose genetic code is not the standard one.

gffread -y translates with the standard code only. This script writes what
`gffread -g genome -y out -S` writes (one record per transcript with CDS,
named by the transcript ID, the stop codon at the end left out, internal
stop codons as *), with the codons read in the NCBI translation table given.

Transcript of a CDS line: Parent (GFF3), else transcript_id (GTF). The phase
of the first CDS segment in transcription order sets the reading frame.

Usage: translate_gff.py --gff in.gff3 --genome genome.fa --translation-table 6 --out proteins.fa
"""
import argparse
import os
import re
import sys
from collections import OrderedDict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gff3_lib import GFF3Error, codon_table, read_fasta, subseq, translate  # noqa: E402


def transcript_of(attributes: str):
    """Transcript IDs of a CDS line (a GFF3 CDS may have several parents)."""
    m = re.search(r'transcript_id "([^"]+)"', attributes)
    if m:
        return [m.group(1)]
    m = re.search(r"(?:^|;)\s*Parent=([^;]+)", attributes)
    return m.group(1).strip().split(",") if m else []


def read_cds(path: str):
    """{transcript: (seqid, strand, [(start, end, phase)])} in file order."""
    cds = OrderedDict()
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            cols = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(cols) < 9 or cols[2] != "CDS":
                continue
            for tid in transcript_of(cols[8]):
                entry = cds.setdefault(tid, (cols[0], cols[6], []))
                entry[2].append((int(cols[3]), int(cols[4]), cols[7]))
    return cds


def protein(seqid: str, strand: str, parts, genome, table: int) -> str:
    parts = sorted(parts, reverse=strand == "-")
    seq = "".join(subseq(genome, seqid, start, end, strand) for start, end, _ in parts).upper()
    phase = parts[0][2]
    aa = translate(seq[int(phase) if phase in ("1", "2") else 0:], table)
    return aa[:-1] if aa.endswith("*") else aa


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff", required=True, help="GFF3 or GTF with CDS lines")
    parser.add_argument("--genome", required=True, help="genome FASTA")
    parser.add_argument("--translation-table", type=int, required=True, help="NCBI translation table")
    parser.add_argument("--out", required=True, help="protein FASTA")
    args = parser.parse_args(argv)
    try:
        codon_table(args.translation_table)
        cds = read_cds(args.gff)
        genome = read_fasta(args.genome, {entry[0] for entry in cds.values()})
        with open(args.out, "w", encoding="utf-8") as out:
            for tid, (seqid, strand, parts) in cds.items():
                if seqid not in genome:
                    raise GFF3Error(f"sequence {seqid} of transcript {tid} is not in the genome")
                aa = protein(seqid, strand, parts, genome, args.translation_table)
                if aa:
                    out.write(f">{tid}\n")
                    for i in range(0, len(aa), 70):
                        out.write(aa[i:i + 70] + "\n")
    except GFF3Error as exc:
        print(f"translate_gff.py: {args.gff}: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
