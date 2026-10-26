#!/usr/bin/env python3
"""Convert the revised TD2 peptides (CDS_CLASSIFY_AND_REVISE output) back to a
genome-coordinate GFF3 that combine_genesets.py can read as --td2.

Two steps:
 1. revised_candidates.pep + stringtie.gtf -> transcript-coordinate GFF3 of
    the CDS of each ORF (and its UTRs on the transcript).
 2. cdna_alignment_orf_to_genome_orf.pl maps those CDS blocks onto the genome
    using stringtie.gff3 and the transcript FASTA.

Keeps every isoform; combine_genesets.py does its own ranking. Logic is a
trimmed copy of the orf-to-genome part of bin/hc_module.py (from Tiberius).
"""
import argparse
import logging
import re
import subprocess
import sys

from Bio import SeqIO


ORF_COORDS = re.compile(r":(\d+)-(\d+)\((\+|-)\)")
GTF_TID = re.compile(r'transcript_id "([^"]+)"')


def pep_to_transcript_gff3(orf_pep, stringtie_gtf, out_gff3):
    transcript_length = {}
    with open(stringtie_gtf) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'exon':
                continue
            m = GTF_TID.search(f[8])
            if not m:
                continue
            transcript_length[m.group(1)] = transcript_length.get(m.group(1), 0) + (int(f[4]) - int(f[3]) + 1)

    stringtie_src, td_src = 'StringTie', 'TransDecoder'
    with open(out_gff3, 'w') as out:
        for record in SeqIO.parse(orf_pep, 'fasta'):
            tid = record.id.split('.p')[0]
            tlen = transcript_length.get(tid)
            if tlen is None:
                logging.warning('transcript %s not in GTF, skipping', tid)
                continue
            m = ORF_COORDS.search(record.description)
            if not m:
                logging.warning('no ORF coordinates on %s header, skipping', record.id)
                continue
            start, stop, strand = int(m.group(1)), int(m.group(2)), m.group(3)
            parts = record.description.split()
            gene_id = f'{parts[0]}.gene'
            gene_name = parts[1] if len(parts) > 1 else parts[0]
            out.write(f'{tid}\t{stringtie_src}\tgene\t1\t{tlen}\t.\t{strand}\t.\tID={gene_id};Name={gene_name}\n')
            out.write(f'{tid}\t{stringtie_src}\tmRNA\t1\t{tlen}\t.\t{strand}\t.\tID={record.id};Parent={gene_id};Name={gene_name}\n')
            out.write(f'{tid}\t{stringtie_src}\texon\t1\t{tlen}\t.\t{strand}\t.\tID={record.id}.exon1;Parent={record.id}\n')
            out.write(f'{tid}\t{td_src}\tCDS\t{start}\t{stop}\t.\t{strand}\t0\tID=cds.{record.id};Parent={record.id}\n')
            if strand == '+':
                if start > 1:
                    out.write(f'{tid}\t{td_src}\tfive_prime_UTR\t1\t{start - 1}\t.\t{strand}\t.\tID={record.id}.utr5p1;Parent={record.id}\n')
                if stop < tlen:
                    out.write(f'{tid}\t{td_src}\tthree_prime_UTR\t{stop + 1}\t{tlen}\t.\t{strand}\t.\tID={record.id}.utr3p1;Parent={record.id}\n')
            else:
                if stop < tlen:
                    out.write(f'{tid}\t{td_src}\tfive_prime_UTR\t{stop + 1}\t{tlen}\t.\t{strand}\t.\tID={record.id}.utr5p1;Parent={record.id}\n')
                if start > 1:
                    out.write(f'{tid}\t{td_src}\tthree_prime_UTR\t1\t{start - 1}\t.\t{strand}\t.\tID={record.id}.utr3p1;Parent={record.id}\n')
            out.write('\n')


def transcript_to_genome(cds_gff3, stringtie_gff3, transcripts_fasta, out_gff3, transdecoder_util):
    cmd = [transdecoder_util, cds_gff3, stringtie_gff3, transcripts_fasta]
    with open(out_gff3, 'w') as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True, check=False)
    if r.returncode != 0:
        sys.stderr.write(r.stderr)
        sys.exit(f'cdna_alignment_orf_to_genome_orf.pl failed with code {r.returncode}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--revised-pep', required=True, help='revised_candidates.pep (CDS_CLASSIFY_AND_REVISE output)')
    ap.add_argument('--stringtie-gtf', required=True)
    ap.add_argument('--stringtie-gff3', required=True)
    ap.add_argument('--transcripts-fasta', required=True, help='transcripts.fasta (cdna of TD_ALL)')
    ap.add_argument('--transdecoder-util', required=True,
                    help='cdna_alignment_orf_to_genome_orf.pl command or path')
    ap.add_argument('-o', '--out', required=True, help='genome-coordinate GFF3')
    a = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    tmp = a.out + '.transcript.gff3'
    pep_to_transcript_gff3(a.revised_pep, a.stringtie_gtf, tmp)
    transcript_to_genome(tmp, a.stringtie_gff3, a.transcripts_fasta, a.out, a.transdecoder_util)


if __name__ == '__main__':
    main()
