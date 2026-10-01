#!/usr/bin/env python3
"""High-confidence genes: the TransDecoder ORFs that pass the strict DIAMOND
filter, one ORF per StringTie transcript, in genome coordinates (training.gff).

The Tiberius pipeline also ran an intrinsic stage on this protein-supported
set and appended its result to the same file, so it never changed
training.gff; it is not part of Paludamentum. choose_one_isoform keeps the
longest ORF per StringTie transcript, not per gene, so the alternative
isoforms reach the final merge.
"""
import argparse
import logging

from hc_module import (
    getting_hc_supported_by_proteins,
    choose_one_isoform, from_pep_file_to_gff3, from_transcript_to_genome)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--diamond_revised", required=True)
    ap.add_argument("--revised_pep", required=True)
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--stringtie_gtf", required=True)
    ap.add_argument("--stringtie_gff3", required=True)
    ap.add_argument("--transcripts_fasta", required=True)
    ap.add_argument("--transdecoder_util", required=True,
                    help="cdna_alignment_orf_to_genome_orf.pl of TransDecoder (command or path)")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    hc_pep = getting_hc_supported_by_proteins(args.diamond_revised, args.revised_pep, args.proteins,
                                              f"{args.outdir}/hc_genes.pep")
    hc_single_pep = choose_one_isoform(hc_pep, f"{args.outdir}/hc_one_isoform.pep")
    hc_single_gff = from_pep_file_to_gff3(hc_single_pep, args.stringtie_gtf, f"{args.outdir}/hc_one_isoform.gff3")
    from_transcript_to_genome(hc_single_gff, args.stringtie_gff3, args.transcripts_fasta,
                              f"{args.outdir}/training.gff", args.transdecoder_util)


if __name__ == "__main__":
    main()
