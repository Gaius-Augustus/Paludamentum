# From Tiberius (aae4a8a), tiberius/scripts/hc_module.py.
# Copyright (c) 2023 Lars Gabriel. MIT License, see LICENSE-Tiberius.
# Changed in Paludamentum: the intrinsic HC stage and the alignment helpers
# are removed (the pipeline never used their result), gene and mRNA lines of
# the transcript GFF3 get distinct IDs, and the TransDecoder utility is the
# executable that the caller passes.
import logging, os, re
from Bio import SeqIO
from utility import run_subprocess


#Authors: "Amrei Knuth", "Lars Gabriel"
#Credits: "Katharina Hoff"
#Email:"lars.gabriel@uni-greifswald.de"
#Date: "Janurary 2025"

logger = logging.getLogger(__name__)


def getting_hc_supported_by_proteins(
    diamond_tsv,
    transdecoder_pep,
    protein_file,
    output_path="hc_genes_strict.pep",
    min_pident=50.0,    # stricter identity
    min_qcov=0.9,       # >=90% of ORF covered
    min_tcov=0.9,       # >=90% of DB protein covered
    max_evalue=1e-5,   # very strong hits only
    require_complete=True,
):
    """
    Select a *very high-precision* set of genes supported by protein evidence.

    Assumes DIAMOND was run with default outfmt 6:
      qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

    This is intentionally strict and will have low recall.
    """

    logging.info("Selecting VERY STRICT high-confidence genes supported by protein evidence...")

    out_dir = os.path.dirname(output_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Target protein lengths (from DB)
    t_length_dict = {
        record.id: len(record.seq)
        for record in SeqIO.parse(protein_file, "fasta")
    }

    # Query ORFs
    q_dict = {
        record.id: record
        for record in SeqIO.parse(transdecoder_pep, "fasta")
    }
    logging.info("Loaded %d TransDecoder ORFs", len(q_dict))

    seen_hc_queries = set()

    with open(diamond_tsv) as tsv, open(output_path, "w") as output:
        for line in tsv:
            if not line.strip() or line.startswith("#"):
                continue

            part = line.rstrip("\n").split("\t")
            if len(part) < 12:
                logging.warning("Skipping malformed DIAMOND line: %s", line.strip())
                continue

            query_id = part[0]
            # DIAMOND escapes embedded whitespace in subject ids as the
            # literal characters '\t' / '\n' / '\\'; collapse to the
            # accession (the part before the first such escape or real
            # whitespace) so it matches Biopython's record.id from the
            # FASTA, which splits at the first whitespace.
            protein_id = re.split(r"\\[tn\\]|\s", part[1], maxsplit=1)[0]
            aaident = float(part[2])
            align_length = int(part[3])
            q_start = int(part[6])
            q_end = int(part[7])
            t_start = int(part[8])
            t_end = int(part[9])
            evalue = float(part[10])
            bitscore = float(part[11])

            # Skip queries that are not from TransDecoder
            if query_id not in q_dict:
                continue

            # Avoid writing same query multiple times
            if query_id in seen_hc_queries:
                continue

            record = q_dict[query_id]

            # Only complete ORFs (optional)
            if require_complete and "type:complete" not in record.description:
                continue

            # Query length without trailing '*'
            seq_str = str(record.seq)
            if seq_str.endswith("*"):
                q_length = len(seq_str) - 1
            else:
                q_length = len(seq_str)

            if q_length <= 0:
                continue

            t_length = t_length_dict.get(protein_id, 0)
            if t_length <= 0:
                # No length info -> skip for strict mode
                continue

            # Coverage
            q_cov = align_length / q_length
            t_cov = align_length / t_length

            # OPTIONAL: also require near-N-terminal alignment
            # This is VERY strict; uncomment if you want it:
            # if abs(q_start - 1) > 5 or abs(t_start - 1) > 5:
            #     continue

            # Strict criteria
            if (
                evalue <= max_evalue
                and aaident >= min_pident
                and q_cov >= min_qcov
                and t_cov >= min_tcov
            ):
                SeqIO.write(record, output, "fasta")
                seen_hc_queries.add(query_id)
                # we don't break here because each line is a different query or skipped by seen_hc_queries

    logging.info(
        "VERY STRICT high-confidence selection completed: %d ORFs written to %s",
        len(seen_hc_queries),
        output_path,
    )
    return output_path

def choose_one_isoform(input_pep, output_pep):
    """
    Selects one isoform per gene based on length.

    If multiple isoforms exist, the longest is selected. If ties occur, the first 
    encountered isoform is chosen.

    :param input_pep: Path to the input PEP file.
    :param output_pep: Path to the output PEP file with one isoform per gene.
    """
    logging.info("Selecting one isoform per gene...")

    isoform_dict = {}

    # Parse sequences and group by gene ID
    for record in SeqIO.parse(input_pep, "fasta"):
        gene_id = record.id.split(".p")[0]
        cds_coords = re.search(r":(\d+)-(\d+)\((\+|-)\)", record.description)
        start, stop = int(cds_coords.group(1)), int(cds_coords.group(2))

        if gene_id not in isoform_dict:
            isoform_dict[gene_id] = [(start, stop, record)]
        else:
            isoform_dict[gene_id].append((start, stop, record))

    # Select and write longest isoform per gene
    with open(output_pep, "w") as output:
        for gene_id, isoforms in isoform_dict.items():
            record = max(isoforms, key=lambda x: x[1] - x[0])[2]  # Select longest
            SeqIO.write(record, output, "fasta")

    logging.info("Isoform selection complete. Output saved to %s", output_pep)

    return output_pep


def from_transcript_to_genome(cds_gff3, transcripts_gff3, transcripts_fasta, output_path, 
                        transdecoder_util):
    """
    Convert transcript coordinates to genome coordinates using TransDecoder.

    :param cds_gff3: Path to CDS annotation in GFF3 format in transcript coordinates.
    :param transcripts_gff3: Path to transcript annotation in GFF3 format in genome coordinates.
    :param transcripts_fasta: Path to the transcript sequences in FASTA format.
    :param output_path: Path to the output genome-coordinate GFF3 file.
    :param transdecoder_util: cdna_alignment_orf_to_genome_orf.pl of TransDecoder:
        a command on the PATH or a path (params.tools.transdecoder_util_orf2genome).
    """
    logging.info("Converting transcript coordinates to genome coordinates...")

    command = [
        transdecoder_util,
        cds_gff3,
        transcripts_gff3,
        transcripts_fasta
    ]

    with open(output_path, "w") as output:
        run_subprocess(command, stdout=output, capture_output=False,
                error_message="TransDecoder coordinate transformation failed.")

    logging.info("Transcript-to-genome coordinate transformation complete.")

def from_pep_file_to_gff3(orf_pep, transcript_gtf, output_gff3):
    """
    Convert a TransDecoder peptide FASTA file and StringTie transcript annotations (GTF) into a GFF3 file.

    This function calculates transcript lengths from the GTF file and generates a GFF3 file containing 
    mRNA, exon, CDS, and UTR annotations.

    :param orf_pep: Path to the TransDecoder peptide FASTA file.
    :param transcript_gtf: Path to the StringTie transcript annotation file (GTF).
    :param output_gff3: Path to the output GFF3 file.
    :return: Path to the generated GFF3 file.
    """
    logging.info("Generating GFF3 file from TransDecoder PEP and StringTie GTF...")

    # Parse transcript lengths from GTF
    transcript_lengths = {}

    with open(transcript_gtf, "r") as transcript_file:
        for line in transcript_file:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if parts[2] == "exon":
                start, stop = int(parts[3]), int(parts[4])
                length = stop - start + 1
                transcript_id_match = re.search(r'transcript_id "([^"]+)"', parts[8])

                if transcript_id_match:
                    transcript_id = transcript_id_match.group(1)
                    transcript_lengths[transcript_id] = transcript_lengths.get(transcript_id, 0) + length

    # Create the GFF3 file
    tool_stringtie = "StringTie"
    tool_transdecoder = "TransDecoder"

    with open(output_gff3, "w") as output:
        for record in SeqIO.parse(orf_pep, "fasta"):
            transcript_id = record.id.split(".p")[0]

            if transcript_id not in transcript_lengths:
                logging.warning(f"Transcript ID {transcript_id} not found in GTF, skipping.")
                continue

            transcript_length = transcript_lengths[transcript_id]
            description = record.description

            # Extract gene ID, gene name, ORF coordinates, and strand
            # gene_id_match = re.search(r"gene=([^ ]+)", description)
            # gene_name_match = re.search(r"Name=([^ ]+)", description)
            orf_coords_match = re.search(r":(\d+)-(\d+)\((\+|-)\)", description)

            if not orf_coords_match:
                logging.warning(f"Skipping {record.id} due to missing annotation.")
                continue

            description_parts = description.split()
            # The gene needs an ID of its own: with the ORF ID for both, the
            # mRNA would be its own Parent, which is invalid GFF3.
            gene_id = f"{description_parts[0]}.gene"
            gene_name = description_parts[1]
            orf_start, orf_stop, strand = int(orf_coords_match.group(1)), int(orf_coords_match.group(2)), orf_coords_match.group(3)

            # Write GFF3 records
            output.write(f"{transcript_id}\t{tool_stringtie}\tgene\t1\t{transcript_length}\t.\t{strand}\t.\tID={gene_id};Name={gene_name}\n")
            output.write(f"{transcript_id}\t{tool_stringtie}\tmRNA\t1\t{transcript_length}\t.\t{strand}\t.\tID={record.id};Parent={gene_id};Name={gene_name}\n")
            output.write(f"{transcript_id}\t{tool_stringtie}\texon\t1\t{transcript_length}\t.\t{strand}\t.\tID={record.id}.exon1;Parent={record.id}\n")
            output.write(f"{transcript_id}\t{tool_transdecoder}\tCDS\t{orf_start}\t{orf_stop}\t.\t{strand}\t0\tID=cds.{record.id};Parent={record.id}\n")

            # Add UTR annotations if applicable
            if strand == "+":
                if orf_start > 1:
                    output.write(f"{transcript_id}\t{tool_transdecoder}\tfive_prime_UTR\t1\t{orf_start-1}\t.\t{strand}\t.\tID={record.id}.utr5p1;Parent={record.id}\n")
                if orf_stop < transcript_length:
                    output.write(f"{transcript_id}\t{tool_transdecoder}\tthree_prime_UTR\t{orf_stop+1}\t{transcript_length}\t.\t{strand}\t.\tID={record.id}.utr3p1;Parent={record.id}\n")
            elif strand == "-":
                if orf_stop < transcript_length:
                    output.write(f"{transcript_id}\t{tool_transdecoder}\tfive_prime_UTR\t{orf_stop+1}\t{transcript_length}\t.\t{strand}\t.\tID={record.id}.utr5p1;Parent={record.id}\n")
                if orf_start > 1:
                    output.write(f"{transcript_id}\t{tool_transdecoder}\tthree_prime_UTR\t1\t{orf_start-1}\t.\t{strand}\t.\tID={record.id}.utr3p1;Parent={record.id}\n")

            output.write("\n")

    logging.info("GFF3 file created successfully: %s", output_gff3)
    return output_gff3
