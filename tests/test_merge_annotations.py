"""bin/merge_annotations.py: GFF3 feature types as in NCBI/Ensembl files (Annotrieve)."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "bin" / "merge_annotations.py"

CODING_GTF = (
    'chr1\tTiberius\texon\t100\t300\t.\t+\t.\tgene_id "g1"; transcript_id "g1.t1";\n'
    'chr1\tTiberius\tCDS\t100\t300\t.\t+\t0\tgene_id "g1"; transcript_id "g1.t1";\n'
)
NONCODING_GTF = (
    'chr2\tStringTie\texon\t500\t900\t.\t-\t.\tgene_id "s1"; transcript_id "s1.t1";\n'
)


def merge(tmp_path: Path, *contents: str) -> list[list[str]]:
    inputs = []
    for i, text in enumerate(contents):
        path = tmp_path / f"in{i}.gtf"
        path.write_text(text)
        inputs.append(str(path))
    result = subprocess.run([sys.executable, str(SCRIPT), "--mode", "full", *inputs],
                            capture_output=True, text=True, check=True)
    return [line.split("\t") for line in result.stdout.splitlines()
            if line and not line.startswith("#")]


def features(rows: list[list[str]], seqid: str) -> dict[str, str]:
    return {row[2]: row[8] for row in rows if row[0] == seqid}


def test_coding_transcript_is_mrna_and_gene_is_protein_coding(tmp_path):
    coding = features(merge(tmp_path, CODING_GTF, NONCODING_GTF), "chr1")
    assert "transcript" not in coding
    assert coding["mRNA"] == "ID=gene_000001.t1;Parent=gene_000001"
    assert coding["gene"] == "ID=gene_000001;gene_biotype=protein_coding"


def test_noncoding_transcript_keeps_type_and_has_no_biotype(tmp_path):
    noncoding = features(merge(tmp_path, CODING_GTF, NONCODING_GTF), "chr2")
    assert "mRNA" not in noncoding
    assert noncoding["transcript"] == "ID=gene_000002.t1;Parent=gene_000002"
    assert noncoding["gene"] == "ID=gene_000002"


def test_merged_output_can_be_merged_again(tmp_path):
    """The ab initio GFF3 is merged again with the HC genes (MERGE_GENEFINDER_TRAIN)."""
    first = tmp_path / "first.gff3"
    first.write_text("\n".join("\t".join(r) for r in merge(tmp_path, CODING_GTF)) + "\n")
    again = features(merge(tmp_path, first.read_text()), "chr1")
    assert again["gene"] == "ID=gene_000001;gene_biotype=protein_coding"
    assert "mRNA" in again and "CDS" in again


# A transcript as Tiberius and Vipsania write it themselves (bricks2marble):
# GTF with gene, transcript, codon and intron lines, GFF3 without a gene line.
NATIVE_GTF = (
    'chr1\tTiberius\tgene\t100\t900\t.\t+\t.\tgene_id "g1";\n'
    'chr1\tTiberius\ttranscript\t100\t900\t.\t+\t.\tgene_id "g1"; transcript_id "tx1";\n'
    'chr1\tTiberius\tstart_codon\t100\t102\t.\t+\t0\tgene_id "g1"; transcript_id "tx1";\n'
    'chr1\tTiberius\tstop_codon\t898\t900\t.\t+\t0\tgene_id "g1"; transcript_id "tx1";\n'
    'chr1\tTiberius\tCDS\t100\t300\t.\t+\t0\tgene_id "g1"; transcript_id "tx1";\n'
    'chr1\tTiberius\texon\t100\t300\t.\t+\t0\tgene_id "g1"; transcript_id "tx1";\n'
    'chr1\tTiberius\tintron\t301\t499\t.\t+\t0\tgene_id "g1"; transcript_id "tx1";\n'
    'chr1\tTiberius\tCDS\t500\t900\t.\t+\t0\tgene_id "g1"; transcript_id "tx1";\n'
    'chr1\tTiberius\texon\t500\t900\t.\t+\t0\tgene_id "g1"; transcript_id "tx1";\n'
)
NATIVE_GFF3 = (
    '##gff-version 3\n'
    'chr1\tTiberius\tmRNA\t100\t900\t.\t+\t.\tID=tx1;Name=tx1\n'
    'chr1\tTiberius\tCDS\t100\t300\t.\t+\t0\tID=tx1.CDS.1;Parent=tx1\n'
    'chr1\tTiberius\texon\t100\t300\t.\t+\t0\tID=tx1.exon.1;Parent=tx1\n'
    'chr1\tTiberius\tCDS\t500\t900\t.\t+\t0\tID=tx1.CDS.2;Parent=tx1\n'
    'chr1\tTiberius\texon\t500\t900\t.\t+\t0\tID=tx1.exon.2;Parent=tx1\n'
)


def test_gene_finder_gff3_gives_the_same_annotation_as_its_gtf(tmp_path):
    """An existing prediction (tiberius.result, vipsania.result) may be the GFF3 of the gene finder."""
    from_gff3 = merge(tmp_path, NATIVE_GFF3)
    assert from_gff3 == merge(tmp_path, NATIVE_GTF)
    assert [row[2] for row in from_gff3] == ["gene", "mRNA", "exon", "exon", "CDS", "CDS"]
