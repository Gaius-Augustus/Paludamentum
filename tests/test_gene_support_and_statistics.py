"""Tests of bin/gene_support.py and bin/gene_set_statistics.py on toy GFF3 and hint files."""
from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "bin"


def load(name: str):
    spec = importlib.util.spec_from_file_location(name, BIN / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run(script: str, *args: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(BIN / script), *args], capture_output=True, text=True, cwd=cwd)


def gff(*lines: str) -> str:
    return "##gff-version 3\n" + "".join("\t".join(line.split()) + "\n" for line in lines)


def support_rows(path: Path) -> dict:
    lines = [line for line in path.read_text().splitlines() if line and not line.startswith("#")]
    header = lines[0].split("\t")
    return {row["transcript_id"]: row for row in (dict(zip(header, line.split("\t"))) for line in lines[1:])}


# ---------------------------------------------------------------- gene_support.py

GENES = gff(
    # g1: + strand, 3 CDS segments, introns 201-300 (RNA-Seq) and 401-500 (RNA-Seq and protein)
    "chr1 T gene 101 600 . + . ID=g1;gene_biotype=protein_coding",
    "chr1 T mRNA 101 600 . + . ID=g1.t1;Parent=g1",
    "chr1 T exon 101 200 . + . ID=g1.t1.exon1;Parent=g1.t1",
    "chr1 T exon 301 400 . + . ID=g1.t1.exon2;Parent=g1.t1",
    "chr1 T exon 501 600 . + . ID=g1.t1.exon3;Parent=g1.t1",
    "chr1 T CDS 101 200 . + 0 ID=g1.t1.cds;Parent=g1.t1",
    "chr1 T CDS 301 400 . + 2 ID=g1.t1.cds;Parent=g1.t1",
    "chr1 T CDS 501 600 . + 1 ID=g1.t1.cds;Parent=g1.t1",
    # g2: - strand, the hints of its intron and CDS are on + (strand mismatch)
    "chr1 T gene 1001 1300 . - . ID=g2;gene_biotype=protein_coding",
    "chr1 T mRNA 1001 1300 . - . ID=g2.t1;Parent=g2",
    "chr1 T exon 1001 1100 . - . ID=g2.t1.exon1;Parent=g2.t1",
    "chr1 T exon 1201 1300 . - . ID=g2.t1.exon2;Parent=g2.t1",
    "chr1 T CDS 1001 1100 . - 2 ID=g2.t1.cds;Parent=g2.t1",
    "chr1 T CDS 1201 1300 . - 0 ID=g2.t1.cds;Parent=g2.t1",
    # g3: intron 201-300 by protein only; CDS 1 overlapped by a CDSpart hint, CDS 2 by an
    # RNA-Seq exonpart hint without strand
    "chr2 T gene 101 400 . + . ID=g3;gene_biotype=protein_coding",
    "chr2 T mRNA 101 400 . + . ID=g3.t1;Parent=g3",
    "chr2 T exon 101 200 . + . ID=g3.t1.exon1;Parent=g3.t1",
    "chr2 T exon 301 400 . + . ID=g3.t1.exon2;Parent=g3.t1",
    "chr2 T CDS 101 200 . + 0 ID=g3.t1.cds;Parent=g3.t1",
    "chr2 T CDS 301 400 . + 2 ID=g3.t1.cds;Parent=g3.t1",
    # g4: single-exon
    "chr2 T gene 1001 1300 . + . ID=g4;gene_biotype=protein_coding",
    "chr2 T mRNA 1001 1300 . + . ID=g4.t1;Parent=g4",
    "chr2 T exon 1001 1300 . + . ID=g4.t1.exon1;Parent=g4.t1",
    "chr2 T CDS 1001 1300 . + 0 ID=g4.t1.cds;Parent=g4.t1",
    # g5: transcript without CDS
    "chr2 T gene 2001 2300 . + . ID=g5",
    "chr2 T transcript 2001 2300 . + . ID=g5.t1;Parent=g5",
    "chr2 T exon 2001 2300 . + . ID=g5.t1.exon1;Parent=g5.t1",
)

HINTS = gff(
    "chr1 b2h intron 201 300 5 + . mult=5;pri=4;src=E",
    "chr1 b2h intron 401 500 2 + . mult=2;pri=4;src=E",
    "chr1 miniprot2h intron 401 500 0 + . src=P;grp=p1;pri=4",
    "chr1 miniprot2h intron 1101 1200 0 + . src=P;grp=p2;pri=4",
    "chr1 miniprot2h CDSpart 1020 1080 0 + 0 src=P;grp=p2;pri=4",
    "chr2 miniprot2h intron 201 300 0 + . src=P;grp=p3;pri=4",
    "chr2 miniprot2h CDSpart 120 180 0 + 0 src=P;grp=p3;pri=4",
    "chr2 x exonpart 350 360 1 . . src=E;pri=4",
    "chr2 x intron 1101 1200 1 + . src=X",                   # other source, ignored
    "chr2 x CDSpart 1001 1300 1 + . src=M",                  # other source, ignored
)


@pytest.fixture
def support(tmp_path: Path):
    (tmp_path / "genes.gff3").write_text(GENES)
    (tmp_path / "hints.gff").write_text(HINTS)
    proc = run("gene_support.py", "--gff3", "genes.gff3", "--hints", "hints.gff", "--out", "support.tsv",
               cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    return proc, support_rows(tmp_path / "support.tsv"), (tmp_path / "support.tsv").read_text()


def test_columns_and_skipped_noncoding(support):
    _proc, rows, text = support
    header = [line for line in text.splitlines() if not line.startswith("#")][0].split("\t")
    assert header == load("gene_support").COLUMNS
    assert header[:5] == ["gene_id", "transcript_id", "chrom", "strand", "num_introns"]
    assert list(rows) == ["g1.t1", "g2.t1", "g3.t1", "g4.t1"]        # g5.t1 has no CDS
    assert rows["g1.t1"]["gene_id"] == "g1"


def test_intron_support_rnaseq_and_both(support):
    _proc, rows, _text = support
    r = rows["g1.t1"]
    assert r["num_introns"] == "2"
    assert (r["introns_sup_rnaseq"], r["introns_sup_protein"], r["introns_sup_any"]) == ("2", "1", "2")
    assert (r["pct_introns_sup_rnaseq"], r["pct_introns_sup_protein"], r["pct_introns_sup_any"]) == \
        ("100.0", "50.0", "100.0")
    assert r["num_exons"] == "3" and r["exons_sup_any"] == "0"


def test_intron_support_protein_and_cds_overlap(support):
    _proc, rows, _text = support
    r = rows["g3.t1"]
    assert (r["introns_sup_rnaseq"], r["introns_sup_protein"], r["introns_sup_any"]) == ("0", "1", "1")
    # CDS 101-200 overlapped by CDSpart (protein), CDS 301-400 by an unstranded exonpart (RNA-Seq)
    assert (r["exons_sup_rnaseq"], r["exons_sup_protein"], r["exons_sup_any"]) == ("1", "1", "2")
    assert r["pct_exons_sup_any"] == "100.0"


def test_strand_mismatch_not_counted(support):
    _proc, rows, _text = support
    r = rows["g2.t1"]
    assert r["strand"] == "-"
    assert (r["num_introns"], r["introns_sup_any"], r["pct_introns_sup_any"]) == ("1", "0", "0.0")
    assert r["exons_sup_any"] == "0"


def test_single_exon_and_other_sources(support):
    proc, rows, text = support
    r = rows["g4.t1"]
    assert (r["num_introns"], r["introns_sup_any"], r["pct_introns_sup_any"]) == ("0", "0", "NA")
    assert (r["num_exons"], r["exons_sup_any"]) == ("1", "0")          # src=M CDSpart ignored
    assert "#   Single-exon (no introns to evaluate): 1 (25.0%)" in text
    assert "#   All introns supported: 2 (50.0%)" in text
    assert "#   No intron support: 1 (25.0%)" in text
    assert "# Introns: 4" in text
    assert "#   Supported by both: 1 (25.0%)" in text
    assert len(proc.stderr.strip().splitlines()) == 1
    assert "4 coding transcripts" in proc.stderr and "2 hints of other sources ignored" in proc.stderr


def test_empty_hints(tmp_path: Path):
    (tmp_path / "genes.gff3").write_text(GENES)
    (tmp_path / "hints.gff").write_text("")
    proc = run("gene_support.py", "--gff3", "genes.gff3", "--hints", "hints.gff", "--out", "s.tsv", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    rows = support_rows(tmp_path / "s.tsv")
    assert len(rows) == 4
    for r in rows.values():
        assert r["introns_sup_any"] == "0" and r["exons_sup_any"] == "0"
    assert "#   Supported by any: 0 (0.0%)" in (tmp_path / "s.tsv").read_text()


def test_combined_source_and_overlap_index():
    gs = load("gene_support")
    assert gs.source_class("src=C;pri=4") == "combined"
    assert gs.classify({"combined"}) == (True, True)
    assert gs.source_class("grp=x;pri=4") == "other"
    index = gs.OverlapIndex([(10, 20), (100, 500), (30, 40)])
    assert index.overlaps(20, 25) and index.overlaps(450, 600) and index.overlaps(35, 99)
    assert not index.overlaps(21, 29) and not index.overlaps(41, 99) and not index.overlaps(501, 900)
    assert not index.overlaps(1, 9)
    assert gs.cds_introns([(1, 10), (11, 20), (31, 40)]) == [(21, 30)]


# ---------------------------------------------------------------- gene_set_statistics.py

STATS_GFF = gff(
    # G1: two transcripts, 3 and 2 exons; CDS 200 and 100 bp; span 500 each
    "c1 T gene 1 500 . + . ID=G1;gene_biotype=protein_coding",
    "c1 T mRNA 1 500 . + . ID=G1.t1;Parent=G1",
    "c1 T exon 1 100 . + . ID=G1.t1.exon1;Parent=G1.t1",
    "c1 T exon 201 300 . + . ID=G1.t1.exon2;Parent=G1.t1",
    "c1 T exon 401 500 . + . ID=G1.t1.exon3;Parent=G1.t1",
    "c1 T five_prime_UTR 1 50 . + . ID=G1.t1.utr5p1;Parent=G1.t1",
    "c1 T CDS 51 100 . + 0 ID=G1.t1.cds;Parent=G1.t1",
    "c1 T CDS 201 300 . + 1 ID=G1.t1.cds;Parent=G1.t1",
    "c1 T CDS 401 450 . + 0 ID=G1.t1.cds;Parent=G1.t1",
    "c1 T three_prime_UTR 451 500 . + . ID=G1.t1.utr3p1;Parent=G1.t1",
    "c1 T mRNA 1 500 . + . ID=G1.t2;Parent=G1",
    "c1 T exon 1 100 . + . ID=G1.t2.exon1;Parent=G1.t2",
    "c1 T exon 401 500 . + . ID=G1.t2.exon2;Parent=G1.t2",
    "c1 T CDS 51 100 . + 0 ID=G1.t2.cds;Parent=G1.t2",
    "c1 T CDS 401 450 . + 1 ID=G1.t2.cds;Parent=G1.t2",
    # G2: single exon, CDS 600 bp
    "c1 T gene 1001 1600 . + . ID=G2;gene_biotype=protein_coding",
    "c1 T mRNA 1001 1600 . + . ID=G2.t1;Parent=G2",
    "c1 T exon 1001 1600 . + . ID=G2.t1.exon1;Parent=G2.t1",
    "c1 T CDS 1001 1600 . + 0 ID=G2.t1.cds;Parent=G2.t1",
    # G3: minus strand, 2 exons, CDS 200 bp, span 300
    "c1 T gene 2001 2300 . - . ID=G3;gene_biotype=protein_coding",
    "c1 T mRNA 2001 2300 . - . ID=G3.t1;Parent=G3",
    "c1 T exon 2001 2100 . - . ID=G3.t1.exon1;Parent=G3.t1",
    "c1 T exon 2201 2300 . - . ID=G3.t1.exon2;Parent=G3.t1",
    "c1 T CDS 2001 2100 . - 2 ID=G3.t1.cds;Parent=G3.t1",
    "c1 T CDS 2201 2300 . - 0 ID=G3.t1.cds;Parent=G3.t1",
)

# Hand-counted: 3 genes, 4 transcripts (4/3 = 1.33 per gene), exons 3+2+1+2 = 8 (2.00 per
# transcript), 1 single-exon (25.0 %) and 3 multi-exon (75.0 %); CDS lengths 200, 100, 600,
# 200 (mean 275, upper median 200); spans 500, 500, 600, 300 (mean 475, upper median 500);
# introns per gene (2+1)//2 = 1, 0, 1.
EXPECTED_TEXT = (
    "Paludamentum Gene Set Statistics\n"
    + "=" * 50 + "\n\n"
    + "Genes:                         3\n"
    + "Transcripts:                   4\n"
    + "Mean transcripts/gene:      1.33\n"
    + "Mean exons/transcript:      2.00\n"
    + "\n"
    + "Single-exon transcripts:       1 (25.0%)\n"
    + "Multi-exon transcripts:        3 (75.0%)\n"
    + "\n"
    + "Mean CDS length:             275 bp\n"
    + "Mean genomic span:           475 bp\n"
    + "Median CDS length:           200 bp\n"
    + "Median genomic span:         500 bp\n"
)


def test_statistics_text(tmp_path: Path):
    (tmp_path / "g.gff3").write_text(STATS_GFF)
    proc = run("gene_set_statistics.py", "--gff3", "g.gff3", "--out-dir", "out", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert (tmp_path / "out" / "gene_set_statistics.txt").read_text() == EXPECTED_TEXT


def test_statistics_values():
    gss = load("gene_set_statistics")
    gff3_lib = load("gff3_lib")
    import io
    features = list(gff3_lib.read_features(io.StringIO(STATS_GFF)))
    stats = gss.compute_statistics(gff3_lib.build_annotation(features))
    assert stats["tx_per_gene"] == [2, 1, 1]
    assert stats["exon_counts"] == [3, 2, 1, 2]
    assert stats["cds_lengths"] == [200, 100, 600, 200]
    assert stats["genomic_spans"] == [500, 500, 600, 300]
    assert stats["introns_per_gene"] == [1, 0, 1]


def test_statistics_empty_gene_set(tmp_path: Path):
    (tmp_path / "g.gff3").write_text("##gff-version 3\n")
    proc = run("gene_set_statistics.py", "--gff3", "g.gff3", "--out-dir", "out", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    text = (tmp_path / "out" / "gene_set_statistics.txt").read_text()
    assert "Genes:                         0\n" in text and "(NA)" in text


def test_plots_written(tmp_path: Path):
    pytest.importorskip("matplotlib")
    (tmp_path / "g.gff3").write_text(STATS_GFF)
    (tmp_path / "genes.gff3").write_text(GENES)
    (tmp_path / "hints.gff").write_text(HINTS)
    assert run("gene_support.py", "--gff3", "genes.gff3", "--hints", "hints.gff", "--out", "s.tsv",
               cwd=tmp_path).returncode == 0
    proc = run("gene_set_statistics.py", "--gff3", "g.gff3", "--out-dir", "with", "--support", "s.tsv",
               cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    names = {"gene_set_statistics.txt", "isoform_and_exon_structure.png", "transcript_lengths.png",
             "introns_per_gene.png", "evidence_support.png"}
    assert {p.name for p in (tmp_path / "with").iterdir()} == names
    for name in names - {"gene_set_statistics.txt"}:
        assert (tmp_path / "with" / name).read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"

    proc = run("gene_set_statistics.py", "--gff3", "g.gff3", "--out-dir", "without", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert {p.name for p in (tmp_path / "without").iterdir()} == names - {"evidence_support.png"}

    # an empty gene set and a support table without rows still give every plot
    (tmp_path / "empty.gff3").write_text("##gff-version 3\n")
    (tmp_path / "hints0.gff").write_text("")
    assert run("gene_support.py", "--gff3", "empty.gff3", "--hints", "hints0.gff", "--out", "s0.tsv",
               cwd=tmp_path).returncode == 0
    proc = run("gene_set_statistics.py", "--gff3", "empty.gff3", "--out-dir", "empty", "--support", "s0.tsv",
               cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert {p.name for p in (tmp_path / "empty").iterdir()} == names


def test_evidence_fractions_skip_single_exon():
    gss = load("gene_set_statistics")
    data = [{"num_introns": "0", "introns_sup_any": "0"}, {"num_introns": "4", "introns_sup_any": "1"},
            {"num_introns": "2", "introns_sup_any": "2"}]
    assert gss.support_fractions(data) == [0.25, 1.0]
    with pytest.raises(ValueError):
        gss.support_fractions([{"transcript_id": "x"}])
