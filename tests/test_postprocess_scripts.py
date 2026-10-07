"""bin/gff3_lib.py, normalize_gff3.py, sanity_filter_gff3.py, longest_isoform.py
and validate_gff3.sh on toy data (docs/postprocessing.md, the GFF3 contract)."""
from __future__ import annotations

import importlib.util
import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "bin"
GT = os.environ.get("GT_BIN") or shutil.which("gt")
needs_gt = pytest.mark.skipif(GT is None, reason="GenomeTools gt not found (set GT_BIN)")

sys.path.insert(0, str(BIN))
import gff3_lib  # noqa: E402


def run(script: str, *args: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(BIN / script), *args], capture_output=True, text=True, cwd=cwd)


def rows(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text().splitlines() if line and not line.startswith("#")]


# ---------------------------------------------------------------- toy genome

SENSE = [c for c in gff3_lib.CODON_TABLE if gff3_lib.CODON_TABLE[c] != "*" and c != "ATG"]


def orf(rng: random.Random, n_codons: int, stop: str = "TAA", internal_stop_at: int | None = None) -> str:
    """ATG, n_codons sense codons, a stop codon; optionally a TGA at codon internal_stop_at."""
    codons = [rng.choice(SENSE) for _ in range(n_codons)]
    if internal_stop_at is not None:
        codons[internal_stop_at] = "TGA"
    return "ATG" + "".join(codons) + stop


def toy_genome() -> str:
    rng = random.Random(7)
    seq = [rng.choice("ACGT") for _ in range(3000)]

    def put(start: int, text: str) -> None:
        seq[start - 1:start - 1 + len(text)] = list(text)

    a = orf(rng, 39)                              # 123 nt: CDS 111-160 + 301-373
    put(111, a[:50])
    put(301, a[50:])
    put(601, gff3_lib.revcomp(orf(rng, 28, internal_stop_at=14)))   # B: internal stop, minus strand
    d = orf(rng, 19, stop="")                     # D: 60 nt without stop, TAG after it
    put(1001, d + "TAG")
    f1 = orf(rng, 28)                             # F.t1 1401-1490
    put(1401, f1)
    put(1601, orf(rng, 18, internal_stop_at=5))   # F.t3 1601-1660
    g = orf(rng, 28)                              # G: minus strand, CDS 2101-2150 + 2001-2040
    put(2101, gff3_lib.revcomp(g[:50]))
    put(2001, gff3_lib.revcomp(g[50:]))
    return "".join(seq)


TOY_GFF3 = """##gff-version 3
c1\tx\tgene\t101\t400\t.\t+\t.\tID=gA
c1\tx\tmRNA\t101\t400\t.\t+\t.\tID=gA.t1;Parent=gA
c1\tx\texon\t101\t160\t.\t+\t.\tID=e1;Parent=gA.t1
c1\tx\texon\t301\t400\t.\t+\t.\tID=e2;Parent=gA.t1
c1\tx\tCDS\t111\t160\t.\t+\t0\tID=c1;Parent=gA.t1
c1\tx\tCDS\t301\t373\t.\t+\t.\tID=c2;Parent=gA.t1
c1\tx\tUTR\t101\t110\t.\t+\t.\tParent=gA.t1
c1\tx\tUTR\t374\t400\t.\t+\t.\tParent=gA.t1
c1\tx\tgene\t601\t690\t.\t-\t.\tID=gB
c1\tx\tmRNA\t601\t690\t.\t-\t.\tID=gB.t1;Parent=gB
c1\tx\tCDS\t601\t690\t.\t-\t0\tParent=gB.t1
c1\tx\tgene\t801\t810\t.\t+\t.\tID=gC
c1\tx\tmRNA\t801\t810\t.\t+\t.\tID=gC.t1;Parent=gC
c1\tx\tCDS\t801\t810\t.\t+\t0\tParent=gC.t1
c1\tx\tgene\t1001\t1060\t.\t+\t.\tID=gD
c1\tx\tmRNA\t1001\t1060\t.\t+\t.\tID=gD.t1;Parent=gD
c1\tx\texon\t1001\t1060\t.\t+\t.\tParent=gD.t1
c1\tx\tCDS\t1001\t1060\t.\t+\t0\tParent=gD.t1
c1\tx\tgene\t1201\t1300\t.\t+\t.\tID=gE
c1\tx\ttranscript\t1201\t1300\t.\t+\t.\tID=gE.t1;Parent=gE
c1\tx\texon\t1201\t1300\t.\t+\t.\tParent=gE.t1
c1\tx\tgene\t1401\t1660\t.\t+\t.\tID=gF
c1\tx\tmRNA\t1401\t1490\t.\t+\t.\tID=gF.t1;Parent=gF
c1\tx\tCDS\t1401\t1490\t.\t+\t0\tParent=gF.t1
c1\tx\tmRNA\t1404\t1490\t.\t+\t.\tID=gF.t2;Parent=gF
c1\tx\tCDS\t1404\t1490\t.\t+\t0\tParent=gF.t2
c1\tx\tmRNA\t1601\t1660\t.\t+\t.\tID=gF.t3;Parent=gF
c1\tx\tCDS\t1601\t1660\t.\t+\t0\tParent=gF.t3
c1\tx\tgene\t2001\t2150\t.\t-\t.\tID=gG
c1\tx\tmRNA\t2001\t2150\t.\t-\t.\tID=gG.t1;Parent=gG
c1\tx\tCDS\t2001\t2040\t.\t-\t2\tParent=gG.t1
c1\tx\tCDS\t2101\t2150\t.\t-\t0\tParent=gG.t1
c1\tx\tstop_codon\t2001\t2003\t.\t-\t0\tParent=gG.t1
"""


@pytest.fixture
def toy(tmp_path: Path) -> Path:
    (tmp_path / "genome.fa").write_text(">c1 toy\n" + toy_genome() + "\n")
    (tmp_path / "in.gff3").write_text(TOY_GFF3)
    return tmp_path


# ---------------------------------------------------------------- gff3_lib

def test_attributes_are_decoded_and_encoded():
    attrs = gff3_lib.parse_attributes("ID=a%3Bb;Parent=p1,p2;Note=x y%2Cz;Ontology_term=GO:1,GO:2")
    assert attrs == {"ID": ["a;b"], "Parent": ["p1", "p2"], "Note": ["x y,z"], "Ontology_term": ["GO:1", "GO:2"]}
    assert gff3_lib.format_attributes(attrs) == "ID=a%3Bb;Parent=p1,p2;Note=x y%2Cz;Ontology_term=GO:1,GO:2"


@pytest.mark.parametrize("phase,length,expected", [(0, 50, 1), (1, 73, 0), (0, 40, 2), (2, 40, 1), (0, 3, 0)])
def test_next_phase(phase, length, expected):
    assert gff3_lib.next_phase(phase, length) == expected


def test_missing_parent_and_duplicate_id_are_errors(tmp_path: Path):
    path = tmp_path / "bad.gff3"
    path.write_text("##gff-version 3\nc1\tx\tgene\t1\t9\t.\t+\t.\tID=g\nc1\tx\tCDS\t1\t9\t.\t+\t0\tParent=nosuch\n")
    with pytest.raises(gff3_lib.GFF3Error, match="Parent nosuch of line 3: CDS"):
        gff3_lib.read_gff3(str(path))
    path.write_text("##gff-version 3\nc1\tx\tgene\t1\t9\t.\t+\t.\tID=g\nc1\tx\tgene\t11\t19\t.\t+\t.\tID=g\n")
    with pytest.raises(gff3_lib.GFF3Error, match="duplicate ID g"):
        gff3_lib.read_gff3(str(path))
    assert len(gff3_lib.read_gff3(str(path), strict=False).genes) == 1


def test_parts_of_one_cds_share_the_id_and_a_lonely_transcript_gets_a_gene(tmp_path: Path):
    path = tmp_path / "in.gff3"
    path.write_text("c1\tx\tmRNA\t1\t30\t.\t+\t.\tID=t\n"
                    "c1\tx\tCDS\t1\t9\t.\t+\t0\tID=t.cds;Parent=t\nc1\tx\tCDS\t21\t30\t.\t+\t0\tID=t.cds;Parent=t\n")
    ann = gff3_lib.read_gff3(str(path))
    assert [g.id for g in ann.genes] == ["t.gene"]
    assert [c.start for c in ann.genes[0].transcripts[0].cds] == [1, 21]


# ---------------------------------------------------------------- normalize_gff3.py

UNNORMALISED = """##gff-version 3
# a comment
c2\tsrc\tgene\t500\t900\t.\t-\t.\tID=g2;Name=two%2Cthree
c2\tsrc\ttranscript\t500\t900\t.\t-\t.\tID=g2.t1;Parent=g2
c2\tsrc\tCDS\t500\t600\t.\t-\t.\tID=x1;Parent=g2.t1
c2\tsrc\tCDS\t700\t850\t.\t-\t.\tID=x2;Parent=g2.t1
c2\tsrc\tUTR\t851\t900\t.\t-\t.\tParent=g2.t1
c2\tsrc\tintron\t601\t699\t.\t-\t.\tParent=g2.t1
c1\tsrc\tgene\t10\t99\t.\t+\t.\tID=g1;Note=a b
c1\tsrc\tmRNA\t10\t99\t.\t+\t.\tID=g1.t1;Parent=g1
c1\tsrc\texon\t10\t99\t.\t+\t.\tID=whatever;Parent=g1.t1
c1\tsrc\tCDS\t10\t99\t.\t+\t0\tID=whatever2;Parent=g1.t1
"""

NORMALISED = """##gff-version 3
c1\tsrc\tgene\t10\t99\t.\t+\t.\tID=g1;Note=a b;gene_biotype=protein_coding
c1\tsrc\tmRNA\t10\t99\t.\t+\t.\tID=g1.t1;Parent=g1
c1\tsrc\texon\t10\t99\t.\t+\t.\tID=g1.t1.exon1;Parent=g1.t1
c1\tsrc\tCDS\t10\t99\t.\t+\t0\tID=g1.t1.cds;Parent=g1.t1
c2\tsrc\tgene\t500\t900\t.\t-\t.\tID=g2;Name=two%2Cthree;gene_biotype=protein_coding
c2\tsrc\tmRNA\t500\t900\t.\t-\t.\tID=g2.t1;Parent=g2
c2\tsrc\texon\t500\t600\t.\t-\t.\tID=g2.t1.exon1;Parent=g2.t1
c2\tsrc\tCDS\t500\t600\t.\t-\t2\tID=g2.t1.cds;Parent=g2.t1
c2\tsrc\texon\t700\t900\t.\t-\t.\tID=g2.t1.exon2;Parent=g2.t1
c2\tsrc\tCDS\t700\t850\t.\t-\t0\tID=g2.t1.cds;Parent=g2.t1
c2\tsrc\tfive_prime_UTR\t851\t900\t.\t-\t.\tID=g2.t1.utr5p1;Parent=g2.t1
"""


def test_normalize_writes_the_contract(tmp_path: Path):
    (tmp_path / "in.gff3").write_text(UNNORMALISED)
    proc = run("normalize_gff3.py", "--gff3", "in.gff3", "--out", "out.gff3", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert (tmp_path / "out.gff3").read_text() == NORMALISED
    # idempotent
    proc = run("normalize_gff3.py", "--gff3", "out.gff3", "--out", "again.gff3", cwd=tmp_path)
    assert (tmp_path / "again.gff3").read_text() == NORMALISED


@pytest.mark.parametrize("bad,message", [
    ("c1\tx\tCDS\t5\t99\t.\t+\t0\tParent=g1.t1\n", "lies outside the exons"),
    ("c1\tx\texon\t50\t120\t.\t+\t.\tParent=g1.t1\n", "exons of g1.t1 overlap"),
    ("c1\tx\tCDS\t200\t210\t.\t+\t0\tParent=nope\n", "Parent nope"),
])
def test_normalize_names_the_broken_feature(bad: str, message: str, tmp_path: Path):
    (tmp_path / "in.gff3").write_text(NORMALISED + bad)
    proc = run("normalize_gff3.py", "--gff3", "in.gff3", "--out", "out.gff3", cwd=tmp_path)
    assert proc.returncode == 1
    assert message in proc.stderr


@needs_gt
def test_merge_annotations_output_passes_gt_after_normalisation(tmp_path: Path):
    from test_merge_annotations import CODING_GTF, NATIVE_GTF, NONCODING_GTF
    inputs = []
    for i, text in enumerate((CODING_GTF, NONCODING_GTF, NATIVE_GTF.replace("chr1", "chr3"))):
        (tmp_path / f"in{i}.gtf").write_text(text)
        inputs.append(str(tmp_path / f"in{i}.gtf"))
    merged = subprocess.run([sys.executable, str(BIN / "merge_annotations.py"), "--mode", "full", *inputs],
                            capture_output=True, text=True, check=True).stdout
    (tmp_path / "merged.gff3").write_text(merged)
    proc = run("normalize_gff3.py", "--gff3", "merged.gff3", "--out", "out.gff3", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    gt = subprocess.run([GT, "gff3validator", str(tmp_path / "out.gff3")], capture_output=True, text=True)
    assert gt.returncode == 0, gt.stderr


# ---------------------------------------------------------------- sanity_filter_gff3.py

def report(path: Path) -> dict[tuple[str, str], str]:
    out = {}
    for row in rows(path)[1:]:
        out[(row[1], row[3])] = row[2]
    return out


def test_sanity_filter(toy: Path):
    proc = run("sanity_filter_gff3.py", "--gff3", "in.gff3", "--genome", "genome.fa",
               "--out", "out.gff3", "--report", "r.tsv", cwd=toy)
    assert proc.returncode == 0, proc.stderr
    assert report(toy / "r.tsv") == {
        ("gB.t1", "internal_stop"): "removed",
        ("gC.t1", "cds_length_mod3"): "removed",
        ("gD.t1", "extended_stop"): "kept",
        ("gE.t1", "no_cds"): "removed",
        ("gF.t2", "non_atg_start"): "kept",
        ("gF.t3", "internal_stop"): "removed",
    }
    assert "# transcripts removed: 4; genes left: 4" in (toy / "r.tsv").read_text()
    out = rows(toy / "out.gff3")
    assert [r[8].split(";")[0] for r in out if r[2] == "gene"] == ["ID=gA", "ID=gD", "ID=gF", "ID=gG"]
    assert [r[8].split(";")[0] for r in out if r[2] == "mRNA"] == \
        ["ID=gA.t1", "ID=gD.t1", "ID=gF.t1", "ID=gF.t2", "ID=gG.t1"]
    # the stop codon of D is added to its CDS, exon, mRNA and gene
    d = [r for r in out if "gD" in r[8]]
    assert {(r[2], r[3], r[4]) for r in d} == {("gene", "1001", "1063"), ("mRNA", "1001", "1063"),
                                              ("exon", "1001", "1063"), ("CDS", "1001", "1063")}
    # phases recomputed; UTRs typed; F spans only its kept transcripts
    a_cds = [(r[3], r[7]) for r in out if r[2] == "CDS" and "gA.t1" in r[8]]
    assert a_cds == [("111", "0"), ("301", "1")]
    assert [r[2] for r in out if "gA.t1.utr" in r[8]] == ["five_prime_UTR", "three_prime_UTR"]
    assert [(r[3], r[4]) for r in out if r[2] == "gene" and "ID=gF" in r[8]] == [("1401", "1490")]
    g_cds = [(r[3], r[7]) for r in out if r[2] == "CDS" and "gG.t1" in r[8]]
    assert g_cds == [("2001", "1"), ("2101", "0")]
    assert not [r for r in out if r[2] == "stop_codon"]
    # every kept transcript translates without internal stop and ends with one
    ann = gff3_lib.read_gff3(str(toy / "out.gff3"))
    genome = gff3_lib.read_fasta(str(toy / "genome.fa"))
    for _gene, tx in ann.transcripts():
        protein = gff3_lib.translate(gff3_lib.spliced_cds(tx, genome))
        assert protein.endswith("*") and "*" not in protein[:-1], tx.id


def test_sanity_filter_keep_reports_only(toy: Path):
    proc = run("sanity_filter_gff3.py", "--gff3", "in.gff3", "--genome", "genome.fa",
               "--out", "out.gff3", "--report", "r.tsv", "--keep", cwd=toy)
    assert proc.returncode == 0, proc.stderr
    flagged = report(toy / "r.tsv")
    assert flagged[("gB.t1", "internal_stop")] == "flagged"
    assert flagged[("gD.t1", "extended_stop")] == "flagged"
    out = rows(toy / "out.gff3")
    assert len([r for r in out if r[2] in ("mRNA", "transcript")]) == 9
    assert ("CDS", "1001", "1060") in {(r[2], r[3], r[4]) for r in out}


@needs_gt
def test_sanity_filter_output_passes_the_validator(toy: Path):
    run("sanity_filter_gff3.py", "--gff3", "in.gff3", "--genome", "genome.fa",
        "--out", "out.gff3", "--report", "r.tsv", cwd=toy)
    gt = subprocess.run([GT, "gff3validator", str(toy / "out.gff3")], capture_output=True, text=True)
    assert gt.returncode == 0, gt.stderr


def test_five_prime_partial_phase_counts_in_the_mod3_check(tmp_path: Path):
    """A 5' partial transcript starts mid-codon: the phase of the first CDS is
    not coding sequence (spliced_cds strips it) and must not fail check a."""
    rng = random.Random(3)
    coding = "".join(rng.choice(SENSE) for _ in range(19)) + "TAA"   # 60 nt, in-frame stop
    plus = "GG" + coding + "ACGT" * 5                                # CDS 1-62, phase 2
    minus = "ACGT" * 5 + gff3_lib.revcomp("GG" + coding)             # CDS 21-82, phase 2
    (tmp_path / "genome.fa").write_text(f">c1\n{plus}\n>c2\n{minus}\n")
    (tmp_path / "in.gff3").write_text(
        "##gff-version 3\n"
        "c1\tx\tgene\t1\t62\t.\t+\t.\tID=gP\n"
        "c1\tx\tmRNA\t1\t62\t.\t+\t.\tID=gP.t1;Parent=gP\n"
        "c1\tx\tCDS\t1\t62\t.\t+\t2\tParent=gP.t1\n"
        "c2\tx\tgene\t21\t82\t.\t-\t.\tID=gM\n"
        "c2\tx\tmRNA\t21\t82\t.\t-\t.\tID=gM.t1;Parent=gM\n"
        "c2\tx\tCDS\t21\t82\t.\t-\t2\tParent=gM.t1\n"
        "c1\tx\tgene\t101\t162\t.\t+\t.\tID=gQ\n"                    # 62 - phase 1 = 61 nt: out of frame
        "c1\tx\tmRNA\t101\t162\t.\t+\t.\tID=gQ.t1;Parent=gQ\n"
        "c1\tx\tCDS\t101\t162\t.\t+\t1\tParent=gQ.t1\n")
    proc = run("sanity_filter_gff3.py", "--gff3", "in.gff3", "--genome", "genome.fa",
               "--out", "out.gff3", "--report", "r.tsv", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert report(tmp_path / "r.tsv") == {("gQ.t1", "cds_length_mod3"): "removed"}
    out = rows(tmp_path / "out.gff3")
    assert [r[8].split(";")[0] for r in out if r[2] == "mRNA"] == ["ID=gP.t1", "ID=gM.t1"]


def test_stop_codon_split_by_an_intron_is_not_added(tmp_path: Path):
    """A CDS that ends at an exon end followed by another exon keeps no_stop."""
    rng = random.Random(1)
    body = orf(rng, 19, stop="")       # 60 nt, no stop
    seq = list("".join(rng.choice("ACGT") for _ in range(400)))
    seq[100:160] = list(body)
    seq[160:163] = list("TAA")         # a stop right after the exon end, inside the intron
    (tmp_path / "genome.fa").write_text(">c1\n" + "".join(seq) + "\n")
    (tmp_path / "in.gff3").write_text(
        "c1\tx\tgene\t101\t300\t.\t+\t.\tID=g\nc1\tx\tmRNA\t101\t300\t.\t+\t.\tID=g.t1;Parent=g\n"
        "c1\tx\texon\t101\t160\t.\t+\t.\tParent=g.t1\nc1\tx\texon\t251\t300\t.\t+\t.\tParent=g.t1\n"
        "c1\tx\tCDS\t101\t160\t.\t+\t0\tParent=g.t1\nc1\tx\tthree_prime_UTR\t251\t300\t.\t+\t.\tParent=g.t1\n")
    proc = run("sanity_filter_gff3.py", "--gff3", "in.gff3", "--genome", "genome.fa",
               "--out", "out.gff3", "--report", "r.tsv", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert report(tmp_path / "r.tsv") == {("g.t1", "no_stop"): "kept"}


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("slack, extended", [(1, False), (2, False), (3, True)])
def test_stop_codon_reaching_into_an_intron_is_not_added(tmp_path: Path, strand: str, slack: int, extended: bool):
    """The exon holds `slack` bases after the CDS end and another exon follows.
    The stop codon read across the exon end is only added when all three bases
    are in the exon (slack 3). With slack 1 or 2 the codon is TA|G or T|GA with
    the intron's first bases, a stop that is not in the transcript."""
    rng = random.Random(slack)
    seq = list("".join(rng.choice("ACGT") for _ in range(400)))
    seq[100:160] = list(orf(rng, 19, stop=""))        # CDS 101..160, no stop
    seq[160:163] = list("TAG" if slack != 1 else "TGA")  # exon end 160+slack, then the intron
    seq[163:165] = list("GT")
    plus = "".join(seq)
    exon1 = (101, 160 + slack)
    cds = (101, 160)
    exon2 = (251, 300)
    if strand == "-":
        plus = gff3_lib.revcomp(plus)
        exon1, cds, exon2 = [(401 - b, 401 - a) for a, b in (exon1, cds, exon2)]
    (tmp_path / "genome.fa").write_text(">c1\n" + plus + "\n")
    (tmp_path / "in.gff3").write_text(
        f"c1\tx\tgene\t101\t300\t.\t{strand}\t.\tID=g\nc1\tx\tmRNA\t101\t300\t.\t{strand}\t.\tID=g.t1;Parent=g\n"
        f"c1\tx\texon\t{exon1[0]}\t{exon1[1]}\t.\t{strand}\t.\tParent=g.t1\n"
        f"c1\tx\texon\t{exon2[0]}\t{exon2[1]}\t.\t{strand}\t.\tParent=g.t1\n"
        f"c1\tx\tCDS\t{cds[0]}\t{cds[1]}\t.\t{strand}\t0\tParent=g.t1\n")
    proc = run("sanity_filter_gff3.py", "--gff3", "in.gff3", "--genome", "genome.fa",
               "--out", "out.gff3", "--report", "r.tsv", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    out = {(r[2], r[3], r[4]) for r in rows(tmp_path / "out.gff3")}
    if extended:
        assert report(tmp_path / "r.tsv") == {("g.t1", "extended_stop"): "kept"}
        new_cds = (cds[0], cds[1] + 3) if strand == "+" else (cds[0] - 3, cds[1])
        assert ("CDS", str(new_cds[0]), str(new_cds[1])) in out
    else:
        assert report(tmp_path / "r.tsv") == {("g.t1", "no_stop"): "kept"}
        assert ("CDS", str(cds[0]), str(cds[1])) in out
    assert ("exon", str(exon1[0]), str(exon1[1])) in out   # the exon never grows into the intron


def test_extended_stop_shortens_the_three_prime_utr(tmp_path: Path):
    rng = random.Random(2)
    seq = list("".join(rng.choice("ACGT") for _ in range(400)))
    seq[100:160] = list(orf(rng, 19, stop=""))
    seq[160:163] = list("TGA")
    (tmp_path / "genome.fa").write_text(">c1\n" + "".join(seq) + "\n")
    (tmp_path / "in.gff3").write_text(
        "c1\tx\tgene\t101\t200\t.\t+\t.\tID=g\nc1\tx\tmRNA\t101\t200\t.\t+\t.\tID=g.t1;Parent=g\n"
        "c1\tx\texon\t101\t200\t.\t+\t.\tParent=g.t1\nc1\tx\tCDS\t101\t160\t.\t+\t0\tParent=g.t1\n"
        "c1\tx\tthree_prime_UTR\t161\t200\t.\t+\t.\tParent=g.t1\n")
    proc = run("sanity_filter_gff3.py", "--gff3", "in.gff3", "--genome", "genome.fa",
               "--out", "out.gff3", "--report", "r.tsv", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    out = {(r[2], r[3], r[4]) for r in rows(tmp_path / "out.gff3")}
    assert ("CDS", "101", "163") in out and ("three_prime_UTR", "164", "200") in out
    assert ("exon", "101", "200") in out


# ---------------------------------------------------------------- longest_isoform.py

def test_longest_isoform(toy: Path):
    run("sanity_filter_gff3.py", "--gff3", "in.gff3", "--genome", "genome.fa",
        "--out", "filtered.gff3", "--report", "r.tsv", cwd=toy)
    proc = run("longest_isoform.py", "--gff3", "filtered.gff3", "--out", "li.gff3", cwd=toy)
    assert proc.returncode == 0, proc.stderr
    out = rows(toy / "li.gff3")
    assert [r[8].split(";")[0] for r in out if r[2] == "mRNA"] == ["ID=gA.t1", "ID=gD.t1", "ID=gF.t1", "ID=gG.t1"]
    assert [(r[3], r[4]) for r in out if r[2] == "gene" and "ID=gF" in r[8]] == [("1401", "1490")]


def test_longest_isoform_drops_noncoding_genes_and_keeps_the_first_on_ties(tmp_path: Path):
    (tmp_path / "in.gff3").write_text(
        "c1\tx\tgene\t1\t90\t.\t+\t.\tID=g\n"
        "c1\tx\tmRNA\t1\t60\t.\t+\t.\tID=g.t1;Parent=g\nc1\tx\tCDS\t1\t60\t.\t+\t0\tParent=g.t1\n"
        "c1\tx\tmRNA\t31\t90\t.\t+\t.\tID=g.t2;Parent=g\nc1\tx\tCDS\t31\t90\t.\t+\t0\tParent=g.t2\n"
        "c1\tx\tgene\t200\t300\t.\t+\t.\tID=n\n"
        "c1\tx\tlnc_RNA\t200\t300\t.\t+\t.\tID=n.t1;Parent=n\nc1\tx\texon\t200\t300\t.\t+\t.\tParent=n.t1\n")
    run("longest_isoform.py", "--gff3", "in.gff3", "--out", "li.gff3", cwd=tmp_path)
    out = rows(tmp_path / "li.gff3")
    assert [r[8] for r in out if r[2] == "mRNA"] == ["ID=g.t1;Parent=g"]
    assert not [r for r in out if r[2] == "lnc_RNA"]


# ---------------------------------------------------------------- validate_gff3.sh

@needs_gt
def test_validate_gff3(tmp_path: Path):
    if not shutil.which("gffread"):
        pytest.skip("gffread not found")
    env = dict(os.environ, PATH=f"{Path(GT).parent}{os.pathsep}{os.environ['PATH']}")
    (tmp_path / "good.gff3").write_text(NORMALISED)
    good = subprocess.run([str(BIN / "validate_gff3.sh"), "good.gff3"], cwd=tmp_path, env=env,
                          capture_output=True, text=True)
    assert good.returncode == 0, good.stderr
    (tmp_path / "bad.gff3").write_text(NORMALISED + "c1\tx\tCDS\t200\t210\t.\t+\t0\tParent=nope\n")
    bad = subprocess.run([str(BIN / "validate_gff3.sh"), "bad.gff3"], cwd=tmp_path, env=env,
                         capture_output=True, text=True)
    assert bad.returncode == 1 and "gt gff3validator rejects" in bad.stderr


# ---------------------------------------------------------------- prefix_trna_ids.sh

def test_trna_ids_get_the_stem_as_prefix(tmp_path: Path):
    (tmp_path / "t.gff3").write_text(
        "##gff-version 3\n"
        "chr1\ttRNAscan-SE\ttRNA\t10\t80\t60.1\t+\t.\tID=chr1.trna1;Name=chr1.tRNA1-AlaAGC;isotype=Ala\n"
        "chr1\ttRNAscan-SE\texon\t10\t80\t.\t+\t.\tID=chr1.trna1.exon1; Parent=chr1.trna1\n")
    out = subprocess.run([str(BIN / "prefix_trna_ids.sh"), "tiberius_evidence", "t.gff3"], cwd=tmp_path,
                         capture_output=True, text=True, check=True).stdout.splitlines()
    assert out[0] == "##gff-version 3"
    assert out[1].endswith("ID=tiberius_evidence-chr1.trna1;Name=chr1.tRNA1-AlaAGC;isotype=Ala")
    assert out[2].endswith("ID=tiberius_evidence-chr1.trna1.exon1;Parent=tiberius_evidence-chr1.trna1")
