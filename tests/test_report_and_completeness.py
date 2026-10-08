"""Tests of completeness_summary.py, filter_proteins.py and paludamentum_report.py.

Fixtures in tests/data/postprocess/: BUSCO 6.1.0 short summaries (genome with
E:, proteins), compleasm summaries (genome with I:, proteins without), a
gffcompare.stats and a gene_set_statistics.txt from BRAKER4/GALBA2 test runs.
"""
from __future__ import annotations

import importlib.util
import json
import re
import subprocess
import sys
from html.parser import HTMLParser
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "bin"
DATA = ROOT / "tests" / "data" / "postprocess"


def load(name: str):
    spec = importlib.util.spec_from_file_location(name, BIN / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run(script: str, *args: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(BIN / script), *args], capture_output=True, text=True, cwd=cwd)


def table_rows(text: str) -> list[list[str]]:
    return [line.split("\t") for line in text.splitlines() if line and not line.startswith("#")]


# ---------------------------------------------------------------- completeness_summary.py

def test_completeness_of_real_summaries(tmp_path: Path):
    out = tmp_path / "completeness.tsv"
    proc = run("completeness_summary.py", "--out", str(out),
               "--busco-genome", str(DATA / "busco_genome_short_summary.txt"),
               "--busco-proteins", str(DATA / "busco_proteins_short_summary.txt"),
               "--compleasm-genome", str(DATA / "compleasm_genome_summary.txt"),
               "--compleasm-proteins", str(DATA / "compleasm_proteins_summary.txt"))
    assert proc.returncode == 0, proc.stderr
    text = out.read_text()
    assert table_rows(text) == [
        ["assessment", "genome_or_proteome", "complete", "single", "duplicated", "fragmented", "missing", "n"],
        ["BUSCO", "genome", "1.6", "1.6", "0", "0", "98.4", "129"],
        ["BUSCO", "proteome", "1.6", "1.6", "0", "0.8", "97.7", "129"],
        ["compleasm", "genome", "1.55", "1.55", "0", "0.78", "97.67", "129"],
        ["compleasm", "proteome", "2.33", "2.33", "0", "0.78", "96.9", "129"],
    ]
    comments = [line for line in text.splitlines() if line.startswith("#")]
    assert "# lineage BUSCO genome: eukaryota_odb12" in comments
    assert "# lineage compleasm proteome: eukaryota_odb12" in comments
    assert any("E (complete BUSCOs with internal stop codons) 50%" in c for c in comments)


def test_compleasm_complete_is_s_plus_d_and_fragmented_is_f_plus_i(tmp_path: Path):
    cs = load("completeness_summary")
    f = tmp_path / "summary.txt"
    f.write_text("## lineage: vertebrata_odb12\nS:90.00%, 90\nD:2.50%, 2\nF:1.25%, 1\nI:0.75%, 1\n"
                 "M:5.50%, 6\nN:100\n")
    res = cs.parse_compleasm(str(f))
    assert res["complete"] == pytest.approx(92.5)
    assert res["fragmented"] == pytest.approx(2.0)
    assert res["n"] == 100 and res["lineage"] == "vertebrata_odb12"


def test_busco5_line_without_e(tmp_path: Path):
    cs = load("completeness_summary")
    f = tmp_path / "short.txt"
    f.write_text("\tC:95.2%[S:93.1%,D:2.1%],F:1.2%,M:3.6%,n:255\n")
    res = cs.parse_busco(str(f))
    assert (res["complete"], res["single"], res["duplicated"], res["fragmented"], res["missing"], res["n"]) == \
        (95.2, 93.1, 2.1, 1.2, 3.6, 255)
    assert res["E"] is None and res["lineage"] is None


def test_placeholder_and_missing_inputs_are_skipped(tmp_path: Path):
    empty = tmp_path / "empty.txt"
    empty.write_text("")
    out = tmp_path / "c.tsv"
    proc = run("completeness_summary.py", "--out", str(out), "--busco-genome", str(empty),
               "--busco-proteins", str(tmp_path / "nonexistent.txt"),
               "--compleasm-proteins", str(DATA / "compleasm_proteins_summary.txt"))
    assert proc.returncode == 0, proc.stderr
    assert table_rows(out.read_text())[1:] == [["compleasm", "proteome", "2.33", "2.33", "0", "0.78", "96.9", "129"]]

    proc = run("completeness_summary.py", "--out", str(out))
    assert proc.returncode == 0
    assert table_rows(out.read_text()) == [["assessment", "genome_or_proteome", "complete", "single", "duplicated",
                                            "fragmented", "missing", "n"]]


def test_unparsable_summary_fails(tmp_path: Path):
    bad = tmp_path / "bad.txt"
    bad.write_text("BUSCO crashed\n")
    proc = run("completeness_summary.py", "--out", str(tmp_path / "c.tsv"), "--busco-genome", str(bad))
    assert proc.returncode == 1
    assert "bad.txt" in proc.stderr


# ---------------------------------------------------------------- filter_proteins.py

def test_filter_proteins_drops_long_and_strips_stop(tmp_path: Path):
    (tmp_path / "in.fa").write_text(">p1 desc\nMKV\nLL*\n>long\n" + "A" * 60 + "\n" + "A" * 40 + "*\n"
                                    ">p3\nMK*L.Q*\n")
    proc = run("filter_proteins.py", "--in", "in.fa", "--out", "out.fa", "--max-length", "100", "--strip-stop",
               cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert (tmp_path / "out.fa").read_text() == ">p1 desc\nMKVLL\n>p3\nMKXLXQ\n"
    assert "dropped 1" in proc.stderr and "long" in proc.stderr


def test_filter_proteins_keeps_stop_without_flag(tmp_path: Path):
    (tmp_path / "in.fa").write_text(">p1\nMK\nV*\n>p2\n" + "A" * 99 + "*\n")
    proc = run("filter_proteins.py", "--in", "in.fa", "--out", "out.fa", "--max-length", "100", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert (tmp_path / "out.fa").read_text() == ">p1\nMKV*\n>p2\n" + "A" * 99 + "*\n"   # 99 residues + stop: kept


# ---------------------------------------------------------------- paludamentum_report.py

class _Checker(HTMLParser):
    """Collects tags, checks nesting of the non-void elements."""
    VOID = {"meta", "img", "br", "hr", "link", "input"}

    def __init__(self):
        super().__init__()
        self.stack: list[str] = []
        self.tags: list[str] = []
        self.errors: list[str] = []

    def handle_starttag(self, tag, attrs):
        self.tags.append(tag)
        if tag not in self.VOID:
            self.stack.append(tag)

    def handle_endtag(self, tag):
        if not self.stack or self.stack[-1] != tag:
            self.errors.append(f"unexpected </{tag}>, open: {self.stack[-3:]}")
        else:
            self.stack.pop()


def check_html(text: str) -> _Checker:
    assert text.startswith("<!DOCTYPE html>")
    checker = _Checker()
    checker.feed(text)
    assert checker.errors == []
    assert checker.stack == []
    for tag in ("html", "head", "title", "style", "body"):
        assert tag in checker.tags
    return checker


def no_remote_loads(text: str) -> None:
    """src= never remote; href=http only inside the references list."""
    assert not re.search(r"""src\s*=\s*["']?\s*(https?:)?//""", text)
    assert "@import" not in text and "url(http" not in text
    refs = text.split('<section id="references">', 1)
    before = refs[0]
    assert 'href="http' not in before
    if len(refs) == 2:
        after_refs = refs[1].split("</section>", 1)[1]
        assert 'href="http' not in after_refs


def test_report_on_empty_dir(tmp_path: Path):
    staged = tmp_path / "staged"
    staged.mkdir()
    proc = run("paludamentum_report.py", "--dir", str(staged), "--out", str(tmp_path / "report.html"))
    assert proc.returncode == 0, proc.stderr
    text = (tmp_path / "report.html").read_text()
    check_html(text)
    assert "<title>Paludamentum report</title>" in text
    assert "<section" not in text
    assert not (tmp_path / "completeness.png").exists()
    no_remote_loads(text)


def test_report_on_missing_dir_and_bad_run_info(tmp_path: Path):
    (tmp_path / "run_info.json").write_text("{not json")
    proc = run("paludamentum_report.py", "--dir", str(tmp_path / "nope"), "--out", str(tmp_path / "r.html"),
               "--run-info", str(tmp_path / "run_info.json"))
    assert proc.returncode == 0, proc.stderr
    check_html((tmp_path / "r.html").read_text())


def test_report_with_empty_placeholder_files(tmp_path: Path):
    staged = tmp_path / "staged"
    for rel in ("x_evidence.gff3", "citations.md", "qc/sanity_filter.tsv", "qc/utr_report.tsv",
                "qc/completeness.tsv", "qc/gene_set_statistics.txt", "qc/gene_support.tsv",
                "qc/software_versions.tsv", "qc/omark_summary.txt", "qc/gffcompare.stats",
                "qc/evidence_support.png", "qc/fantasia/fantasia_summary.txt", "ncrna/rRNA.gff3"):
        (staged / rel).parent.mkdir(parents=True, exist_ok=True)
        (staged / rel).write_text("")
    proc = run("paludamentum_report.py", "--dir", str(staged), "--out", str(tmp_path / "r.html"))
    assert proc.returncode == 0, proc.stderr
    text = (tmp_path / "r.html").read_text()
    check_html(text)
    sections = re.findall(r'<section id="([^"]+)"', text)
    assert sections == ["output-files"]          # only the file listing; every other section skipped


def make_png(path: Path) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(1, 1), dpi=20)
    ax.plot([0, 1], [0, 1])
    fig.savefig(path)
    plt.close(fig)


def full_staged_dir(tmp_path: Path) -> Path:
    staged = tmp_path / "staged"
    qc = staged / "qc"
    (qc / "fantasia").mkdir(parents=True)
    (staged / "ncrna").mkdir()
    stem = "tiberius_evidence"
    gff = ("##gff-version 3\nc1\tTiberius\tgene\t1\t90\t.\t+\t.\tID=g1\n"
           "c1\tTiberius\tmRNA\t1\t90\t.\t+\t.\tID=g1.t1;Parent=g1\n")
    (staged / f"{stem}.gff3").write_text(gff)
    (staged / f"{stem}.gtf").write_text('c1\tTiberius\tCDS\t1\t90\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n')
    (staged / f"{stem}_proteins.fa").write_text(">g1.t1\nMKV\n")
    (staged / f"{stem}_cds.fa").write_text(">g1.t1\nATGAAAGTT\n")
    for suffix in ("_with_ncRNA.gff3", "_go.gff3", "_with_ncRNA_go.gff3"):
        (staged / f"{stem}{suffix}").write_text(gff)
    (staged / "hintsfile.gff").write_text("c1\tProtHint\tintron\t10\t20\t1\t+\t.\tsrc=P;pri=4\n")
    (staged / "params.yaml").write_text("genome: /abs/genome.fa\nmode: rnaseq\n")
    (staged / "citations.md").write_text(
        "# References of this Paludamentum run\n\nPaludamentum 0.5.0, mode `rnaseq`, gene finder: tiberius.\n\n"
        "This file lists the software <and> data.\nPlease cite them.\n\n"
        "- **Tiberius** (gene prediction). Gabriel L. Accurate *ab initio* gene prediction. "
        "https://doi.org/10.1093/bioinformatics/btae685\n"
        "- **Nextflow** (workflow engine). Di Tommaso P. https://doi.org/10.1038/nbt.3820\n")
    (qc / "sanity_filter.tsv").write_text(
        "# transcripts in: 10\n# removed: 2\n# kept: 8\ngene_id\ttranscript_id\taction\treason\n"
        "g2\tg2.t1\tremoved\tinternal_stop\ng3\tg3.t1\tremoved\tinternal_stop\ng4\tg4.t1\tkept\textended_stop\n")
    (qc / "utr_report.tsv").write_text(
        "# transcripts without UTRs: 3\ntranscript_id\tmatched_by\tfive_prime_utr_bp\tthree_prime_utr_bp\n"
        "t1\tshort\t100\t200\nt2\tlong\t50\t0\nt3\tnone\t0\t0\n")
    (qc / "completeness.tsv").write_text(
        "# lineage BUSCO genome: eukaryota_odb12\n"
        "assessment\tgenome_or_proteome\tcomplete\tsingle\tduplicated\tfragmented\tmissing\tn\n"
        "BUSCO\tgenome\t95.2\t93.1\t2.1\t1.2\t3.6\t255\n"
        "compleasm\tproteome\t90\t88\t2\t4\t6\t255\n")
    (qc / "gene_set_statistics.txt").write_text((DATA / "gene_set_statistics.txt").read_text() + "Note: 287 <x>\n")
    make_png(tmp_path / "tiny.png")
    for name in ("isoform_and_exon_structure.png", "transcript_lengths.png", "introns_per_gene.png",
                 "evidence_support.png", "fantasia/fantasia_go_categories.png"):
        (qc / name).write_bytes((tmp_path / "tiny.png").read_bytes())
    (qc / "gene_support.tsv").write_text(   # format of bin/gene_support.py (BRAKER4 column names)
        "# Gene support summary\n# Transcripts: 3\n# Introns: 5\n#   Supported by both: 0 (0.0%)\n"
        "gene_id\ttranscript_id\tchrom\tstrand\tnum_introns\tintrons_sup_rnaseq\tpct_introns_sup_rnaseq\t"
        "introns_sup_protein\tpct_introns_sup_protein\tsupport_class\n"
        "g1\tt1\tc1\t+\t3\t3\t100.0\t0\t0.0\tfull\ng2\tt2\tc1\t+\t2\t0\t0.0\t0\t0.0\tnone\n"
        "g3\tt3\tc1\t-\t0\t0\tNA\t1\tNA\tpartial\n")
    (staged / "methods.md").write_text(
        "# Methods of this Paludamentum run\n\nGenes were predicted with **Tiberius** (Gabriel et al., 2024).\n"
        "Second line of the *ab initio* paragraph.\n\nA second paragraph.\n")
    # OMArk printed its usage over several lines in older runs: those lines are left out
    (qc / "software_versions.tsv").write_text(
        "tool\tversion\timage\nPaludamentum\t0.5.0\t-\n-d/--database\t/img/omark_0.4.1.sif\nETE_NCBI_DB]\n"
        "busco\t6.1.0\tezlabgva/busco:v6.1.0_cv2\nOMArk\t0.4.1\t/home/u/images/omark_0.4.1.sif\n")
    (qc / "omark_summary.txt").write_text("#The selected clade was Vertebrata\nS:90%\n")
    (qc / "gffcompare.stats").write_text((DATA / "gffcompare.stats").read_text())
    (qc / "fantasia" / "fantasia_summary.txt").write_text("Proteins with GO terms: 250\n")
    (staged / "ncrna" / "rRNA.gff3").write_text(
        "##gff-version 3\nc1\tbarrnap:0.9\trRNA\t1\t100\t1e-10\t+\t.\tName=5S_rRNA;product=5S ribosomal RNA\n"
        "c1\tbarrnap:0.9\trRNA\t200\t300\t1e-10\t+\t.\tName=5S_rRNA;product=5S ribosomal RNA\n"
        "c2\tbarrnap:0.9\trRNA\t1\t1800\t0\t-\t.\tName=18S_rRNA;product=18S ribosomal RNA\n")
    (staged / "ncrna" / "tRNAs.gff3").write_text(
        "##gff-version 3\nc1\ttRNAscan-SE\tgene\t1\t72\t.\t+\t.\tID=trna1\n"
        "c1\ttRNAscan-SE\ttRNA\t1\t72\t60\t+\t.\tID=trna1.t;Parent=trna1\n"
        "c1\ttRNAscan-SE\texon\t1\t72\t.\t+\t.\tParent=trna1.t\n")
    (staged / "ncrna" / "ncRNAs_infernal.gff3").write_text(
        "c1\tinfernal\tgene\t1\t100\t.\t+\t.\tID=n1\nc1\tinfernal\tsnoRNA\t1\t100\t.\t+\t.\tID=n1.t;Parent=n1\n")
    (staged / "ncrna" / "lncRNAs.gff3").write_text("")
    return staged


def test_report_with_every_file(tmp_path: Path):
    staged = full_staged_dir(tmp_path)
    (tmp_path / "run_info.json").write_text(json.dumps({
        "version": "0.5.0", "mode": "rnaseq", "genefinder": "tiberius", "model": "vertebrates", "hc": "drusilla",
        "stem": "tiberius_evidence", "busco_lineage": "eukaryota_odb12", "outdir": "/abs/out"}))
    out = tmp_path / "out" / "report.html"
    out.parent.mkdir()
    proc = run("paludamentum_report.py", "--dir", str(staged), "--out", str(out),
               "--run-info", str(tmp_path / "run_info.json"))
    assert proc.returncode == 0, proc.stderr
    text = out.read_text()
    check_html(text)
    no_remote_loads(text)
    assert "<title>Paludamentum report – tiberius_evidence</title>" in text
    assert re.findall(r'<section id="([^"]+)"', text) == [
        "run-summary", "methods", "output-files", "gene-set-statistics", "completeness", "evidence-support",
        "sanity-filter", "utrs", "ncrna", "omark", "gffcompare", "fantasia", "software-versions", "references"]
    assert (out.parent / "completeness.png").stat().st_size > 0
    # logo, 3 statistics plots, completeness, evidence support, ncRNA, FANTASIA
    assert text.count('src="data:image/png;base64,') == 8
    assert '<div class="logo"><img src="data:image/png;base64,' in text
    # methods: paragraphs, inline Markdown, no second title
    methods = text.split('id="methods"')[1].split("</section>")[0]
    assert methods.count("<p>") == 2 and "<strong>Tiberius</strong> (Gabriel et al., 2024)" in methods
    assert "<em>ab initio</em>" in methods and "Methods of this Paludamentum run" not in methods
    # escaping of free text
    assert "Mean exons/transcript:      5.33" in text
    assert "287 &lt;x&gt;" in text and "software &lt;and&gt; data" in text
    # run summary, output files: the staged files, and the chart the report wrote (published to qc/)
    assert "<td>drusilla</td>" in text and "<td>tiberius_evidence_with_ncRNA_go.gff3</td>" in text
    files = text.split('id="output-files"')[1].split("</section>")[0]
    listed = re.findall(r"<tr><td>([^<]+)</td>", files)
    assert listed[:4] == ["tiberius_evidence.gff3", "tiberius_evidence.gtf", "tiberius_evidence_proteins.fa",
                          "tiberius_evidence_cds.fa"]
    assert {"hintsfile.gff", "params.yaml", "qc/completeness.png", "qc/completeness.tsv"} <= set(listed)
    size = load("paludamentum_report").human_size((out.parent / "completeness.png").stat().st_size)
    assert f'<td>qc/completeness.png</td><td>completeness chart</td><td class="num">{size}</td>' in files
    # sanity filter count table
    assert "<td>removed</td><td>internal_stop</td><td class=\"num\">2</td>" in text
    # UTRs: 2 of 3 matched, mean 5' UTR (100 + 50) / 2, mean 3' UTR (200 + 0) / 2
    assert "2 of 3 transcripts matched" in text
    assert "<td>5&#x27; UTR</td><td class=\"num\">2</td><td class=\"num\">75.0</td>" in text
    assert "<td>3&#x27; UTR</td><td class=\"num\">1</td><td class=\"num\">100.0</td>" in text
    # gene support: 3 rows; introns_supported_rnaseq > 0 in 1 row; support_class not numeric
    assert "3 transcripts in" in text
    support = text.split('id="evidence-support"')[1].split("</section>")[0]
    assert "<td>introns_sup_rnaseq</td><td class=\"num\">1</td>" in support
    assert "<td>introns_sup_protein</td><td class=\"num\">1</td>" in support
    assert "pct_introns" not in support and "<td>support_class</td>" not in support
    assert "<pre>Gene support summary\nTranscripts: 3\nIntrons: 5\n  Supported by both: 0 (0.0%)</pre>" in support
    # ncRNA: totals per tool and a chart instead of a table
    ncrna = text.split('id="ncrna"')[1].split("</section>")[0]
    assert "barrnap: 3; tRNAscan-SE: 1; Infernal (Rfam): 1." in ncrna
    assert "<table>" not in ncrna and ncrna.count("data:image/png;base64,") == 1
    # software versions: Paludamentum first, malformed lines left out, image files by name
    versions = text.split('id="software-versions"')[1].split("</section>")[0]
    assert "<th>Container image</th>" in versions
    assert re.findall(r"<tr><td>([^<]+)</td>", versions) == ["Paludamentum", "busco", "OMArk"]
    assert '<td title="/home/u/images/omark_0.4.1.sif"><code>omark_0.4.1.sif</code></td>' in versions
    # gffcompare table + full file
    assert "<td>Intron chain level</td><td class=\"num\">48.0</td><td class=\"num\">64.8</td>" in text
    assert "<details>" in text
    # references
    assert "<strong>Tiberius</strong>" in text and "<em>ab initio</em>" in text and "<code>rnaseq</code>" in text
    assert '<a href="https://doi.org/10.1038/nbt.3820">' in text
    assert "References of this Paludamentum run" not in text   # the section has its own heading


def test_completeness_png_path_option(tmp_path: Path):
    pytest.importorskip("matplotlib")
    staged = tmp_path / "staged"
    (staged / "qc").mkdir(parents=True)
    (staged / "qc" / "completeness.tsv").write_text(
        "assessment\tgenome_or_proteome\tcomplete\tsingle\tduplicated\tfragmented\tmissing\tn\n"
        "BUSCO\tproteome\t1.6\t1.6\t0\t0.8\t97.7\t129\n")
    png = tmp_path / "plots" / "c.png"
    proc = run("paludamentum_report.py", "--dir", str(staged), "--out", str(tmp_path / "r.html"),
               "--completeness-png", str(png))
    assert proc.returncode == 0, proc.stderr
    assert png.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    assert not (tmp_path / "completeness.png").exists()
    # the output files list the chart under its published name, with the size of the written file
    text = (tmp_path / "r.html").read_text()
    files = text.split('id="output-files"')[1].split("</section>")[0]
    assert "<td>qc/completeness.png</td><td>completeness chart</td>" in files
    assert "<td>qc/completeness.tsv</td>" in files


def test_output_files_without_chart(tmp_path: Path):
    """No completeness.tsv rows: no chart, and no qc/completeness.png row."""
    staged = tmp_path / "staged"
    (staged / "qc").mkdir(parents=True)
    (staged / "qc" / "completeness.tsv").write_text("# nothing assessed\n")
    (staged / "hintsfile.gff").write_text("")
    proc = run("paludamentum_report.py", "--dir", str(staged), "--out", str(tmp_path / "r.html"))
    assert proc.returncode == 0, proc.stderr
    text = (tmp_path / "r.html").read_text()
    assert re.findall(r'<section id="([^"]+)"', text) == ["output-files"]
    assert "qc/completeness.png" not in text and "<td>hintsfile.gff</td>" in text
    assert not (tmp_path / "completeness.png").exists()


def test_ncrna_counts_by_type_and_rrna_subunit(tmp_path: Path):
    rep = load("paludamentum_report")
    gff = tmp_path / "rRNA.gff3"
    gff.write_text("##gff-version 3\n"
                   "c1\tbarrnap\tgene\t1\t100\t.\t+\t.\tID=g1\n"
                   "c1\tbarrnap\trRNA\t1\t100\t.\t+\t.\tID=r1;Parent=g1;Name=x_evidence-rRNA_1_5S_rRNA\n"
                   "c1\tbarrnap\trRNA\t1\t100\t.\t+\t.\tID=r2;Name=x_evidence-rRNA_2_5.8S_rRNA\n"
                   "c1\tbarrnap\trRNA\t1\t100\t.\t+\t.\tID=r3;Name=x_evidence-rRNA_3_5S_rRNA\n"
                   "c1\tbarrnap\trRNA\t1\t100\t.\t+\t.\tID=r4\n"
                   "c1\tbarrnap\texon\t1\t100\t.\t+\t.\tParent=r4\n")
    assert rep.count_ncrna(str(gff)) == [("5S rRNA", 2), ("5.8S rRNA", 1), ("rRNA", 1)]


def test_png_shown_at_physical_size(tmp_path: Path):
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rep = load("paludamentum_report")
    fig, ax = plt.subplots(figsize=(5, 2))
    fig.savefig(tmp_path / "p.png", dpi=150)
    plt.close(fig)
    data = (tmp_path / "p.png").read_bytes()
    assert rep.png_width_inches(data) == pytest.approx(5, abs=0.01)
    assert 'style="width:494px"' in rep.png_figure(data, "x")     # 5 in * 96 px + 14 px padding and border
    assert "style=" not in rep.png_figure(b"not a png", "x")


def test_citations_rendering():
    rep = load("paludamentum_report")
    html = rep.render_citations("- **A** (x). Some *ab initio* ref. https://example.org/a_b.\n")
    assert html == ('<ul class="refs"><li><strong>A</strong> (x). Some <em>ab initio</em> ref. '
                    '<a href="https://example.org/a_b">https://example.org/a_b</a>.</li></ul>')


# ---------------------------------------------------------------- fantasia_summary.py

FANTASIA_HEADER = "query_accession,go_id,reliability_index,distance,go_description,category\n"


@pytest.mark.parametrize("rows", [
    "",                                                          # results.csv with a header only
    "g1.t1,GO:0005524,0.21,0.9,ATP binding,molecular_function\n"  # one term below the cutoff
    "g2.t1,GO:0005634,0.40,0.8,nucleus,cellular_component\n",
], ids=["no-rows", "all-below-cutoff"])
def test_fantasia_summary_writes_a_placeholder_png_when_no_go_term_passes(tmp_path: Path, rows: str):
    """fantasia_go_categories.png is a required output of FANTASIA_SUMMARY; it must exist
    even when nothing survives --min-score (a placeholder figure, not a crash)."""
    results = tmp_path / "results.csv"
    results.write_text(FANTASIA_HEADER + rows)
    out = tmp_path / "out"
    proc = run("fantasia_summary.py", "--results", str(results), "--out-dir", str(out), "--min-score", "0.5")
    assert proc.returncode == 0, proc.stderr
    png = out / "fantasia_go_categories.png"
    assert png.is_file() and png.stat().st_size > 0
    assert png.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    assert (out / "fantasia_go_terms.tsv").read_text() == "transcript_id\tgo_id\tgo_name\tgo_namespace\treliability_index\n"
    assert "... above cutoff:                 0" in (out / "fantasia_summary.txt").read_text()


def test_fantasia_summary_chart_when_terms_pass(tmp_path: Path):
    results = tmp_path / "results.csv"
    results.write_text(FANTASIA_HEADER
                       + "g1.t1,GO:0005524,0.91,0.1,ATP binding,molecular_function\n"
                       + "g2.t1,GO:0005634,0.77,0.2,nucleus,cellular_component\n")
    out = tmp_path / "out"
    proc = run("fantasia_summary.py", "--results", str(results), "--out-dir", str(out))
    assert proc.returncode == 0, proc.stderr
    assert (out / "fantasia_go_categories.png").stat().st_size > 0
    assert len((out / "fantasia_go_terms.tsv").read_text().splitlines()) == 3
