"""bin/cap_isoforms_per_gene.py on toy GTF and GFF3 (regression: GTF lost every
transcript, exon and CDS line because only key=value attributes were parsed)."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "bin" / "cap_isoforms_per_gene.py"


def cap(max_n: int, text: str, tmp_path: Path, stdin: bool = False) -> subprocess.CompletedProcess:
    if stdin:
        return subprocess.run([sys.executable, str(SCRIPT), "--max", str(max_n), "-"],
                              input=text, capture_output=True, text=True)
    path = tmp_path / "in.annot"
    path.write_text(text)
    return subprocess.run([sys.executable, str(SCRIPT), "--max", str(max_n), str(path)],
                          capture_output=True, text=True)


def line(ftype: str, start: int, end: int, attrs: str) -> str:
    phase = "0" if ftype == "CDS" else "."
    return f"chr1\tsrc\t{ftype}\t{start}\t{end}\t.\t+\t{phase}\t{attrs}\n"


def gtf_tx(tid: str, cds_len: int, with_transcript_line: bool = True, gid: str = "g") -> str:
    a = f'gene_id "{gid}"; transcript_id "{tid}";'
    out = line("transcript", 1, 1000, a) if with_transcript_line else ""
    return out + line("exon", 1, 1000, a) + line("CDS", 1, cds_len, a)


GTF_GENE = '##comment\nchr1\tsrc\tgene\t1\t1000\t.\t+\t.\tgene_id "g";\n'


def gtf_three(tmp_path: Path) -> str:
    return (GTF_GENE + gtf_tx("g.t1", 300) + gtf_tx("g.t2", 200, with_transcript_line=False)
            + gtf_tx("g.t3", 100))


def test_gtf_unchanged_below_cap(tmp_path: Path) -> None:
    text = gtf_three(tmp_path)
    proc = cap(5, text, tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout == text


def test_gtf_keeps_two_longest(tmp_path: Path) -> None:
    proc = cap(2, gtf_three(tmp_path), tmp_path)
    assert proc.returncode == 0, proc.stderr
    expected = GTF_GENE + gtf_tx("g.t1", 300) + gtf_tx("g.t2", 200, with_transcript_line=False)
    assert proc.stdout == expected
    assert "g.t3" not in proc.stdout


def test_gtf_transcript_without_transcript_line_is_dropped_by_cds_length(tmp_path: Path) -> None:
    text = GTF_GENE + gtf_tx("g.t1", 100) + gtf_tx("g.t2", 300, with_transcript_line=False)
    proc = cap(1, text, tmp_path)
    assert proc.stdout == GTF_GENE + gtf_tx("g.t2", 300, with_transcript_line=False)


def test_gtf_gene_and_comment_lines_always_written(tmp_path: Path) -> None:
    text = "# header\n" + GTF_GENE + gtf_tx("g.t1", 300) + gtf_tx("g.t2", 100)
    out = cap(1, text, tmp_path).stdout
    assert out.startswith("# header\n##comment\n")
    assert "\tgene\t" in out
    assert "g.t2" not in out and "g.t1" in out


def test_gtf_genes_capped_independently(tmp_path: Path) -> None:
    text = (gtf_tx("a.t1", 300, gid="a") + gtf_tx("a.t2", 100, gid="a")
            + gtf_tx("b.t1", 50, gid="b") + gtf_tx("b.t2", 40, gid="b"))
    out = cap(1, text, tmp_path).stdout
    assert out == gtf_tx("a.t1", 300, gid="a") + gtf_tx("b.t1", 50, gid="b")


GFF_GENE = "##gff-version 3\nchr1\tsrc\tgene\t1\t1000\t.\t+\t.\tID=g\n"


def gff_tx(tid: str, cds_len: int) -> str:
    return (line("mRNA", 1, 1000, f"ID={tid};Parent=g") + line("exon", 1, 1000, f"Parent={tid}")
            + line("CDS", 1, cds_len, f"Parent={tid}"))


def test_gff3_caps_by_cds_length(tmp_path: Path) -> None:
    text = GFF_GENE + gff_tx("t1", 100) + gff_tx("t2", 300) + gff_tx("t3", 200)
    proc = cap(2, text, tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout == GFF_GENE + gff_tx("t1", 100)[:0] + gff_tx("t2", 300) + gff_tx("t3", 200)
    assert "ID=t1" not in proc.stdout and "Parent=t1" not in proc.stdout


def test_gff3_unchanged_below_cap(tmp_path: Path) -> None:
    text = GFF_GENE + gff_tx("t1", 100) + gff_tx("t2", 300)
    assert cap(10, text, tmp_path).stdout == text


def test_reads_stdin(tmp_path: Path) -> None:
    proc = cap(2, gtf_three(tmp_path), tmp_path, stdin=True)
    assert proc.returncode == 0, proc.stderr
    assert "g.t3" not in proc.stdout and "g.t1" in proc.stdout and "g.t2" in proc.stdout
