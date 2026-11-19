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


# --- malformed transcripts: dropped with a warning on stderr -----------------

def gtf_parts(tid: str, cds: list[tuple[int, int]], gid: str = "g", extra: str = "") -> str:
    a = f'gene_id "{gid}"; transcript_id "{tid}";'
    out = line("transcript", 1, 1000, a) + extra.replace("ATTRS", a)
    for s, e in cds:
        out += line("CDS", s, e, a)
    return out


def gff_parts(tid: str, cds: list[tuple[int, int]]) -> str:
    out = line("mRNA", 1, 1000, f"ID={tid};Parent=g")
    for s, e in cds:
        out += line("CDS", s, e, f"Parent={tid}")
    return out


def test_gtf_overlapping_cds_dropped_with_warning(tmp_path: Path) -> None:
    good = gtf_parts("h.t1", [(1, 100), (201, 300)], gid="h")
    bad = gtf_parts("g.bad", [(1, 100), (50, 150)])
    proc = cap(10, GTF_GENE + bad + good, tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout == GTF_GENE + good
    assert "g.bad" not in proc.stdout
    assert "dropped 1 malformed transcript(s)" in proc.stderr
    assert "g.bad" in proc.stderr


def test_gtf_inverted_intron_row_dropped(tmp_path: Path) -> None:
    a = 'gene_id "g"; transcript_id "g.inv";'
    bad = gtf_parts("g.inv", [(1, 100), (201, 300)], extra=line("intron", 5276, 389, "ATTRS"))
    assert a in bad
    good = gtf_parts("g.ok", [(1, 100)])
    proc = cap(10, GTF_GENE + bad + good, tmp_path)
    assert proc.stdout == GTF_GENE + good
    assert "g.inv" in proc.stderr and "g.ok" not in proc.stderr


def test_gff3_overlapping_cds_drops_mrna_and_children(tmp_path: Path) -> None:
    good = gff_parts("t_ok", [(1, 100), (201, 300)])
    bad = gff_parts("t_bad", [(1, 100), (100, 200)])
    proc = cap(10, GFF_GENE + bad + good, tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout == GFF_GENE + good
    assert "t_bad" not in proc.stdout
    assert "t_bad" in proc.stderr


def test_adjacent_cds_parts_not_malformed(tmp_path: Path) -> None:
    text = GTF_GENE + gtf_parts("g.t1", [(1, 100), (101, 200)])
    proc = cap(10, text, tmp_path)
    assert proc.stdout == text
    assert proc.stderr == ""


def test_max_picks_among_well_formed(tmp_path: Path) -> None:
    bad = gtf_parts("g.long", [(1, 500), (400, 900)])
    short = gtf_parts("g.short", [(1, 50)])
    proc = cap(1, GTF_GENE + bad + short, tmp_path)
    assert proc.stdout == GTF_GENE + short
    assert "g.long" in proc.stderr


def test_warning_truncated_after_ten_ids(tmp_path: Path) -> None:
    text = GTF_GENE + "".join(gtf_parts(f"g.b{i:02d}", [(1, 100), (50, 150)]) for i in range(12))
    proc = cap(50, text, tmp_path)
    assert proc.stdout == GTF_GENE
    err = proc.stderr.strip()
    assert "dropped 12 malformed transcript(s)" in err
    assert err.endswith(" ...")
    assert "g.b09" in err and "g.b10" not in err
