"""Tests of bin/add_utrs_from_stringtie.py on toy GFF3 and StringTie GTF inputs."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "bin"
sys.path.insert(0, str(BIN))
from gff3_lib import read_gff3  # noqa: E402


def run(script: str, *args: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(BIN / script), *args], capture_output=True, text=True, cwd=cwd)


def gene(gid: str, strand: str, cds: list[tuple[int, int]], utr5: list[tuple[int, int]] = (),
         utr3: list[tuple[int, int]] = (), seqid: str = "chr1") -> str:
    """GFF3 of a gene with one mRNA <gid>.t1; exons = union of CDS and UTRs."""
    tx = f"{gid}.t1"
    parts = sorted(list(cds) + list(utr5) + list(utr3))
    exons: list[list[int]] = []
    for s, e in parts:
        if exons and s <= exons[-1][1] + 1:
            exons[-1][1] = max(exons[-1][1], e)
        else:
            exons.append([s, e])
    lo, hi = exons[0][0], exons[-1][1]
    lines = [f"{seqid}\tTiberius\tgene\t{lo}\t{hi}\t.\t{strand}\t.\tID={gid}",
             f"{seqid}\tTiberius\tmRNA\t{lo}\t{hi}\t.\t{strand}\t.\tID={tx};Parent={gid}"]
    for n, (s, e) in enumerate(exons, 1):
        lines.append(f"{seqid}\tTiberius\texon\t{s}\t{e}\t.\t{strand}\t.\tID={tx}.exon{n};Parent={tx}")
    phase = 0
    for s, e in (cds if strand == "+" else sorted(cds, reverse=True)):
        lines.append(f"{seqid}\tTiberius\tCDS\t{s}\t{e}\t.\t{strand}\t{phase}\tID={tx}.cds;Parent={tx}")
        phase = (3 - ((e - s + 1 - phase) % 3)) % 3
    for kind, tag, utrs in (("five_prime_UTR", "utr5p", utr5), ("three_prime_UTR", "utr3p", utr3)):
        for n, (s, e) in enumerate(utrs, 1):
            lines.append(f"{seqid}\tTiberius\t{kind}\t{s}\t{e}\t.\t{strand}\t.\tID={tx}.{tag}{n};Parent={tx}")
    return "\n".join(lines) + "\n"


def stringtie(tid: str, strand: str, exons: list[tuple[int, int]], seqid: str = "chr1",
              transcript_line: bool = True) -> str:
    attrs = f'gene_id "{tid.rsplit(".", 1)[0]}"; transcript_id "{tid}";'
    lines = []
    if transcript_line:
        lines.append(f"{seqid}\tStringTie\ttranscript\t{exons[0][0]}\t{exons[-1][1]}\t1000\t{strand}\t.\t"
                     f'{attrs} cov "12.3"; FPKM "4.5"; TPM "6.7";')
    for n, (s, e) in enumerate(exons, 1):
        lines.append(f'{seqid}\tStringTie\texon\t{s}\t{e}\t1000\t{strand}\t.\t{attrs} exon_number "{n}"; cov "12.3";')
    return "\n".join(lines) + "\n"


class Result:
    def __init__(self, tmp: Path, proc: subprocess.CompletedProcess):
        assert proc.returncode == 0, proc.stderr
        self.proc = proc
        self.path = tmp / "out.gff3"
        self.ann = read_gff3(str(self.path))
        self.tx = {t.id: t for _g, t in self.ann.transcripts()}
        self.genes = {g.id: g for g in self.ann.genes}
        text = (tmp / "report.tsv").read_text()
        self.comments = [line for line in text.splitlines() if line.startswith("#")]
        body = [line.split("\t") for line in text.splitlines() if line and not line.startswith("#")]
        assert body[0] == ["transcript_id", "matched_by", "five_prime_utr_bp", "three_prime_utr_bp"]
        self.report = {r[0]: (r[1], int(r[2]), int(r[3])) for r in body[1:]}

    def spans(self, tx: str, kind: str) -> list[tuple[int, int]]:
        return [(f.start, f.end) for f in self.tx[tx].of_type(kind)]


def annotate(tmp: Path, gff: str, short: list[str] = (), long: list[str] = (), ext: int = 5000) -> Result:
    (tmp / "in.gff3").write_text("##gff-version 3\n" + gff)
    args = ["--gff3", "in.gff3", "--out", "out.gff3", "--report", "report.tsv", "--max-utr-extension", str(ext)]
    if short:
        for n, text in enumerate(short):
            (tmp / f"short{n}.gtf").write_text(text)
        args += ["--stringtie"] + [f"short{n}.gtf" for n in range(len(short))]
    if long:
        for n, text in enumerate(long):
            (tmp / f"long{n}.gtf").write_text(text)
        args += ["--longread"] + [f"long{n}.gtf" for n in range(len(long))]
    return Result(tmp, run("add_utrs_from_stringtie.py", *args, cwd=tmp))


def cds_lines(path: Path) -> list[str]:
    return sorted(line for line in path.read_text().splitlines()
                  if not line.startswith("#") and line.split("\t")[2] == "CDS")


MULTI = [(1000, 1200), (1500, 1700), (2000, 2100)]


# ---------------------------------------------------------------- matching

def test_multi_exon_gene_matched_by_intron_subset(tmp_path: Path):
    # StringTie has an extra 5' exon and intron: the CDS introns are a subset
    st = stringtie("MSTRG.1.1", "+", [(300, 500), (800, 1200), (1500, 1700), (2000, 2400)])
    res = annotate(tmp_path, gene("g1", "+", MULTI), short=[st])
    assert res.spans("g1.t1", "five_prime_UTR") == [(300, 500), (800, 999)]
    assert res.spans("g1.t1", "three_prime_UTR") == [(2101, 2400)]
    assert res.spans("g1.t1", "exon") == [(300, 500), (800, 1200), (1500, 1700), (2000, 2400)]
    assert res.spans("g1.t1", "CDS") == MULTI
    assert (res.tx["g1.t1"].feature.start, res.tx["g1.t1"].feature.end) == (300, 2400)
    assert (res.genes["g1"].feature.start, res.genes["g1"].feature.end) == (300, 2400)
    assert res.report["g1.t1"] == ("short", 201 + 200, 300)
    utr = res.tx["g1.t1"].of_type("five_prime_UTR")[0]
    assert utr.source == "stringtie2utr" and utr.id == "g1.t1.utr5p1"
    assert [e.id for e in res.tx["g1.t1"].exons] == [f"g1.t1.exon{n}" for n in range(1, 5)]


def test_stringtie_transcript_missing_an_intron_does_not_match(tmp_path: Path):
    st = stringtie("MSTRG.1.1", "+", [(800, 1200), (1500, 2400)])      # intron 1701-1999 retained
    res = annotate(tmp_path, gene("g1", "+", MULTI), short=[st])
    assert res.report["g1.t1"] == ("none", 0, 0)
    assert res.tx["g1.t1"].utrs == []
    assert res.spans("g1.t1", "exon") == MULTI


def test_other_strand_and_unstranded_transcripts_do_not_match(tmp_path: Path):
    exons = [(800, 1200), (1500, 1700), (2000, 2400)]
    st = stringtie("MSTRG.1.1", "-", exons) + stringtie("MSTRG.2.1", ".", exons)
    res = annotate(tmp_path, gene("g1", "+", MULTI), short=[st])
    assert res.report["g1.t1"] == ("none", 0, 0)


def test_all_stringtie_transcripts_per_intron_longest_wins(tmp_path: Path):
    # both StringTie transcripts and both genes share intron 1201-1499; the
    # original kept only one transcript per intron
    st = (stringtie("MSTRG.1.1", "+", [(900, 1200), (1500, 1700), (2000, 2600)])
          + stringtie("MSTRG.1.2", "+", [(950, 1200), (1500, 1700), (2000, 2200)]))
    gff = gene("g1", "+", MULTI) + gene("g2", "+", [(1000, 1200), (1500, 1700)])
    res = annotate(tmp_path, gff, short=[st])
    assert res.spans("g1.t1", "five_prime_UTR") == [(900, 999)]
    assert res.spans("g1.t1", "three_prime_UTR") == [(2101, 2600)]
    assert res.report["g2.t1"] == ("short", 100, 601)


def test_single_exon_gene_matched_by_overlap(tmp_path: Path):
    st = stringtie("MSTRG.1.1", "+", [(2800, 3800)], transcript_line=False)
    res = annotate(tmp_path, gene("g1", "+", [(3000, 3600)]), short=[st])
    assert res.spans("g1.t1", "five_prime_UTR") == [(2800, 2999)]
    assert res.spans("g1.t1", "three_prime_UTR") == [(3601, 3800)]
    assert res.spans("g1.t1", "exon") == [(2800, 3800)]


def test_single_exon_gene_rejected_by_intron_in_cds(tmp_path: Path):
    st = stringtie("MSTRG.1.1", "+", [(2800, 3100), (3300, 3800)])
    res = annotate(tmp_path, gene("g1", "+", [(3000, 3600)]), short=[st])
    assert res.report["g1.t1"] == ("none", 0, 0)
    assert res.spans("g1.t1", "exon") == [(3000, 3600)]


def test_single_exon_stringtie_exon_starting_inside_cds_gives_no_utr_over_the_cds(tmp_path: Path):
    st = stringtie("MSTRG.1.1", "+", [(3200, 4000)])     # longer than the CDS, starts inside it
    res = annotate(tmp_path, gene("g1", "+", [(3000, 3600)]), short=[st])
    assert res.spans("g1.t1", "five_prime_UTR") == []
    assert res.spans("g1.t1", "three_prime_UTR") == [(3601, 4000)]
    assert res.spans("g1.t1", "exon") == [(3000, 4000)]


def test_stringtie_exon_shorter_than_overlapped_cds_segment_is_not_used(tmp_path: Path):
    # the StringTie transcript starts inside the first CDS segment: no 5' UTR
    st = stringtie("MSTRG.1.1", "+", [(1100, 1200), (1500, 1700), (2000, 2400)])
    res = annotate(tmp_path, gene("g1", "+", MULTI), short=[st])
    assert res.report["g1.t1"] == ("short", 0, 300)


def test_utr_beyond_a_dropped_stringtie_exon_is_not_used(tmp_path: Path):
    # the StringTie exon over the first CDS segment is shorter than it and dropped; the exon
    # further upstream would be joined to the CDS by an intron (701-999) StringTie does not have
    st = stringtie("MSTRG.1.1", "+", [(500, 700), (1100, 1200), (1500, 1700), (2000, 2400)])
    res = annotate(tmp_path, gene("g1", "+", MULTI), short=[st])
    assert res.report["g1.t1"] == ("short", 0, 300)
    assert res.spans("g1.t1", "five_prime_UTR") == []
    assert res.spans("g1.t1", "exon") == [(1000, 1200), (1500, 1700), (2000, 2400)]
    # same on the other side, minus strand: the 5' side of the transcript is the CDS end
    st = stringtie("MSTRG.2.1", "-", [(500, 900), (1000, 1200), (1500, 1700), (2000, 2050), (2300, 2400)])
    res = annotate(tmp_path, gene("g2", "-", MULTI), short=[st])
    assert res.report["g2.t1"] == ("short", 0, 401)
    assert res.spans("g2.t1", "three_prime_UTR") == [(500, 900)]
    assert res.spans("g2.t1", "exon") == [(500, 900), (1000, 1200), (1500, 1700), (2000, 2100)]


# ---------------------------------------------------------------- clipping

def test_same_strand_matched_neighbour_clips_the_extension(tmp_path: Path):
    # a read-through transcript joins g1 and g2; both match it
    st = stringtie("MSTRG.1.1", "+", [(900, 2600)])
    gff = gene("g1", "+", [(1000, 1500)]) + gene("g2", "+", [(2000, 2500)])
    res = annotate(tmp_path, gff, short=[st])
    assert res.spans("g1.t1", "three_prime_UTR") == [(1501, 1999)]    # stops before g2
    assert res.spans("g2.t1", "five_prime_UTR") == [(1501, 1999)]     # starts after g1
    assert res.spans("g2.t1", "three_prime_UTR") == [(2501, 2600)]
    assert any("stopped at a neighbour gene: 2" in line for line in res.comments)


def test_unmatched_neighbour_does_not_clip(tmp_path: Path):
    # g2 has no match (StringTie lacks its intron), so it does not stop g1's UTR
    st = stringtie("MSTRG.1.1", "+", [(900, 2600)])
    gff = gene("g1", "+", [(1000, 1500)]) + gene("g2", "+", [(2000, 2100), (2300, 2400)])
    res = annotate(tmp_path, gff, short=[st])
    assert res.report["g2.t1"][0] == "none"
    assert res.spans("g1.t1", "three_prime_UTR") == [(1501, 2600)]


def test_gene_with_utrs_clips_the_extension_on_both_sides(tmp_path: Path):
    # the HC genes g1 and g3 already have UTRs: they are no candidates, but they are barriers
    st = stringtie("MSTRG.1.1", "+", [(200, 2400)])
    gff = (gene("g1", "+", [(100, 400)], utr3=[(401, 500)])
           + gene("g2", "+", [(1000, 1500)])
           + gene("g3", "+", [(2100, 2500)], utr5=[(2000, 2099)]))
    res = annotate(tmp_path, gff, short=[st])
    assert res.report["g1.t1"][0] == "has_utr" and res.report["g3.t1"][0] == "has_utr"
    assert res.spans("g2.t1", "five_prime_UTR") == [(501, 999)]        # starts after g1's UTR
    assert res.spans("g2.t1", "three_prime_UTR") == [(1501, 1999)]     # stops before g3's UTR
    assert res.spans("g2.t1", "exon") == [(501, 1999)]
    assert any("stopped at a neighbour gene: 1" in line for line in res.comments)


def test_max_utr_extension_clips(tmp_path: Path):
    st = stringtie("MSTRG.1.1", "+", [(1000, 9000)])
    res = annotate(tmp_path, gene("g1", "+", [(5000, 5300)]), short=[st], ext=500)
    assert res.spans("g1.t1", "five_prime_UTR") == [(4500, 4999)]
    assert res.spans("g1.t1", "three_prime_UTR") == [(5301, 5800)]


def test_max_utr_extension_drops_distant_exons(tmp_path: Path):
    st = stringtie("MSTRG.1.1", "+", [(100, 200), (800, 1200), (1500, 1700), (2000, 2100), (2500, 2600)])
    res = annotate(tmp_path, gene("g1", "+", MULTI), short=[st], ext=300)
    assert res.spans("g1.t1", "five_prime_UTR") == [(800, 999)]
    assert res.spans("g1.t1", "three_prime_UTR") == []


# ---------------------------------------------------------------- long reads, strand, existing UTRs

def test_long_read_match_overrides_longer_short_read_match(tmp_path: Path):
    short = stringtie("MSTRG.1.1", "+", [(500, 1200), (1500, 1700), (2000, 3000)])
    long = stringtie("MSTRG.1.1", "+", [(900, 1200), (1500, 1700), (2000, 2200)])
    res = annotate(tmp_path, gene("g1", "+", MULTI), short=[short], long=[long])
    assert res.report["g1.t1"] == ("long", 100, 100)
    assert res.spans("g1.t1", "exon") == [(900, 1200), (1500, 1700), (2000, 2200)]


def test_several_short_read_files_longest_wins(tmp_path: Path):
    a = stringtie("MSTRG.1.1", "+", [(900, 1200), (1500, 1700), (2000, 2200)])
    b = stringtie("MSTRG.1.1", "+", [(700, 1200), (1500, 1700), (2000, 2200)])
    res = annotate(tmp_path, gene("g1", "+", MULTI), short=[a, b, ""])     # the empty file is fine
    assert res.report["g1.t1"] == ("short", 300, 100)


def test_minus_strand_utr_types(tmp_path: Path):
    st = stringtie("MSTRG.1.1", "-", [(900, 1200), (1500, 1900)])
    res = annotate(tmp_path, gene("g1", "-", [(1000, 1200), (1500, 1700)]), short=[st])
    assert res.spans("g1.t1", "five_prime_UTR") == [(1701, 1900)]
    assert res.spans("g1.t1", "three_prime_UTR") == [(900, 999)]
    assert res.report["g1.t1"] == ("short", 200, 100)


def test_transcript_with_utrs_is_unchanged(tmp_path: Path):
    gff = gene("g1", "+", MULTI, utr5=[(950, 999)], utr3=[(2101, 2150)])
    st = stringtie("MSTRG.1.1", "+", [(500, 1200), (1500, 1700), (2000, 3000)])
    res = annotate(tmp_path, gff, short=[st])
    assert res.report["g1.t1"] == ("has_utr", 50, 50)
    want = sorted(line for line in gff.splitlines())
    got = sorted(line for line in res.path.read_text().splitlines() if not line.startswith("#"))
    assert got == want


def test_stringtie_ending_at_the_stop_codon_gives_no_three_prime_utr(tmp_path: Path):
    # the CDS includes the stop codon (1700-1702); StringTie ends exactly there
    cds = [(1000, 1200), (1500, 1702)]
    st = stringtie("MSTRG.1.1", "+", [(800, 1200), (1500, 1702)])
    res = annotate(tmp_path, gene("g1", "+", cds), short=[st])
    assert res.spans("g1.t1", "three_prime_UTR") == []
    assert res.spans("g1.t1", "five_prime_UTR") == [(800, 999)]
    assert all(f.length != 3 for f in res.tx["g1.t1"].children if f.type != "CDS")
    assert res.report["g1.t1"] == ("short", 200, 0)


def test_no_stringtie_input_writes_the_input(tmp_path: Path):
    gff = gene("g1", "+", MULTI) + gene("g2", "-", [(5000, 5600)])
    res = annotate(tmp_path, gff)
    assert res.report == {"g1.t1": ("none", 0, 0), "g2.t1": ("none", 0, 0)}
    got = sorted(line for line in res.path.read_text().splitlines() if not line.startswith("#"))
    assert got == sorted(gff.splitlines())


# ---------------------------------------------------------------- invariants and validity

def mixed_case() -> tuple[str, str, str]:
    gff = (gene("g1", "+", MULTI)
           + gene("g2", "-", [(5000, 5200), (5500, 5700)])
           + gene("g3", "+", [(8000, 8600)])
           + gene("g4", "+", [(10000, 10300), (10600, 10900)], utr5=[(9900, 9999)])
           + gene("g5", "+", [(20000, 20300)], seqid="chr2"))
    short = (stringtie("MSTRG.1.1", "+", [(300, 500), (800, 1200), (1500, 1700), (2000, 2400)])
             + stringtie("MSTRG.2.1", "-", [(4800, 5200), (5500, 6000)])
             + stringtie("MSTRG.3.1", "+", [(7900, 8700)])
             + stringtie("MSTRG.4.1", "+", [(9000, 10300), (10600, 11500)])
             + stringtie("MSTRG.1.1", "+", [(19000, 21000)], seqid="chr2"))   # ID reused on another sequence
    long = stringtie("PB.1.1", "+", [(7950, 8650)])
    return gff, short, long


def test_cds_unchanged(tmp_path: Path):
    gff, short, long = mixed_case()
    res = annotate(tmp_path, gff, short=[short], long=[long])
    assert cds_lines(tmp_path / "in.gff3") == cds_lines(res.path)
    assert res.report["g3.t1"] == ("long", 50, 50)
    assert res.report["g4.t1"][0] == "has_utr"
    assert res.report["g5.t1"] == ("short", 1000, 700)
    assert any(line == "# matched long: 1" for line in res.comments)
    assert "add_utrs_from_stringtie.py:" in res.proc.stderr and len(res.proc.stderr.splitlines()) == 1


def gt_binary() -> str | None:
    return os.environ.get("GT_BIN") or shutil.which("gt")


def test_output_passes_normalize_and_gt(tmp_path: Path):
    gff, short, long = mixed_case()
    res = annotate(tmp_path, gff, short=[short], long=[long])
    proc = run("normalize_gff3.py", "--gff3", str(res.path), "--out", "norm.gff3", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert cds_lines(tmp_path / "in.gff3") != []
    norm = read_gff3(str(tmp_path / "norm.gff3"))
    assert sum(1 for _g, t in norm.transcripts() if t.utrs) == 5
    gt = gt_binary()
    if not gt:
        pytest.skip("GenomeTools gt not available (set GT_BIN)")
    proc = subprocess.run([gt, "gff3validator", str(tmp_path / "norm.gff3")], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
