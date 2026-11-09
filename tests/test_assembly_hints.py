"""bin/assembly_hints.py: intron hints and a coverage track from a StringTie assembly (a `stringtie` input
without reads), read back with the parsers of combine_genesets.py."""
from __future__ import annotations

import gzip
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "bin"
SCRIPT = BIN / "assembly_hints.py"

sys.path.insert(0, str(BIN))
import combine_genesets as cg  # noqa: E402

GTF = """\
# stringtie -p 4 -o stringtie.gtf reads.bam
chr1\tStringTie\ttranscript\t101\t400\t1000\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.1"; cov "10.0"; FPKM "1"; TPM "2";
chr1\tStringTie\texon\t101\t200\t1000\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1"; cov "10.0";
chr1\tStringTie\texon\t301\t400\t1000\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "2"; cov "10.0";
chr1\tStringTie\ttranscript\t151\t600\t1000\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.2"; cov "4.4"; FPKM "1"; TPM "1";
chr1\tStringTie\texon\t151\t200\t1000\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.2"; exon_number "1"; cov "4.4";
chr1\tStringTie\texon\t301\t350\t1000\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.2"; exon_number "2"; cov "4.4";
chr1\tStringTie\texon\t501\t600\t1000\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.2"; exon_number "3"; cov "4.4";
chr1\tStringTie\ttranscript\t1001\t1100\t1000\t-\t.\tgene_id "STRG.2"; transcript_id "STRG.2.1"; cov "2.5"; FPKM "1"; TPM "1";
chr1\tStringTie\texon\t1001\t1100\t1000\t-\t.\tgene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "1"; cov "2.5";
chr2\tStringTie\ttranscript\t11\t80\t1000\t-\t.\tgene_id "STRG.3"; transcript_id "STRG.3.1"; cov "1.0"; FPKM "1"; TPM "1";
chr2\tStringTie\texon\t11\t20\t1000\t-\t.\tgene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "1"; cov "1.0";
chr2\tStringTie\texon\t41\t80\t1000\t-\t.\tgene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "2"; cov "1.0";
"""


def run(tmp_path: Path, *extra: str) -> tuple[Path, Path, subprocess.CompletedProcess]:
    gtf = tmp_path / "stringtie.gtf"
    gtf.write_text(GTF)
    hints, bg = tmp_path / "introns.gff", tmp_path / "cov.bg.gz"
    r = subprocess.run([sys.executable, str(SCRIPT), str(gtf), "--introns", str(hints), "--bedgraph", str(bg), *extra],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return hints, bg, r


def test_introns_with_summed_cov_as_multiplicity(tmp_path: Path):
    hints, _, r = run(tmp_path)
    assert "4 transcripts, 3 distinct introns, mult = cov" in r.stderr
    lines = [l.split("\t") for l in hints.read_text().splitlines()]
    # bam2hints layout: first to last intron base, score = mult, mult= attribute, src=E
    assert [(l[0], l[2], l[3], l[4], l[5], l[6]) for l in lines] == [
        ("chr1", "intron", "201", "300", "14", "+"),      # shared by STRG.1.1 (10.0) and STRG.1.2 (4.4): 14.4 -> 14
        ("chr1", "intron", "351", "500", "4", "+"),       # STRG.1.2 alone: 4.4 -> 4
        ("chr2", "intron", "21", "40", "1", "-"),
    ]
    assert all(l[1] == "asm2h" and l[8] == f"mult={l[5]};pri=4;src=E" for l in lines)
    # the single-exon transcript gives no intron; the parser of combine_genesets.py reads the multiplicities
    mult = cg.read_intron_hints([str(hints)])
    assert mult == {("chr1", 201, 300): 14, ("chr1", 351, 500): 4, ("chr2", 21, 40): 1}


def test_introns_counted_per_transcript(tmp_path: Path):
    hints, _, r = run(tmp_path, "--mult", "count")
    assert "mult = count" in r.stderr
    mult = cg.read_intron_hints([str(hints)])
    assert mult == {("chr1", 201, 300): 2, ("chr1", 351, 500): 1, ("chr2", 21, 40): 1}


def test_coverage_track_is_the_largest_exon_cov(tmp_path: Path):
    _, bg, r = run(tmp_path)
    assert "coverage track: 6 segments" in r.stderr
    with gzip.open(bg, "rt") as fh:
        rows = [l.split("\t") for l in fh.read().splitlines()]
    # 0-based half-open, zero gaps left out, the largest cov where exons overlap (10 over 4.4 in 150-200 and
    # 300-350), neighbouring segments of equal depth merged (100-150 + 150-200, 300-350 + 350-400)
    assert rows == [["chr1", "100", "200", "10"], ["chr1", "300", "400", "10"], ["chr1", "500", "600", "4.4"],
                    ["chr1", "1000", "1100", "2.5"], ["chr2", "10", "20", "1"], ["chr2", "40", "80", "1"]]
    cov = cg.read_bedgraph(str(bg), 2)       # depth >= 2 as combine_genesets.py --min-depth
    assert cg.covered_bases(cov, "chr1", 101, 400) == 200      # both exons of STRG.1.1, not the intron
    assert cg.covered_bases(cov, "chr1", 501, 600) == 100
    assert cg.covered_bases(cov, "chr1", 1001, 1100) == 100    # 2.5 >= 2
    assert cg.covered_bases(cov, "chr2", 11, 80) == 0          # depth 1 < 2


def test_transcripts_without_cov_count_as_one(tmp_path: Path):
    gtf = tmp_path / "merged.gtf"
    gtf.write_text("\n".join(l.replace(' cov "10.0";', "").replace(' cov "4.4";', "") for l in GTF.splitlines()) + "\n")
    hints = tmp_path / "introns.gff"
    r = subprocess.run([sys.executable, str(SCRIPT), str(gtf), "--introns", str(hints)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "2 transcripts without cov attribute" in r.stderr
    assert cg.read_intron_hints([str(hints)]) == {("chr1", 201, 300): 2, ("chr1", 351, 500): 1, ("chr2", 21, 40): 1}


def test_help_runs():
    r = subprocess.run([sys.executable, str(SCRIPT), "--help"], capture_output=True, text=True)
    assert r.returncode == 0 and "Intron hints and a coverage track" in r.stdout
