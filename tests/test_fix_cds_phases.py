"""Tests for bin/fix_cds_phases.py: CDS phases of transcripts with a start
codon, counted from the 5' end in transcription order."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parent.parent / "bin"


def run(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(BIN / "fix_cds_phases.py"), *args], capture_output=True, text=True)


def rows(text: str) -> list[list[str]]:
    return [line.split("\t") for line in text.splitlines() if line and not line.startswith("#")]


def test_fix_cds_phases_counts_minus_strand_phases_from_the_start_codon(tmp_path):
    # rescued_t649 of the T. rubripes hint rescue, three of its CDS: lengths 4,
    # 5, 103 in transcription order; the hint rescue Tiberius wrote 1, 0, 1
    a = 'gene_id "g1"; transcript_id "t1";'
    gtf = [
        f"c\tTiberius\ttranscript\t2875770\t2877481\t.\t-\t.\t{a}",
        f"c\tTiberius\tstart_codon\t2877479\t2877481\t.\t-\t0\t{a}",
        f"c\tTiberius\tCDS\t2875770\t2875872\t.\t-\t1\t{a}",
        f"c\tTiberius\tCDS\t2876379\t2876383\t.\t-\t0\t{a}",
        f"c\tTiberius\tCDS\t2877478\t2877481\t.\t-\t1\t{a}",
        # no start codon: phases are kept
        'c\tTiberius\tCDS\t100\t200\t.\t-\t2\tgene_id "g2"; transcript_id "t2";',
    ]
    (tmp_path / "in.gtf").write_text("\n".join(gtf) + "\n")
    r = run(str(tmp_path / "in.gtf"))
    assert r.returncode == 0, r.stderr
    phases = {(row[3], row[8].split('"')[3]): row[7] for row in rows(r.stdout) if row[2] == "CDS"}
    assert phases == {("2877478", "t1"): "0", ("2876379", "t1"): "2", ("2875770", "t1"): "0",
                      ("100", "t2"): "2"}
    assert "rewritten in 1 of 1" in r.stderr
