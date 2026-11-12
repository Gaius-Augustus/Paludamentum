"""bin/concat_assemblies.py: a short-read and a long-read StringTie assembly side by side in one GTF."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent.parent / "bin" / "concat_assemblies.py"

SHORT = ('# stringtie short.bam\n'
         'c1\tStringTie\ttranscript\t1\t500\t1000\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.1"; cov "5.5"; FPKM "1.0"; TPM "2.0";\n'
         'c1\tStringTie\texon\t1\t200\t1000\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1"; cov "5.0";\n'
         'c1\tStringTie\texon\t300\t500\t1000\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "2"; cov "6.0";\n')
LONG = ('c1\tStringTie\ttranscript\t1\t500\t1000\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.1"; cov "2.0"; FPKM "9.0"; TPM "8.0";\n'
        'c1\tStringTie\texon\t1\t500\t1000\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1"; cov "2.0";\n'
        '\n'
        'c2\tStringTie\ttranscript\t7\t90\t1000\t-\t.\tgene_id "STRG.2"; transcript_id "STRG.2.1"; reference_id "t"; cov "1.0";\n'
        'c2\tStringTie\texon\t7\t90\t1000\t-\t.\tgene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "1"; cov "1.0";\n')


def run(tmp_path: Path, short: str, long_: str) -> subprocess.CompletedProcess:
    (tmp_path / "short.gtf").write_text(short)
    (tmp_path / "long.gtf").write_text(long_)
    return subprocess.run([sys.executable, str(SCRIPT), "--short", str(tmp_path / "short.gtf"),
                           "--long", str(tmp_path / "long.gtf")], capture_output=True, text=True)


def test_ids_are_prefixed_and_everything_else_is_kept(tmp_path: Path) -> None:
    proc = run(tmp_path, SHORT, LONG)
    assert proc.returncode == 0, proc.stderr
    lines = [line for line in proc.stdout.splitlines() if not line.startswith("#")]
    assert len(lines) == 7                                   # 3 short + 4 long features, no blank line
    short, long_ = lines[:3], lines[3:]
    assert all('gene_id "sr.STRG.1"; transcript_id "sr.STRG.1.1";' in line for line in short)
    assert all('gene_id "lr.STRG.' in line and 'transcript_id "lr.STRG.' in line for line in long_)
    # the same STRG.1.1 of both assemblies is two transcripts now
    ids = {line.split('transcript_id "')[1].split('"')[0] for line in lines}
    assert ids == {"sr.STRG.1.1", "lr.STRG.1.1", "lr.STRG.2.1"}
    # coordinates, strand and the attributes that the Drusilla pre-filter reads are unchanged
    assert short[0].split("\t")[:8] == ["c1", "StringTie", "transcript", "1", "500", "1000", "+", "."]
    assert 'cov "5.5"; FPKM "1.0"; TPM "2.0";' in short[0]
    assert 'reference_id "t"; cov "1.0";' in long_[2]        # other attributes are not touched


def test_a_file_that_is_not_an_assembly_is_an_error(tmp_path: Path) -> None:
    proc = run(tmp_path, SHORT, "# nothing\n")
    assert proc.returncode != 0 and "no feature lines" in proc.stderr
    proc = run(tmp_path, SHORT, "c1\tx\tgene\t1\t2\t.\t+\t.\tID=g1\n")
    assert proc.returncode != 0 and "without transcript_id" in proc.stderr
    proc = run(tmp_path, "c1 StringTie transcript 1 2\n", LONG)
    assert proc.returncode != 0 and "9 columns" in proc.stderr
