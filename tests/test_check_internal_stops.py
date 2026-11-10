"""bin/check_internal_stops.py on toy protein FASTA files: a warning on stdout
when more than 5% of the proteins have a stop before their last residue,
silence otherwise, exit status 0 on valid usage."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "bin" / "check_internal_stops.py"


def check(tmp_path: Path, sequences: list[str], name: str = "tiberius", raw: str | None = None):
    """Run the script on a FASTA of one record per element of ``sequences``
    (a sequence may contain newlines), or on the text ``raw``."""
    fasta = tmp_path / "proteins.fa"
    fasta.write_text(raw if raw is not None
                     else "".join(f">p{i}\n{seq}\n" for i, seq in enumerate(sequences)))
    return subprocess.run([sys.executable, str(SCRIPT), str(fasta), name],
                          capture_output=True, text=True)


def test_above_threshold_warns_with_percent_counts_and_name(tmp_path: Path) -> None:
    # 2 of 10 = 20.0%
    proc = check(tmp_path, ["MKV*LL"] * 2 + ["MKVLL*"] * 8, name="drusilla")
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.startswith("WARN:")
    assert "20.0%" in proc.stdout
    assert "(2 of 10)" in proc.stdout
    assert "drusilla" in proc.stdout


def test_exactly_five_percent_is_silent(tmp_path: Path) -> None:
    proc = check(tmp_path, ["MKV*LL"] + ["MKVLL"] * 19)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout == ""


def test_just_above_five_percent_warns(tmp_path: Path) -> None:
    # 2 of 39 = 5.1%
    proc = check(tmp_path, ["MKV*LL"] * 2 + ["MKVLL"] * 37)
    assert proc.returncode == 0, proc.stderr
    assert "5.1%" in proc.stdout and "(2 of 39)" in proc.stdout


@pytest.mark.parametrize("stop", ["*", "."])
def test_terminal_stop_is_not_counted(stop: str, tmp_path: Path) -> None:
    proc = check(tmp_path, [f"MKVLL{stop}"] * 10)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout == ""


def test_internal_dot_is_counted(tmp_path: Path) -> None:
    proc = check(tmp_path, ["MKV.LL"] * 10)
    assert proc.returncode == 0, proc.stderr
    assert "100.0%" in proc.stdout and "(10 of 10)" in proc.stdout


def test_stop_in_multiline_sequence_is_counted(tmp_path: Path) -> None:
    proc = check(tmp_path, ["MKVLL\nMK*LL\nMKVLL*"] * 10)
    assert proc.returncode == 0, proc.stderr
    assert "(10 of 10)" in proc.stdout


def test_stop_at_end_of_a_line_that_is_not_the_last_is_counted(tmp_path: Path) -> None:
    proc = check(tmp_path, ["MKVLL*\nMKVLL*"] * 10)
    assert proc.returncode == 0, proc.stderr
    assert "(10 of 10)" in proc.stdout


def test_terminal_stop_on_last_line_of_multiline_sequence_is_not_counted(tmp_path: Path) -> None:
    proc = check(tmp_path, ["MKVLL\nMKVLL*"] * 10)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout == ""


def test_empty_file_is_silent(tmp_path: Path) -> None:
    proc = check(tmp_path, [], raw="")
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout == ""


@pytest.mark.parametrize("argv", [[], ["only_one.fa"], ["a.fa", "name", "extra"]])
def test_wrong_argument_count_exits_non_zero(argv: list[str]) -> None:
    proc = subprocess.run([sys.executable, str(SCRIPT), *argv], capture_output=True, text=True)
    assert proc.returncode != 0
    assert "Usage" in proc.stderr
