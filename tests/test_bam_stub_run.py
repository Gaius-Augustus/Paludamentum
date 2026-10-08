"""Stub runs of the BAM handling of the evidence subworkflows.

A single BAM is passed on as it is (no SAMTOOLS_MERGE), and only user BAMs
(rnaseq_bam) whose header does not say SO:coordinate are sorted. bam2hints and
StringTie read the BAMs without sorting them again.

Needs a ``nextflow`` executable (or NEXTFLOW_BIN), like test_stub_run.py.
"""
from __future__ import annotations

import struct
import zlib
from pathlib import Path

import pytest

from test_stub_run import (GENEFINDER, ISOSEQ, NEXTFLOW, PROTEINS, SHORT_READS, assert_ok,
                           drusilla_params, run_pipeline)

pytestmark = [
    pytest.mark.nextflow,
    pytest.mark.skipif(NEXTFLOW is None, reason="nextflow executable not found"),
]

BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def bgzf_block(data: bytes) -> bytes:
    """One BGZF block: a gzip member with the BC extra field (block size - 1)."""
    deflate = zlib.compressobj(6, zlib.DEFLATED, -15)
    payload = deflate.compress(data) + deflate.flush()
    header = b"\x1f\x8b\x08\x04" + b"\x00" * 4 + b"\x00\xff" + struct.pack("<H", 6)
    extra = b"BC" + struct.pack("<HH", 2, len(payload) + 25)
    return header + extra + payload + struct.pack("<II", zlib.crc32(data), len(data))


def write_bam(path: Path, sort_order: str | None) -> str:
    """A BAM header without reads; @HD with SO:<sort_order>, or no @HD line."""
    text = (f"@HD\tVN:1.6\tSO:{sort_order}\n" if sort_order else "") + "@SQ\tSN:seq1\tLN:100\n"
    name = b"seq1\x00"
    body = (b"BAM\x01" + struct.pack("<i", len(text)) + text.encode() + struct.pack("<i", 1)
            + struct.pack("<i", len(name)) + name + struct.pack("<i", 100))
    path.write_bytes(bgzf_block(body) + BGZF_EOF)
    return str(path)


def tasks(stdout: str, process: str) -> int:
    """Number of tasks of a process (its last name component) in the log."""
    return sum(1 for line in stdout.splitlines()
               if "process >" in line and line.split("process >")[1].split()[0].split(":")[-1] == process)


@pytest.mark.parametrize("sort_order,sorted_by_pipeline", [
    ("coordinate", False), ("unsorted", True), ("queryname", True), (None, True),
], ids=["coordinate", "unsorted", "queryname", "no-hd-line"])
def test_one_user_bam_is_sorted_by_its_header_and_not_merged(sort_order, sorted_by_pipeline,
                                                             tmp_path: Path) -> None:
    bam = write_bam(tmp_path / "lib.bam", sort_order)
    proc, published = run_pipeline(tmp_path, {**GENEFINDER["tiberius"], **PROTEINS, "rnaseq_bam": [bam]})
    assert_ok(proc)
    assert tasks(proc.stdout, "SAMTOOLS_SORT") == int(sorted_by_pipeline), proc.stdout
    assert tasks(proc.stdout, "SAMTOOLS_MERGE") == 0, proc.stdout
    assert tasks(proc.stdout, "BAM2HINTS_RNA") == 1, proc.stdout
    assert "hintsfile.gff" in published


def test_an_empty_file_counts_as_unsorted(tmp_path: Path) -> None:
    bam = tmp_path / "empty.bam"
    bam.write_bytes(b"")
    proc, _ = run_pipeline(tmp_path, {**GENEFINDER["tiberius"], **PROTEINS, "rnaseq_bam": [str(bam)]})
    assert_ok(proc)
    assert tasks(proc.stdout, "SAMTOOLS_SORT") == 1, proc.stdout


def test_two_user_bams_only_the_unsorted_one_is_sorted(tmp_path: Path) -> None:
    bams = [write_bam(tmp_path / "a.bam", "coordinate"), write_bam(tmp_path / "b.bam", "unsorted")]
    proc, _ = run_pipeline(tmp_path, {**GENEFINDER["tiberius"], **PROTEINS, "rnaseq_bam": bams})
    assert_ok(proc)
    assert tasks(proc.stdout, "SAMTOOLS_SORT") == 1, proc.stdout
    assert tasks(proc.stdout, "SAMTOOLS_MERGE") == 1, proc.stdout


def test_bam_and_reads_are_merged(tmp_path: Path) -> None:
    params = {**GENEFINDER["tiberius"], **PROTEINS, **SHORT_READS,
              "rnaseq_bam": [write_bam(tmp_path / "a.bam", "coordinate")]}
    proc, _ = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert tasks(proc.stdout, "SAMTOOLS_SORT") == 0, proc.stdout
    assert tasks(proc.stdout, "SAMTOOLS_MERGE") == 1, proc.stdout


@pytest.mark.parametrize("flow", ["transdecoder", "drusilla"])
def test_one_library_per_read_type_is_not_merged(flow: str, tmp_path: Path) -> None:
    """Mixed mode with one short-read and one Iso-Seq library; the Drusilla flow
    assembles both BAMs as they are (STRINGTIE_ASSEMBLE_MIX)."""
    base = drusilla_params(tmp_path, "tiberius") if flow == "drusilla" else GENEFINDER["tiberius"]
    proc, published = run_pipeline(tmp_path, {**base, **PROTEINS, **SHORT_READS, **ISOSEQ})
    assert_ok(proc)
    assert tasks(proc.stdout, "SAMTOOLS_MERGE") == 0, proc.stdout
    assert tasks(proc.stdout, "SAMTOOLS_SORT") == 0, proc.stdout
    assert tasks(proc.stdout, "BAM2HINTS_RNA") == 1 and tasks(proc.stdout, "BAM2HINTS_ISO") == 1, proc.stdout
    assert (tasks(proc.stdout, "STRINGTIE_ASSEMBLE_MIX") == 1) == (flow == "drusilla"), proc.stdout
    assert "tiberius_evidence.gff3" in published
