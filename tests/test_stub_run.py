"""Nextflow stub runs: check the wiring and the published file names of every mode.

No tools, containers or GPU are needed, only a ``nextflow`` executable (or the
environment variable ``NEXTFLOW_BIN``). The tests are skipped without it.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "tests" / "data"
NEXTFLOW = os.environ.get("NEXTFLOW_BIN") or shutil.which("nextflow")

pytestmark = [
    pytest.mark.nextflow,
    pytest.mark.skipif(NEXTFLOW is None, reason="nextflow executable not found"),
]

EVIDENCE_OUTPUTS = {
    "tiberius_evidence.gff3",
    "tiberius_evidence_proteins.fa",
    "intermediate/tiberius_ab_initio.gff3",
    "intermediate/hc.gff3",
    "hintsfile.gff",
}

PROTEINS = {"proteins": str(DATA / "tiny_proteins.faa")}
SHORT_READS = {"rnaseq_paired": str(DATA / "reads_{1,2}.fastq")}
ISOSEQ = {"isoseq": [str(DATA / "isoseq.fastq")]}

MODES = {
    "abinitio": ({}, {"tiberius_ab_initio.gff3"}),
    "proteins": (PROTEINS, EVIDENCE_OUTPUTS),
    "rnaseq": ({**PROTEINS, **SHORT_READS}, EVIDENCE_OUTPUTS),
    "isoseq": ({**PROTEINS, **ISOSEQ}, EVIDENCE_OUTPUTS),
    "mixed": ({**PROTEINS, **SHORT_READS, **ISOSEQ}, EVIDENCE_OUTPUTS),
}


def published_files(outdir: Path) -> set[str]:
    return {str(p.relative_to(outdir)) for p in outdir.rglob("*") if p.is_file()}


@pytest.mark.parametrize("mode", sorted(MODES))
def test_stub_run_publishes_expected_files(mode: str, tmp_path: Path) -> None:
    evidence, expected = MODES[mode]
    outdir = tmp_path / "out"
    params = {
        "genome": str(DATA / "tiny.fa"),
        "outdir": str(outdir),
        # any existing file works as model config in a stub run
        "tiberius": {"run": True, "model_cfg": str(DATA / "tiny.fa")},
        **evidence,
    }
    params_file = tmp_path / "params.yaml"
    params_file.write_text(yaml.safe_dump(params))

    env = dict(os.environ, NXF_ANSI_LOG="false")
    proc = subprocess.run(
        [
            NEXTFLOW, "run", str(ROOT / "main.nf"), "-stub-run",
            "-params-file", str(params_file),
            "-c", str(ROOT / "tests" / "stub.config"),
            "-work-dir", str(tmp_path / "work"),
        ],
        cwd=tmp_path, env=env, capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stdout[-4000:] + proc.stderr[-2000:]

    published = published_files(outdir)
    # params.yaml etc. may be added by launchers; only require the expected set
    missing = expected - published
    assert not missing, f"missing {missing}; published {sorted(published)}"
    if mode == "abinitio":
        assert "tiberius_evidence.gff3" not in published
