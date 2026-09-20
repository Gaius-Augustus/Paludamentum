"""Nextflow stub runs: check the wiring and the published file names.

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

PROTEINS = {"proteins": str(DATA / "tiny_proteins.faa")}
SHORT_READS = {"rnaseq_paired": str(DATA / "reads_{1,2}.fastq")}
ISOSEQ = {"isoseq": [str(DATA / "isoseq.fastq")]}

EVIDENCE = {
    "abinitio": {},
    "proteins": PROTEINS,
    "rnaseq": {**PROTEINS, **SHORT_READS},
    "isoseq": {**PROTEINS, **ISOSEQ},
    "mixed": {**PROTEINS, **SHORT_READS, **ISOSEQ},
}

GENEFINDER = {
    # any existing file works as Tiberius model config in a stub run
    "tiberius": {"tiberius": {"run": True, "model_cfg": str(DATA / "tiny.fa")}},
    "vipsania": {"vipsania": {"run": True, "model": "Fungi"}},
}


def expected_outputs(tool: str, mode: str) -> set[str]:
    if mode == "abinitio":
        return {f"{tool}_ab_initio.gff3"}
    return {
        f"{tool}_evidence.gff3",
        f"{tool}_evidence_proteins.fa",
        f"intermediate/{tool}_ab_initio.gff3",
        "intermediate/hc.gff3",
        "hintsfile.gff",
    }


def run_pipeline(tmp_path: Path, params: dict) -> tuple[subprocess.CompletedProcess, set[str]]:
    outdir = tmp_path / "out"
    params = {"genome": str(DATA / "tiny.fa"), "outdir": str(outdir), **params}
    params_file = tmp_path / "params.yaml"
    params_file.write_text(yaml.safe_dump(params))
    proc = subprocess.run(
        [
            NEXTFLOW, "run", str(ROOT / "main.nf"), "-stub-run",
            "-params-file", str(params_file),
            "-c", str(ROOT / "tests" / "stub.config"),
            "-work-dir", str(tmp_path / "work"),
        ],
        cwd=tmp_path, env=dict(os.environ, NXF_ANSI_LOG="false"),
        capture_output=True, text=True,
    )
    published = set()
    if outdir.exists():
        published = {str(p.relative_to(outdir)) for p in outdir.rglob("*") if p.is_file()}
    return proc, published


def assert_ok(proc: subprocess.CompletedProcess) -> None:
    assert proc.returncode == 0, proc.stdout[-4000:] + proc.stderr[-2000:]


@pytest.mark.parametrize("mode", sorted(EVIDENCE))
@pytest.mark.parametrize("tool", sorted(GENEFINDER))
def test_published_files(tool: str, mode: str, tmp_path: Path) -> None:
    proc, published = run_pipeline(tmp_path, {**GENEFINDER[tool], **EVIDENCE[mode]})
    assert_ok(proc)
    missing = expected_outputs(tool, mode) - published
    assert not missing, f"missing {missing}; published {sorted(published)}"
    other = "vipsania" if tool == "tiberius" else "tiberius"
    assert not [f for f in published if Path(f).name.startswith(other)], sorted(published)
    if mode == "abinitio":
        assert f"{tool}_evidence.gff3" not in published


def test_tiberius_file_names_are_unchanged(tmp_path: Path) -> None:
    """The names that Tiberius users rely on (before the gene finder abstraction)."""
    proc, published = run_pipeline(tmp_path, {**GENEFINDER["tiberius"], **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    assert published == {
        "tiberius_evidence.gff3",
        "tiberius_evidence_proteins.fa",
        "intermediate/tiberius_ab_initio.gff3",
        "intermediate/hc.gff3",
        "hintsfile.gff",
    }


def test_legacy_mode_name_tiberius(tmp_path: Path) -> None:
    proc, published = run_pipeline(tmp_path, {**GENEFINDER["tiberius"], "mode": "tiberius"})
    assert_ok(proc)
    assert "Running mode: abinitio" in proc.stdout
    assert "tiberius_ab_initio.gff3" in published


def test_vipsania_is_chunked_without_finetuning(tmp_path: Path) -> None:
    proc, published = run_pipeline(tmp_path, GENEFINDER["vipsania"])
    assert_ok(proc)
    assert "Gene finder : vipsania" in proc.stdout
    logs = {f for f in published if f.startswith("intermediate/vipsania/")}
    # the stub of SPLIT_GENOME makes two chunks
    assert logs == {
        "intermediate/vipsania/vipsania.genome_chunk_1.fa.log",
        "intermediate/vipsania/vipsania.genome_chunk_2.fa.log",
    }


def test_vipsania_finetuning_uses_the_whole_genome(tmp_path: Path) -> None:
    params = {"vipsania": {"run": True, "model": "Fungi", "finetune": True, "finetune_epochs": 1}}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "SPLIT_GENOME" not in proc.stdout
    assert "intermediate/vipsania/vipsania.tiny.fa.log" in published
    assert "vipsania_ab_initio.gff3" in published


def test_vipsania_local_model_dir_skips_the_download(tmp_path: Path) -> None:
    model_dir = tmp_path / "models"
    (model_dir / "abc123").mkdir(parents=True)
    (model_dir / "versions.json").write_text("{}")
    params = {"vipsania": {"run": True, "model": "abc123", "model_dir": str(model_dir)}}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "DOWNLOAD_VIPSANIA_MODEL" not in proc.stdout
    assert "vipsania_ab_initio.gff3" in published


def test_existing_result_is_reused(tmp_path: Path) -> None:
    result = tmp_path / "previous.gtf"
    result.write_text("")
    params = {"vipsania": {"run": True, "result": str(result)}, **PROTEINS}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "RUN_VIPSANIA" not in proc.stdout and "SPLIT_GENOME" not in proc.stdout
    assert "vipsania_evidence.gff3" in published


def test_evidence_without_gene_finder(tmp_path: Path) -> None:
    proc, published = run_pipeline(tmp_path, PROTEINS)
    assert_ok(proc)
    assert "Gene finder : none" in proc.stdout
    assert {"hintsfile.gff", "intermediate/hc.gff3"} <= published
    assert not [f for f in published if "evidence.gff3" in f]


def test_explicit_genefinder_wins(tmp_path: Path) -> None:
    params = {**GENEFINDER["tiberius"], **GENEFINDER["vipsania"], "genefinder": "vipsania"}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "vipsania_ab_initio.gff3" in published


def test_two_gene_finders_are_rejected(tmp_path: Path) -> None:
    proc, _ = run_pipeline(tmp_path, {**GENEFINDER["tiberius"], **GENEFINDER["vipsania"]})
    assert proc.returncode != 0
    assert "More than one gene finder" in proc.stdout + proc.stderr


def test_vipsania_needs_a_model(tmp_path: Path) -> None:
    proc, _ = run_pipeline(tmp_path, {"vipsania": {"run": True}})
    assert proc.returncode != 0
    assert "params.vipsania.model is required" in proc.stdout + proc.stderr
