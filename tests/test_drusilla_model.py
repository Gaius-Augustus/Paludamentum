"""DRUSILLA_ANNOTATE renders `--model <name>`: params.drusilla.model, else drusilla_model.

The process runs for real (not as a stub) in a one-process workflow with a stand-in `drusilla`
that creates the output, like tests/test_no_writes_to_inputs.py; the rendered .command.sh in
the work directory is then inspected. The clade does not choose the model. With shards, the
genome index exists before the shards start (they wrote it all at once before, and one read a
half-written index). Needs a ``nextflow`` executable (or NEXTFLOW_BIN).
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
NEXTFLOW = os.environ.get("NEXTFLOW_BIN") or shutil.which("nextflow")

pytestmark = [
    pytest.mark.nextflow,
    pytest.mark.skipif(NEXTFLOW is None, reason="nextflow executable not found"),
]

FAKE_DRUSILLA = """#!/usr/bin/env bash
# stand-in for `drusilla`: `models download` does nothing, `annotate` writes an orfs.gtf that
# says whether the genome index existed when it started
set -euo pipefail
if [[ $1 == models ]]; then exit 0; fi
[[ $1 == annotate ]] || { echo "unexpected subcommand $1" >&2; exit 64; }
while [[ $# -gt 0 ]]; do
    case $1 in
        --out-dir) out=$2; shift 2 ;;
        --genome) genome=$2; shift 2 ;;
        *) shift ;;
    esac
done
mkdir -p "$out"
if [[ -e $genome.fai ]]; then echo "index yes" > "$out/orfs.gtf"; else echo "index no" > "$out/orfs.gtf"; fi
"""

# stand-in for `gffread`: writes <genome>.fai when none exists, as gffread does
FAKE_GFFREAD = """#!/usr/bin/env bash
set -euo pipefail
while [[ $# -gt 0 ]]; do
    case $1 in
        -g) genome=$2; shift 2 ;;
        *) shift ;;
    esac
done
[[ -e $genome.fai ]] || printf 'c\\t4\\t3\\t4\\t5\\n' > "$genome.fai"
"""

WORKFLOW = f"""nextflow.enable.dsl=2
include {{ DRUSILLA_ANNOTATE }} from '{ROOT.as_posix()}/modules/drusilla.nf'
workflow {{
    DRUSILLA_ANNOTATE(file(params.gtf), file(params.genome), false, file(params.placeholder))
}}
"""

# params.drusilla without a model key; conf/base.config is not loaded
CONFIG = """process {
  executor = 'local'
  cpus     = 1
  memory   = '1 GB'
}
singularity.enabled = false
docker.enabled      = false
params.outdir = '%s'
params.drusilla = [shards: %d, flank_bp: 0, unstranded_both_strands: false%s]
"""

# two genes, one per shard with shards: 2
TWO_GENES = (
    'c\tx\texon\t1\t2\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    'c\tx\texon\t3\t4\t.\t+\t.\tgene_id "g2"; transcript_id "t2";\n'
)


def task_dir(tmp_path: Path, model_setting: str, shards: int = 1, gtf: str = "") -> Path:
    (tmp_path / "in.gtf").write_text(gtf)
    (tmp_path / "genome.fa").write_text(">c\nACGT\n")
    (tmp_path / "placeholder.txt").write_text("")
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    for name, text in (("drusilla", FAKE_DRUSILLA), ("gffread", FAKE_GFFREAD)):
        (fake_bin / name).write_text(text)
        (fake_bin / name).chmod(0o755)
    (tmp_path / "main.nf").write_text(WORKFLOW)
    (tmp_path / "test.config").write_text(CONFIG % (tmp_path / "results", shards, model_setting))
    env = {**os.environ, "PATH": f"{fake_bin}{os.pathsep}{os.environ['PATH']}", "NXF_ANSI_LOG": "false"}
    proc = subprocess.run(
        [NEXTFLOW, "run", str(tmp_path / "main.nf"), "-c", str(tmp_path / "test.config"),
         "--gtf", str(tmp_path / "in.gtf"), "--genome", str(tmp_path / "genome.fa"),
         "--placeholder", str(tmp_path / "placeholder.txt"), "-work-dir", str(tmp_path / "work")],
        cwd=tmp_path, env=env, capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, proc.stdout[-3000:] + proc.stderr[-3000:]
    scripts = list((tmp_path / "work").glob("*/*/.command.sh"))
    assert len(scripts) == 1, scripts
    return scripts[0].parent


def command_sh(tmp_path: Path, model_setting: str) -> str:
    return (task_dir(tmp_path, model_setting) / ".command.sh").read_text()


def test_default_model_is_drusilla_model(tmp_path: Path) -> None:
    script = command_sh(tmp_path, "")
    assert "--model drusilla_model " in script or "--model drusilla_model\n" in script, script
    assert "drusilla models download drusilla_model" in script


def test_params_drusilla_model_is_used(tmp_path: Path) -> None:
    script = command_sh(tmp_path, ", model: 'vertebrates'")
    assert "--model vertebrates" in script and "--model drusilla_model" not in script, script
    assert "drusilla models download vertebrates" in script


def test_shards_start_with_the_genome_index(tmp_path: Path) -> None:
    work = task_dir(tmp_path, "", shards=2, gtf=TWO_GENES)
    assert sorted(p.name for p in work.glob("shard_*.gtf")) == ["shard_0.gtf", "shard_1.gtf"]
    assert (work / "drusilla" / "orfs.gtf").read_text().splitlines() == ["index yes", "index yes"]
