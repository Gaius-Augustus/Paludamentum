"""No task may write to its staged inputs.

With process.scratch = true and the default symlink staging, an input in the
task directory is a symlink to the user's original file, and a write to it
changes the original. STRAND_FIX_BAM once wrote its output to the staged input
path (`${bam.name}` of an input staged as unfixed/* renders as unfixed/<name>,
while `${bam.baseName}` is <name> without the directory) and so truncated a
user's BAM; SAMTOOLS_SORT had the same pattern.

The static test rejects that pattern in every module. The Nextflow tests run
STRAND_FIX_BAM and SAMTOOLS_SORT for real (not as stubs) with stand-in tools
that refuse an output resolving to the input, and check that the original is
unchanged and the output is a new file. They need a ``nextflow`` executable
(or NEXTFLOW_BIN), like tests/test_stub_run.py.
"""
from __future__ import annotations

import hashlib
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
NEXTFLOW = os.environ.get("NEXTFLOW_BIN") or shutil.which("nextflow")
needs_nextflow = pytest.mark.skipif(NEXTFLOW is None, reason="nextflow executable not found")

STAGED_IN_DIR = re.compile(r"""path\s*\(?\s*(\w+)\s*,\s*stageAs:\s*["'][^"'*]*/\*["']""")
PROCESS = re.compile(r"^process\s+(\w+)\s*\{", re.MULTILINE)


def process_blocks():
    for path in sorted((ROOT / "modules").glob("*.nf")):
        text = path.read_text()
        starts = [(m.start(), m.group(1)) for m in PROCESS.finditer(text)]
        for i, (start, name) in enumerate(starts):
            end = starts[i + 1][0] if i + 1 < len(starts) else len(text)
            yield path.name, name, text[start:end]


def test_no_output_is_named_after_a_staged_input():
    """`${x.name}` of an input staged in a directory is the input's own path: an
    output, a redirect or a stub `touch` of that name writes to the original."""
    offenders = []
    for module, process, block in process_blocks():
        lines = [line for line in block.splitlines() if not line.strip().startswith("//")]
        body = "\n".join(lines)
        for var in STAGED_IN_DIR.findall(body):
            # `tag` only labels the task; everything else is a file name
            uses = [line.strip() for line in lines
                    if "${" + var + ".name}" in line and not line.strip().startswith("tag ")]
            if uses:
                offenders.append(f"{module}:{process}: {uses}")
    assert not offenders, offenders


FAKE_DRUSILLA = """#!/usr/bin/env bash
# stand-in for `drusilla strand-fix-bam`: copies the input, refuses to write onto it
set -euo pipefail
[[ $1 == strand-fix-bam ]] || { echo "unexpected subcommand $1" >&2; exit 64; }
shift
while [[ $# -gt 0 ]]; do
    case $1 in
        --in-bam) in=$2; shift 2 ;;
        --out-bam) out=$2; shift 2 ;;
        --report-tsv) report=$2; shift 2 ;;
        --threads|--run-regex) shift 2 ;;
        *) shift ;;
    esac
done
if [[ -e $out && $(readlink -f "$out") == $(readlink -f "$in") ]]; then
    echo "out-bam $out is the input $in" >&2
    exit 65
fi
cp "$in" "$out"
printf 'run\\tprotocol\\nSRR1\\tunstranded\\n' > "$report"
"""

FAKE_SAMTOOLS = """#!/usr/bin/env bash
# stand-in for `samtools sort -@ N -o OUT IN`: copies the input, refuses to write onto it
set -euo pipefail
[[ $1 == sort ]] || { echo "unexpected subcommand $1" >&2; exit 64; }
shift
while [[ $# -gt 0 ]]; do
    case $1 in
        -o) out=$2; shift 2 ;;
        -@) shift 2 ;;
        *) in=$1; shift ;;
    esac
done
if [[ -e $out && $(readlink -f "$out") == $(readlink -f "$in") ]]; then
    echo "sort output $out is the input $in" >&2
    exit 65
fi
cp "$in" "$out"
"""

WORKFLOWS = {
    "STRAND_FIX_BAM": "STRAND_FIX_BAM(channel.of([file(params.bam), '']))\n"
                      "    STRAND_FIX_BAM.out.bam.view { f -> \"OUT_BAM ${f}\" }",
    "SAMTOOLS_SORT": "SAMTOOLS_SORT(channel.of(file(params.bam)))\n"
                     "    SAMTOOLS_SORT.out.bam.view { f -> \"OUT_BAM ${f}\" }",
}

CONFIG = """process {
  executor = 'local'
  scratch  = true
  cpus     = 1
  memory   = '1 GB'
}
singularity.enabled = false
docker.enabled      = false
params.outdir = '%s'
params.tools = [samtools: 'samtools']
"""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_process(tmp_path: Path, process: str, fake_name: str, fake_script: str) -> tuple[Path, str]:
    original = tmp_path / "user_data" / "VARUS.bam"
    original.parent.mkdir()
    original.write_bytes(b"BAM\x01" + os.urandom(4096))
    before = sha256(original)
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    (fake_bin / fake_name).write_text(fake_script)
    (fake_bin / fake_name).chmod(0o755)
    (tmp_path / "main.nf").write_text(
        "nextflow.enable.dsl=2\n"
        f"include {{ {process} }} from '{ROOT}/modules/rnaseq.nf'\n"
        "workflow {\n    " + WORKFLOWS[process] + "\n}\n")
    (tmp_path / "test.config").write_text(CONFIG % (tmp_path / "results"))
    env = {**os.environ, "PATH": f"{fake_bin}{os.pathsep}{os.environ['PATH']}", "NXF_ANSI_LOG": "false"}
    proc = subprocess.run(
        [NEXTFLOW, "run", str(tmp_path / "main.nf"), "-c", str(tmp_path / "test.config"),
         "--bam", str(original), "-work-dir", str(tmp_path / "work")],
        cwd=tmp_path, env=env, capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, proc.stdout[-3000:] + proc.stderr[-3000:]
    assert sha256(original) == before, f"{process} changed the user's original BAM"
    out = [line.split(" ", 1)[1] for line in proc.stdout.splitlines() if line.startswith("OUT_BAM ")]
    assert len(out) == 1, proc.stdout
    out_bam = Path(out[0])
    assert not out_bam.is_symlink() and out_bam.resolve() != original.resolve(), out_bam
    assert sha256(out_bam) == before
    return out_bam, proc.stdout


@pytest.mark.nextflow
@needs_nextflow
def test_strand_fix_bam_leaves_the_staged_original_alone(tmp_path: Path) -> None:
    out_bam, _ = run_process(tmp_path, "STRAND_FIX_BAM", "drusilla", FAKE_DRUSILLA)
    assert out_bam.name == "VARUS.bam" and out_bam.parent.name != "unfixed", out_bam


@pytest.mark.nextflow
@needs_nextflow
def test_samtools_sort_leaves_the_staged_original_alone(tmp_path: Path) -> None:
    out_bam, _ = run_process(tmp_path, "SAMTOOLS_SORT", "samtools", FAKE_SAMTOOLS)
    assert out_bam.name == "VARUS.bam" and out_bam.parent.name != "unsorted", out_bam
