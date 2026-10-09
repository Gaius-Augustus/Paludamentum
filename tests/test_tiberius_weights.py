"""DOWNLOAD_TIBERIUS_WEIGHTS and RUN_TIBERIUS with their real scripts.

The download fetches a local archive through a file:// URL. RUN_TIBERIUS runs
with a fake tiberius.py that records its arguments, and fake bin/ scripts for
tiberius_batch_size.py and tiberius_model_args.py.
Needs a ``nextflow`` executable (or NEXTFLOW_BIN), like test_stub_run.py.
"""
from __future__ import annotations

import json
import os
import subprocess
import tarfile
from pathlib import Path

import pytest

from test_stub_run import DATA, NEXTFLOW, ROOT, assert_ok

pytestmark = [
    pytest.mark.nextflow,
    pytest.mark.skipif(NEXTFLOW is None, reason="nextflow executable not found"),
]

PROCESS_TEST = """
include { %s } from '%s/modules/genefinder.nf'
params.outdir = 'out'
workflow {
%s
}
"""


def run_process(tmp_path: Path, process: str, body: str,
                fakes: dict[str, str] | None = None) -> tuple[subprocess.CompletedProcess, Path]:
    """Runs one process of modules/genefinder.nf without -stub-run; returns the outdir."""
    (tmp_path / "bin").mkdir()
    for name, text in (fakes or {}).items():
        (tmp_path / "bin" / name).write_text(text)
        (tmp_path / "bin" / name).chmod(0o755)
    script = tmp_path / "test.nf"
    script.write_text(PROCESS_TEST % (process, ROOT, body))
    config = tmp_path / "test.config"
    config.write_text("process.executor = 'local'\nprocess.container = null\nprocess.cpus = 1\n")
    outdir = tmp_path / "out"
    proc = subprocess.run(
        [NEXTFLOW, "run", str(script), "-c", str(config), "--outdir", str(outdir),
         "-work-dir", str(tmp_path / "work")],
        cwd=tmp_path, env=dict(os.environ, NXF_ANSI_LOG="false"), capture_output=True, text=True,
    )
    return proc, outdir


# ---- DOWNLOAD_TIBERIUS_WEIGHTS ----------------------------------------------------

DOWNLOAD = """
    DOWNLOAD_TIBERIUS_WEIGHTS('%s').flatten()
        .subscribe { d -> d.copyTo("${params.outdir}/${d.name}") }
"""


def archive(tmp_path: Path, name: str, top: str) -> str:
    """A <name>.tar.gz with the directory top; returns its file:// URL."""
    src = tmp_path / "src" / top
    src.mkdir(parents=True)
    (src / "weights.h5").write_text("w\n")
    path = tmp_path / f"{name}.tar.gz"
    with tarfile.open(path, "w:gz") as tar:
        tar.add(src, arcname=top)
    return path.as_uri()


def test_the_weights_are_downloaded_and_extracted(tmp_path: Path) -> None:
    url = archive(tmp_path, "tiny_weights", "tiny_weights")
    proc, out = run_process(tmp_path, "DOWNLOAD_TIBERIUS_WEIGHTS", DOWNLOAD % url)
    assert_ok(proc)
    assert (out / "tiny_weights" / "weights.h5").read_text() == "w\n"
    # the archive itself is not emitted
    assert sorted(p.name for p in out.iterdir()) == ["tiny_weights"]


def test_an_archive_without_the_weights_directory_fails(tmp_path: Path) -> None:
    url = archive(tmp_path, "tiny_weights", "other")
    proc, _ = run_process(tmp_path, "DOWNLOAD_TIBERIUS_WEIGHTS", DOWNLOAD % url)
    assert proc.returncode != 0
    assert "has no directory tiny_weights, found: other" in proc.stdout + proc.stderr


def test_a_failed_download_fails(tmp_path: Path) -> None:
    proc, _ = run_process(tmp_path, "DOWNLOAD_TIBERIUS_WEIGHTS",
                          DOWNLOAD % (tmp_path / "missing.tar.gz").as_uri())
    assert proc.returncode != 0


# ---- RUN_TIBERIUS -----------------------------------------------------------------

FAKES = {
    "tiberius.py": r'''#!/usr/bin/env python3
import json, sys
a = sys.argv[1:]
json.dump(a, open(a[a.index("--out") + 1], "w"))
''',
    "tiberius_batch_size.py": "#!/bin/sh\n",
    # one argument per line, as bin/tiberius_model_args.py
    "tiberius_model_args.py": '#!/bin/sh\nprintf -- "--model\\n%s/tiny_weights\\n--seq_len\\n99990\\n" "$PWD"\n',
}

RUN = """
    RUN_TIBERIUS(file('%(genome)s'), file('%(cfg)s'), %(use)s, %(weights)s)
        .gtf.subscribe { f -> f.copyTo("${params.outdir}/${f.name}") }
"""


@pytest.mark.parametrize("use_weights", [True, False], ids=["staged-weights", "no-weights"])
def test_run_tiberius_model_arguments(use_weights: bool, tmp_path: Path) -> None:
    cfg = tmp_path / "tiny.yaml"
    cfg.write_text('weights_url: "https://example.org/models/tiny_weights.tar.gz"\n')
    (tmp_path / "weights" / "tiny_weights").mkdir(parents=True)
    body = RUN % {
        "genome": DATA / "tiny.fa", "cfg": cfg, "use": str(use_weights).lower(),
        "weights": f"channel.fromPath('{tmp_path}/weights/*', type: 'any').collect()" if use_weights else "[]",
    }
    proc, out = run_process(tmp_path, "RUN_TIBERIUS", body, FAKES)
    assert_ok(proc)
    args = json.loads((out / "tiberius.tiny.fa.gtf").read_text())
    model = args[2:args.index("--out")]
    if use_weights:
        assert model[0] == "--model" and model[1].endswith("/tiny_weights")
        assert model[2:] == ["--seq_len", "99990"]
    else:
        assert model == ["--model_cfg", "tiny.yaml"]
