"""HINT_RESCUE_TIBERIUS with its real script: which model the rescue passes to Tiberius.

A fake tiberius.py records its arguments; fake bin/ scripts stand in for
tiberius_batch_size.py, tiberius_model_args.py, filter_and_merge_rescue_gtf.py
and fix_cds_phases.py.
Needs a ``nextflow`` executable (or NEXTFLOW_BIN), like test_stub_run.py.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from test_drusilla_stub_run import run_process
from test_stub_run import NEXTFLOW, assert_ok

pytestmark = [
    pytest.mark.nextflow,
    pytest.mark.skipif(NEXTFLOW is None, reason="nextflow executable not found"),
]

FAKE_TIBERIUS = r'''#!/usr/bin/env python3
import json, sys
a = sys.argv[1:]
json.dump(a, open(a[a.index("--out") + 1], "w"))
'''

FAKES = {
    # records the model configuration it was given, prints no batch size
    "tiberius_batch_size.py": '#!/bin/sh\necho "$2" > batch_model_cfg.txt\n',
    # one argument per line, as bin/tiberius_model_args.py; the path has a space
    "tiberius_model_args.py": '#!/bin/sh\nprintf -- "--model\\n%s/my weights/w\\n--seq_len\\n99990\\n" "$PWD"\n',
    "filter_and_merge_rescue_gtf.py": '#!/bin/sh\ncp "$1" "$4"\n',
    "fix_cds_phases.py": '#!/bin/sh\ncat "$1"\n',
}

RESCUE = """
    loci = channel.of(tuple(file('%(dir)s/loci.fa'), file('%(dir)s/hints.gff'), file('%(dir)s/manifest.tsv')))
    weights = %(weights)s
    HINT_RESCUE_TIBERIUS(loci, '%(name)s', %(file)s, %(use)s, weights)
"""


def rescue_args(tmp_path: Path, name: str = "", model_file: Path | None = None,
                weights: bool = False) -> list[str]:
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    for f in ("loci.fa", "hints.gff", "manifest.tsv"):
        (inputs / f).write_text("x\n")
    tiberius = tmp_path / "tiberius" / "tiberius.py"
    tiberius.parent.mkdir()
    tiberius.write_text(FAKE_TIBERIUS)
    wdir = tmp_path / "weights"
    (wdir / "w").mkdir(parents=True)
    body = RESCUE % {
        "dir": inputs, "name": name, "use": str(weights).lower(),
        "file": f"file('{model_file}')" if model_file else "[]",
        "weights": f"channel.fromPath('{wdir}/*', type: 'any').collect()" if weights else "[]",
    }
    drusilla = {"rescue_tiberius": str(tiberius), "rescue_seq_len": 99990, "rescue_hint_weight": 2.5}
    proc, out = run_process(tmp_path, "HINT_RESCUE_TIBERIUS", body, drusilla, FAKES)
    assert_ok(proc)
    return json.loads((out / "intermediate" / "hint_rescue.gtf").read_text())


def model_part(args: list[str]) -> list[str]:
    """The arguments between --genome GENOME and --hints."""
    return args[2:args.index("--hints")]


def test_a_model_of_the_image_is_passed_by_name(tmp_path: Path) -> None:
    args = rescue_args(tmp_path, name="vertebrates")
    assert model_part(args) == ["--model_cfg", "vertebrates"]


def test_a_model_file_is_staged_and_passed(tmp_path: Path) -> None:
    cfg = tmp_path / "custom models" / "my vertebrates.yaml"
    cfg.parent.mkdir()
    cfg.write_text('target_species: "Vertebrata"\n')
    args = rescue_args(tmp_path, model_file=cfg)
    assert model_part(args) == ["--model_cfg", "my vertebrates.yaml"]


def test_staged_weights_are_passed_with_model(tmp_path: Path) -> None:
    cfg = tmp_path / "vertebrates.yaml"
    cfg.write_text('target_species: "Vertebrata"\n')
    args = rescue_args(tmp_path, model_file=cfg, weights=True)
    model, path, *rest = model_part(args)
    assert model == "--model" and path.endswith("/my weights/w")
    assert rest == ["--seq_len", "99990"]
