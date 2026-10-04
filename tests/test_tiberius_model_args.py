"""bin/tiberius_model_args.py: run Tiberius with pre-downloaded weights via --model."""
from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "bin" / "tiberius_model_args.py"
spec = importlib.util.spec_from_file_location("tiberius_model_args", SCRIPT)
tma = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tma)

MODEL_CFG = ROOT / "tiberius" / "model_cfg"


def run(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)


def weights(tmp_path: Path, name: str) -> Path:
    path = tmp_path / name
    path.mkdir()
    (path / "saved_model.pb").write_text("")
    return path


@pytest.mark.parametrize("url, name", [
    ("https://x/models/tiberius_weights_v2.tar.gz", "tiberius_weights_v2"),
    ("https://x/models/fungi.tar.gz", "fungi"),
    ("https://x/models//tiberius_denovo_weights_v2.tar.gz", "tiberius_denovo_weights_v2"),
])
def test_weights_name_as_tiberius(url, name):
    assert tma.weights_name(url) == name


@pytest.mark.skipif(not (MODEL_CFG / "mammalia_softmasking_v2.yaml").exists(), reason="tiberius submodule not checked out")
def test_softmasking_model(tmp_path):
    path = weights(tmp_path, "tiberius_weights_v2")
    args = tma.model_args(str(MODEL_CFG / "mammalia_softmasking_v2.yaml"), str(tmp_path))
    assert args == ["--model", str(path), "--seq_len", "400050"]


@pytest.mark.skipif(not (MODEL_CFG / "fungi.yaml").exists(), reason="tiberius submodule not checked out")
def test_unmasked_model_keeps_given_seq_len(tmp_path):
    path = weights(tmp_path, "fungi")
    proc = run("--model_cfg", str(MODEL_CFG / "fungi.yaml"), "--weights_dir", str(tmp_path), "--seq_len", "99999")
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.split() == ["--model", str(path), "--no_softmasking"]


def test_missing_weights_fail(tmp_path):
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text('weights_url: "https://x/fungi.tar.gz"\nsoftmasking: False\ndefault_seq_len: 100_080\n')
    (tmp_path / "w").mkdir()
    weights(tmp_path / "w", "insecta")
    proc = run("--model_cfg", str(cfg), "--weights_dir", str(tmp_path / "w"))
    assert proc.returncode == 1
    assert proc.stdout == ""
    assert "'fungi'" in proc.stderr and "insecta" in proc.stderr


def test_empty_weights_dir_fails(tmp_path):
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text('weights_url: "https://x/fungi.tar.gz"\nsoftmasking: False\ndefault_seq_len: 100_080\n')
    (tmp_path / "fungi").mkdir()
    assert run("--model_cfg", str(cfg), "--weights_dir", str(tmp_path)).returncode == 1


def test_one_argument_per_line_keeps_spaces(tmp_path):
    cfg = tmp_path / "my cfg.yaml"
    cfg.write_text('weights_url: "https://x/fungi.tar.gz"\nsoftmasking: True\ndefault_seq_len: 100_080\n')
    path = weights(tmp_path, "fungi")
    proc = run("--model_cfg", str(cfg), "--weights_dir", str(tmp_path))
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.splitlines() == ["--model", str(path), "--seq_len", "100080"]
