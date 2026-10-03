"""bin/tiberius_batch_size.py: cap Tiberius' automatic batch size where it overflows int32."""
from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "bin" / "tiberius_batch_size.py"
spec = importlib.util.spec_from_file_location("tiberius_batch_size", SCRIPT)
tbs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tbs)

VERTEBRATES_SEQ_LEN = 400_050   # default_seq_len of Tiberius model_cfg/vertebrates.yaml


def run(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)


def test_rtx_pro_6000_96gb_is_capped():
    """A 96 GB GPU (RTX PRO 6000): the automatic 24 x 400050 x 256 overflows int32."""
    assert tbs.tiberius_auto_batch_size(VERTEBRATES_SEQ_LEN, 97887 / 1024) == 24
    assert tbs.batch_size_cap(VERTEBRATES_SEQ_LEN, 97887) == 20
    assert 20 * VERTEBRATES_SEQ_LEN * 256 <= tbs.INT32_MAX < 24 * VERTEBRATES_SEQ_LEN * 256


@pytest.mark.parametrize("gpu_mb", [81920, 40960, 24576, 8192])   # A100 80/40 GB, 24 GB, 8 GB
def test_smaller_gpus_keep_the_automatic_choice(gpu_mb):
    assert tbs.batch_size_cap(VERTEBRATES_SEQ_LEN, gpu_mb) is None


def test_seq_len_is_read_from_the_model_cfg(tmp_path: Path):
    cfg = tmp_path / "m.yaml"
    cfg.write_text('target_species: "x"\ndefault_seq_len: 400_050\n')
    assert tbs.seq_len_from_model_cfg(str(cfg)) == 400_050
    proc = run("--model_cfg", str(cfg), "--gpu_mb", "97887")
    assert proc.returncode == 0 and proc.stdout.strip() == "--batch_size 20"


def test_real_vertebrates_cfg():
    cfg = ROOT / "tiberius" / "model_cfg" / "vertebrates.yaml"
    if not cfg.is_file():
        pytest.skip("Tiberius submodule not checked out")
    assert tbs.seq_len_from_model_cfg(str(cfg)) == VERTEBRATES_SEQ_LEN


def test_explicit_seq_len_wins(tmp_path: Path):
    """batch x seq_len is about constant in Tiberius' estimate: 96 GB overflows at any seq_len."""
    cfg = tmp_path / "m.yaml"
    cfg.write_text("default_seq_len: 400_050\n")
    assert tbs.tiberius_auto_batch_size(99_990, 97887 / 1024) == 96
    assert run("--model_cfg", str(cfg), "--seq_len", "99990", "--gpu_mb", "97887").stdout.strip() == "--batch_size 83"
    assert run("--model_cfg", str(cfg), "--seq_len", "99990", "--gpu_mb", "81920").stdout == ""


def test_never_fails(tmp_path: Path):
    # no nvidia-smi, missing model_cfg: prints nothing, exit 0
    assert run("--gpu_mb", "97887", "--model_cfg", str(tmp_path / "missing.yaml")).returncode == 0
    proc = subprocess.run([sys.executable, str(SCRIPT), "--model_cfg", str(tmp_path / "missing.yaml")],
                          capture_output=True, text=True, env={"PATH": "/nonexistent"})
    assert proc.returncode == 0 and proc.stdout == ""
