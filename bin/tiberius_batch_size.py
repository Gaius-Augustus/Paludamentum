#!/usr/bin/env python3
"""Cap the Tiberius batch size where its automatic choice would crash.

Tiberius sizes its batch from the GPU memory. On GPUs with more than 80 GB
(for example the 96 GB RTX PRO 6000 Blackwell) it picks a batch so large that
batch x seq_len x 256 channels exceeds the int32 range of TensorFlow's GPU
kernel launch configuration, and every task aborts with
"Check failed: work_element_count >= 0".

This script repeats the automatic choice of Tiberius (compute_auto_batch_size
in tiberius/main.py) and prints "--batch_size N" only if that choice exceeds
the int32 limit. Otherwise it prints nothing and Tiberius decides as usual.
It never fails the task: on any problem it prints nothing.

Usage: tiberius_batch_size.py --model_cfg CFG.yaml [--seq_len N] [--gpu_mb MB]
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys

INT32_MAX = 2**31 - 1
CHANNELS = 256                 # widest per-position layer of the Tiberius models
TIBERIUS_DEFAULT_SEQ_LEN = 500_040
PREFERRED_BATCH_SIZES = [1, 2, 4, 8, 16, 24, 32, 64, 96, 128, 160, 192]


def tiberius_auto_batch_size(seq_len: int, gpu_memory_gb: float, safety_factor: float = 0.9) -> int:
    """Same estimate as tiberius.main.compute_auto_batch_size."""
    if gpu_memory_gb >= 70:
        ref_mem, ref_bs = 80.0, 18
    elif gpu_memory_gb >= 16:
        ref_mem, ref_bs = 25.0, 8
    else:
        ref_mem, ref_bs = 8.0, 2
    estimated = (gpu_memory_gb / ref_mem) * ref_bs * (500_004 / seq_len) * safety_factor
    batch_size = max(1, estimated)
    candidates = [b for b in PREFERRED_BATCH_SIZES if abs(b - batch_size) < 3]
    if candidates:
        batch_size = min(candidates, key=lambda b: abs(b - batch_size))
    return int(batch_size + 0.5)


def int32_safe_batch_size(seq_len: int) -> int:
    return max(1, INT32_MAX // (seq_len * CHANNELS))


def seq_len_from_model_cfg(path: str) -> int | None:
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            match = re.match(r"\s*default_seq_len\s*:\s*([0-9_]+)", line)
            if match:
                return int(match.group(1).replace("_", ""))
    return None


def gpu_memory_mb() -> float | None:
    """Memory of the first visible GPU, as Tiberius reads it."""
    try:
        out = subprocess.run(
            ["nvidia-smi", "--query-gpu=memory.total", "--format=csv,noheader,nounits"],
            check=True, capture_output=True, text=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        return None
    lines = [line.strip() for line in out.splitlines() if line.strip()]
    return float(lines[0]) if lines else None


def batch_size_cap(seq_len: int, gpu_mb: float) -> int | None:
    """The batch size to force, or None if the automatic choice is safe."""
    auto = tiberius_auto_batch_size(seq_len, gpu_mb / 1024.0)
    limit = int32_safe_batch_size(seq_len)
    return limit if auto > limit else None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--model_cfg", help="Tiberius model configuration (reads default_seq_len)")
    parser.add_argument("--seq_len", type=int, help="seq_len passed to Tiberius, if any")
    parser.add_argument("--gpu_mb", type=float, help="GPU memory in MiB (default: nvidia-smi)")
    args = parser.parse_args(argv)
    try:
        seq_len = args.seq_len
        if not seq_len and args.model_cfg:
            seq_len = seq_len_from_model_cfg(args.model_cfg)
        seq_len = seq_len or TIBERIUS_DEFAULT_SEQ_LEN
        gpu_mb = args.gpu_mb if args.gpu_mb is not None else gpu_memory_mb()
        if gpu_mb is None:
            return 0
        cap = batch_size_cap(seq_len, gpu_mb)
    except Exception as exc:  # never fail the task because of this helper
        print(f"[tiberius_batch_size] {exc}; leaving the batch size to Tiberius", file=sys.stderr)
        return 0
    if cap is not None:
        print(f"[tiberius_batch_size] GPU {gpu_mb:.0f} MiB, seq_len {seq_len}: automatic batch size "
              f"would overflow int32; using --batch_size {cap}", file=sys.stderr)
        print(f"--batch_size {cap}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
