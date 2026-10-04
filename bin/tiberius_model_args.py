#!/usr/bin/env python3
"""Point Tiberius at pre-downloaded weights instead of letting it download them.

With --model_cfg, Tiberius looks for the extracted weights only in
<tiberius package>/../model_weights, and falls back to the working directory
only if that directory is missing or read-only. The Tiberius image creates it
world-writable, so under Docker the weights staged into the task directory are
ignored and Tiberius downloads them, which fails on nodes without internet.

This script reads the model configuration, finds the extracted weights
directory in WEIGHTS_DIR (the archive name of weights_url without its
extensions, as Tiberius names it) and prints the arguments that run the same
model with --model instead of --model_cfg, one per line (read them with
`mapfile -t`, so paths with spaces stay one argument):

    --model <abs path> [--seq_len <default_seq_len>] [--no_softmasking]

--seq_len is printed only if none is given (Tiberius takes default_seq_len
from the configuration only with --model_cfg). It fails if the weights
directory is missing or empty, so a wrong model_dir never starts a download.

Usage: tiberius_model_args.py --model_cfg CFG.yaml [--weights_dir DIR] [--seq_len N]
"""
from __future__ import annotations

import argparse
import os
import sys

import yaml


def weights_name(weights_url: str) -> str:
    """Directory the weights archive extracts to (tiberius.main.resolve_weight_download)."""
    return weights_url.split("/")[-1].split(".")[0]


def model_args(model_cfg: str, weights_dir: str, seq_len: int | None = None) -> list[str]:
    with open(model_cfg, encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    for key in ("weights_url", "softmasking", "default_seq_len"):
        if key not in config:
            raise ValueError(f"{model_cfg} has no {key}")
    name = weights_name(config["weights_url"])
    path = os.path.join(weights_dir, name)
    if not os.path.isdir(path) or not os.listdir(path):
        present = sorted(os.listdir(weights_dir)) if os.path.isdir(weights_dir) else []
        raise FileNotFoundError(
            f"weights directory '{name}' for {model_cfg} not found in tiberius.model_dir "
            f"(found: {', '.join(present) or 'nothing'}). Download {config['weights_url']} "
            f"and extract it in tiberius.model_dir."
        )
    args = ["--model", os.path.abspath(path)]
    if not seq_len:
        args += ["--seq_len", str(int(config["default_seq_len"]))]
    if not config["softmasking"]:
        args.append("--no_softmasking")
    return args


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--model_cfg", required=True, help="Tiberius model configuration (YAML)")
    parser.add_argument("--weights_dir", default=".", help="directory with the extracted weights (default: .)")
    parser.add_argument("--seq_len", type=int, help="seq_len passed to Tiberius, if any")
    args = parser.parse_args(argv)
    try:
        print("\n".join(model_args(args.model_cfg, args.weights_dir, args.seq_len)))
    except (OSError, ValueError, yaml.YAMLError) as exc:
        print(f"[tiberius_model_args] {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
