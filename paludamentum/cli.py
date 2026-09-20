"""Standalone command line: ``python -m paludamentum``.

The gene finders bring their own, richer command lines (``tiberius.py``,
``vipsania annotate``). This one expects a complete params YAML file.
"""
from __future__ import annotations

import argparse
from typing import Sequence

from . import __version__
from .launcher import GENEFINDER_CLI, run_nextflow_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="paludamentum",
        description="Run the Paludamentum evidence integration pipeline with Nextflow.",
        epilog='Arguments after "--" are forwarded verbatim to Nextflow.',
    )
    parser.add_argument("--version", action="version", version=f"paludamentum {__version__}")
    parser.add_argument("-p", "--params_yaml", required=True,
                        help="Params YAML file, see conf/parameters.yaml.")
    parser.add_argument("-c", "--nf_config", default="base.config",
                        help="Nextflow config: a path, or the name of a config in conf/ "
                             "(e.g. local, slurm_generic). Default: base.config.")
    parser.add_argument("--genefinder", choices=sorted(GENEFINDER_CLI), default="tiberius",
                        help="Gene finder that the params file configures (default: tiberius).")
    parser.add_argument("--profile", help="Nextflow profile(s) to activate (comma-separated).")
    parser.add_argument("--nextflow_bin", default="nextflow", help="Nextflow executable.")
    parser.add_argument("--resume", action="store_true", help="Pass -resume to Nextflow.")
    parser.add_argument("--work_dir", help="Nextflow work directory (-work-dir).")
    parser.add_argument("--check_tools", action="store_true",
                        help="Check that all tool binaries are available (for runs without containers).")
    parser.add_argument("--skip_singularity_check", action="store_true",
                        help="Do not require a singularity executable.")
    parser.add_argument("--dry_run", action="store_true",
                        help="Only validate inputs and executables; do not start Nextflow.")
    parser.add_argument("nextflow_args", nargs=argparse.REMAINDER,
                        help=argparse.SUPPRESS)
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    run_nextflow_pipeline(args, genefinder=args.genefinder)


if __name__ == "__main__":
    main()
