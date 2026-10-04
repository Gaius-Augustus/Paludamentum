"""Command line of Paludamentum: ``paludamentum`` or ``python -m paludamentum``.

Runs the evidence integration pipeline with one of the gene finders that are
submodules of this repository. Inputs come from a params YAML file, from the
command line, or both; command line values override the file. The merged
parameters are written to ``<outdir>/params.yaml``.
"""
from __future__ import annotations

import argparse
from typing import Sequence

from . import __version__
from .launcher import GENEFINDER_CLI, merge_run_params, run_nextflow_pipeline

MODES = ("abinitio", "proteins", "rnaseq", "isoseq", "mixed")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="paludamentum",
        description="Run the Paludamentum evidence integration pipeline with Nextflow.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  Tiberius, ab initio, parallelized over GPUs:\n"
            "    paludamentum --nf_config slurm_generic --genome genome.fa --model_cfg diatoms\n"
            "  Tiberius with evidence from a params file:\n"
            "    paludamentum --params_yaml params.yaml --nf_config slurm_generic\n"
            "  Vipsania with evidence:\n"
            "    paludamentum --genefinder vipsania --nf_config local --genome genome.fa --model Fungi "
            "--proteins proteins.faa\n"
            'Arguments after "--" are forwarded verbatim to Nextflow.'
        ),
    )
    parser.add_argument("--version", action="version", version=f"paludamentum {__version__}")

    general = parser.add_argument_group("General")
    general.add_argument("-p", "--params_yaml",
                         help="Params YAML file, see conf/parameters.yaml. Command line values override it.")
    general.add_argument("-c", "--nf_config", default="base.config",
                         help="Nextflow config: a path, or the name of a config in conf/ "
                              "(e.g. local, slurm_generic). Default: base.config.")
    general.add_argument("--genefinder", choices=sorted(GENEFINDER_CLI),
                         help="Gene finder to run. Default: 'genefinder' or the block with run=true in the "
                              "params file, else vipsania if --model is given, else tiberius.")

    nextflow = parser.add_argument_group("Nextflow")
    nextflow.add_argument("--profile", help="Nextflow profile(s) to activate (comma-separated).")
    nextflow.add_argument("--nextflow_bin", default="nextflow", help="Nextflow executable.")
    nextflow.add_argument("--resume", action="store_true", help="Pass -resume to Nextflow.")
    nextflow.add_argument("--work_dir", help="Nextflow work directory (-work-dir).")
    nextflow.add_argument("--check_tools", action="store_true",
                          help="Check that all tool binaries are available (for runs without containers).")
    nextflow.add_argument("--skip_singularity_check", action="store_true",
                          help="Do not require a singularity executable.")
    nextflow.add_argument("--dry_run", action="store_true",
                          help="Only write the params file and validate inputs and executables; "
                               "do not start Nextflow.")

    inputs = parser.add_argument_group("Inputs and settings (override the params file)")
    inputs.add_argument("--genome", help="Genome FASTA (required unless set in the params file).")
    inputs.add_argument("--outdir", help="Output directory. Default: <genefinder>_results.")
    inputs.add_argument("--threads", type=int, help="Thread count for pipeline processes.")
    inputs.add_argument("--proteins", nargs="*", default=[], help="Protein FASTA file(s).")
    inputs.add_argument("--odb12Partitions", nargs="*", default=[],
                        help="OrthoDB v12 partition name(s) to download and append to the proteins.")
    inputs.add_argument("--rnaseq_single", nargs="*", default=[], help="RNA-Seq single-end FASTQ file(s).")
    inputs.add_argument("--rnaseq_paired", nargs="*", default=[],
                        help='Paired-end RNA-Seq: one quoted glob for all libraries ("RNA/*_{1,2}.fastq.gz"), '
                             "or the two files of one library. Several explicit pairs go into the params file.")
    inputs.add_argument("--rnaseq_bam", nargs="*", default=[], help="Aligned short read BAM file(s).")
    inputs.add_argument("--rnaseq_varus", nargs="*", default=[],
                        help="pyVARUS output directory(ies) of short-read runs, instead of their BAM.")
    inputs.add_argument("--rnaseq_sra_single", nargs="*", default=[], help="RNA-Seq single-end SRA accession(s).")
    inputs.add_argument("--rnaseq_sra_paired", nargs="*", default=[], help="RNA-Seq paired-end SRA accession(s).")
    inputs.add_argument("--isoseq", nargs="*", default=[], help="Iso-Seq FASTQ file(s).")
    inputs.add_argument("--isoseq_sra", nargs="*", default=[], help="Iso-Seq SRA accession(s).")
    inputs.add_argument("--isoseq_varus", nargs="*", default=[],
                        help="pyVARUS output directory(ies) of --longreads runs, instead of their BAM.")
    inputs.add_argument("--mixed_varus",
                        help="Output directory of 'varus assemble --short --long': the assembly of the "
                             "Drusilla flow in mixed mode.")
    inputs.add_argument("--mode", choices=MODES, help="Force the pipeline mode instead of inferring it from the inputs.")
    inputs.add_argument("--scoring_matrix", help="Scoring matrix CSV for miniprot-boundary-scorer.")

    finder = parser.add_argument_group("Gene finder (both)")
    finder.add_argument("--result", help="Existing prediction (GTF/GFF3) to use instead of running the gene finder.")
    finder.add_argument("--min_split_size", type=int, help="Minimal size in bp of a genome chunk.")
    finder.add_argument("--max_files", type=int, help="Maximal number of genome chunks.")
    finder.add_argument("--max_parallel", type=int, help="Cap of concurrently running gene finder tasks.")
    finder.add_argument("--batch_size", type=int, help="Batch size of the gene finder.")

    tiberius = parser.add_argument_group("Tiberius")
    tiberius.add_argument("--model_cfg",
                          help="Model configuration: a path, or a name in model_cfg/ of the Tiberius "
                               "submodule (e.g. diatoms).")
    tiberius.add_argument("--seq_len", type=int, help="Forwarded to tiberius.py --seq_len.")
    finder.add_argument("--model_dir",
                        help="Directory with downloaded models, instead of the download at the start of the run: the "
                             "extracted weights of the Tiberius model configuration, or 'vipsania download -d DIR'.")

    vipsania = parser.add_argument_group("Vipsania")
    vipsania.add_argument("--model", help="Clade name (e.g. Fungi) or model id.")
    vipsania.add_argument("--finetune", action="store_true",
                          help="Finetune on the genome before annotating (one task, no chunks).")
    vipsania.add_argument("--finetune_epochs", type=int, help="Forwarded to vipsania annotate --finetune_epochs.")
    vipsania.add_argument("--context", type=int, help="Forwarded to vipsania annotate -T.")

    parser.add_argument("nextflow_args", nargs=argparse.REMAINDER, help=argparse.SUPPRESS)
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    # params.yaml is written only after the validation passed
    genefinder, params = merge_run_params(args)
    note = "" if params[genefinder].get("run") else " (run: false, evidence only)"
    print(f"[INFO] Gene finder: {genefinder}{note}")
    run_nextflow_pipeline(args, genefinder=genefinder, params=params)


if __name__ == "__main__":
    main()
