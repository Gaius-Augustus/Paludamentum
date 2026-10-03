"""Launcher for the Paludamentum Nextflow pipeline.

Builds the params file of a run from defaults, a params YAML file and command
line values, validates inputs and executables, then runs ``nextflow run
main.nf``. The command line is ``paludamentum`` (``python -m paludamentum``).

The gene finders (Tiberius, Vipsania) and Drusilla are git submodules of this
repository. The launcher uses the Tiberius checkout to resolve model
configuration names and, for runs without containers, to find ``tiberius.py``.
"""
from __future__ import annotations

import glob
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import yaml


GLOB_CHARS = set("*?[]{}")

DEFAULT_TOOL_BINARIES: Dict[str, str] = {
    "hisat2": "hisat2",
    "hisat2_build": "hisat2-build",
    "minimap2": "minimap2",
    "stringtie": "stringtie",
    "samtools": "samtools",
    "transdecoder_longorfs": "TransDecoder.LongOrfs",
    "transdecoder_predict": "TransDecoder.Predict",
    "transdecoder_util_gtf2fa": "gtf_genome_to_cdna_fasta.pl",
    "transdecoder_util_orf2genome": "cdna_alignment_orf_to_genome_orf.pl",
    "transdecoder_gtf2gff": "gtf_to_alignment_gff3.pl",
    "diamond": "diamond",
    "bedtools": "bedtools",
    "miniprot": "miniprot",
    "miniprot_boundary_scorer": "miniprot_boundary_scorer",
    "miniprothint": "miniprothint.py",
    "bam2hints": "bam2hints",
}

TOOL_DESCRIPTIONS: Dict[str, str] = {
    "hisat2": "HISAT2 aligner (params.tools.hisat2)",
    "hisat2_build": "HISAT2 indexer (params.tools.hisat2_build)",
    "minimap2": "Minimap2 aligner (params.tools.minimap2)",
    "stringtie": "StringTie assembler (params.tools.stringtie)",
    "samtools": "Samtools (params.tools.samtools)",
    "transdecoder_longorfs": "TransDecoder.LongOrfs (params.tools.transdecoder_longorfs)",
    "transdecoder_predict": "TransDecoder.Predict (params.tools.transdecoder_predict)",
    "transdecoder_util_gtf2fa": "gtf_genome_to_cdna_fasta.pl (params.tools.transdecoder_util_gtf2fa)",
    "transdecoder_util_orf2genome": "cdna_alignment_orf_to_genome_orf.pl (params.tools.transdecoder_util_orf2genome)",
    "transdecoder_gtf2gff": "gtf_to_alignment_gff3.pl (params.tools.transdecoder_gtf2gff)",
    "diamond": "DIAMOND aligner (params.tools.diamond)",
    "bedtools": "BEDTools (params.tools.bedtools)",
    "miniprot": "MiniProt aligner (params.tools.miniprot)",
    "miniprot_boundary_scorer": "MiniProt boundary scorer (params.tools.miniprot_boundary_scorer)",
    "miniprothint": "MiniProtHint (params.tools.miniprothint)",
    "bam2hints": "BAM2HINTS (params.tools.bam2hints)",
}

GENERAL_COMMANDS = {
    "nextflow": "Nextflow executable used to launch the pipeline",
    "java": "Java runtime (version 17 or newer, required by Nextflow)",
    "singularity": "Singularity/Apptainer runtime",
    "python3": "System Python 3 interpreter",
}
# Either command satisfies the container runtime check.
CONTAINER_COMMANDS = ("singularity", "apptainer")
MIN_JAVA = 17

# Checked only with --check_tools (i.e. when not relying on the container).
OPTIONAL_COMMANDS = {
    "prefetch": "SRA Toolkit prefetch utility (optional but recommended for SRA downloads)",
}

SRA_INPUT_KEYS = ("rnaseq_sra_single", "rnaseq_sra_paired", "isoseq_sra")


# Command that must be on PATH for each gene finder (checked with --check_tools).
GENEFINDER_CLI: Dict[str, str] = {
    "tiberius": "tiberius.py",
    "vipsania": "vipsania",
}

# Git submodules of this repository (directory names relative to the root):
# the gene finders, and Drusilla, the ORF annotator for assembled transcripts.
SUBMODULES: Dict[str, str] = {
    "tiberius": "tiberius",
    "vipsania": "vipsania",
    "drusilla": "drusilla",
}

# Params keys that the command line can set, by params block.
TOP_LEVEL_CLI_KEYS = (
    "threads", "outdir", "genome", "proteins", "odb12Partitions",
    "rnaseq_single", "rnaseq_paired", "rnaseq_bam", "rnaseq_varus",
    "rnaseq_sra_single", "rnaseq_sra_paired", "isoseq", "isoseq_sra", "isoseq_varus", "mixed_varus",
    "mode", "scoring_matrix",
)
GENEFINDER_CLI_KEYS: Dict[str, Tuple[str, ...]] = {
    "tiberius": ("model_cfg", "model_dir", "result", "min_split_size", "max_files", "max_parallel", "batch_size",
                 "seq_len"),
    "vipsania": ("model", "model_dir", "result", "min_split_size", "max_files", "max_parallel", "batch_size",
                 "context", "finetune", "finetune_epochs"),
}

# Params values that are paths: they are written as absolute paths, because
# Nextflow expands neither '~' nor environment variables.
PATH_KEYS = ("genome", "proteins", "rnaseq_single", "rnaseq_paired", "rnaseq_bam", "isoseq", "scoring_matrix",
             "rnaseq_varus", "isoseq_varus", "mixed_varus")
GENEFINDER_PATH_KEYS = ("result", "model_cfg", "model_dir")

ROOT_ENV = "PALUDAMENTUM_ROOT"


def _default_repo_root() -> Path:
    """
    Root of the Paludamentum checkout (the directory that holds main.nf):
    ``$PALUDAMENTUM_ROOT``, else the checkout that this module is imported
    from (``pip install -e .``, ``python -m paludamentum`` in the checkout).
    A plain ``pip install .`` copies only the launcher into site-packages, so
    the pipeline must then be named by the environment variable.
    """
    env_root = os.environ.get(ROOT_ENV)
    if env_root:
        root = Path(env_root).expanduser().resolve()
        if not (root / "main.nf").is_file():
            raise SystemExit(f"{ROOT_ENV}={env_root} is not a Paludamentum checkout (no main.nf in it).")
        return root
    root = Path(__file__).resolve().parent.parent
    if (root / "main.nf").is_file():
        return root
    raise SystemExit(
        "The Paludamentum pipeline (main.nf, conf/, modules/) was not found next to the installed "
        f"launcher ({root}). Install the launcher from the checkout with 'pip install -e .', or set "
        f"{ROOT_ENV} to the directory of the checkout."
    )


def pipeline_paths(root_override: str | Path | None = None) -> Tuple[Path, Path, Path]:
    """Return (pipeline root, main.nf, conf/base.config)."""
    repo_root = Path(root_override).expanduser().resolve() if root_override else _default_repo_root()
    pipeline_main = repo_root / "main.nf"
    base_config = repo_root / "conf" / "base.config"
    return repo_root, pipeline_main, base_config


def resolve_nf_config(value: str | Path, root: str | Path | None = None) -> Path:
    """
    Resolve a Nextflow config given on the command line.

    Accepts an existing path, or the name of a config shipped in ``conf/``
    (``conf/slurm_generic.config``, ``slurm_generic.config`` or ``slurm_generic``).
    """
    repo_root, _, _ = pipeline_paths(root)
    conf_dir = repo_root / "conf"
    raw = Path(os.path.expandvars(str(value))).expanduser()
    if raw.is_file():
        return raw.resolve()
    # A path (with a directory part) that does not exist is an error; only
    # the shorthand conf/<name> and a bare name are looked up in conf/.
    if raw.parent != Path(".") and raw.parent != Path("conf"):
        raise SystemExit(f"Nextflow config not found: {raw}")
    candidates = [conf_dir / raw.name]
    if raw.suffix != ".config":
        candidates.append(conf_dir / f"{raw.name}.config")
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    shipped = ", ".join(sorted(p.stem for p in conf_dir.glob("*.config")))
    raise SystemExit(f"Nextflow config not found: {value} (searched: {raw}, "
                     f"{', '.join(str(c) for c in candidates)}). Shipped configs: {shipped}")


def submodule_root(name: str, root: str | Path | None = None) -> Path:
    """Directory of a submodule (tiberius, vipsania, drusilla) in the checkout."""
    if name not in SUBMODULES:
        raise SystemExit(f"Unknown submodule '{name}'. Known: {', '.join(SUBMODULES)}.")
    repo_root, _, _ = pipeline_paths(root)
    return repo_root / SUBMODULES[name]


def submodule_version(name: str, root: str | Path | None = None) -> str | None:
    """Version in the submodule's pyproject.toml; None if the submodule is not checked out."""
    pyproject = submodule_root(name, root) / "pyproject.toml"
    if not pyproject.is_file():
        return None
    match = re.search(r'^version\s*=\s*"([^"]+)"', pyproject.read_text(encoding="utf-8"), re.MULTILINE)
    return match.group(1) if match else None


def container_tags(base_config: Path) -> Dict[str, str]:
    """Image tags pinned in base.config, keyed by image name (tiberius, vipsania, ...)."""
    tags: Dict[str, str] = {}
    pattern = re.compile(r"""container\s*=\s*["']docker://([^:"']+):([^"']+)["']""")
    for match in pattern.finditer(base_config.read_text(encoding="utf-8")):
        tags[match.group(1).rsplit("/", 1)[-1]] = match.group(2)
    return tags


def version_mismatches(root: str | Path | None = None) -> List[str]:
    """
    Warnings for submodules (gene finders and Drusilla) whose version differs
    from the image tag that conf/base.config runs. The submodule pins the
    version; the image of the same version must be pinned in base.config.
    """
    repo_root, _, base_config = pipeline_paths(root)
    tags = container_tags(base_config)
    problems = []
    for name in SUBMODULES:
        version = submodule_version(name, repo_root)
        tag = tags.get(name)
        if version and tag and version != tag:
            problems.append(
                f"The {name} submodule is version {version}, but conf/base.config runs the image "
                f"{name}:{tag}. Bump the submodule or the image tag."
            )
    return problems


def tiberius_checkout(root: str | Path | None = None) -> Path | None:
    """Directory of the Tiberius submodule if it holds tiberius.py, else None."""
    directory = submodule_root("tiberius", root)
    return directory if (directory / "tiberius.py").is_file() else None


def tiberius_bin_dir(root: str | Path | None = None) -> Path | None:
    """
    Directory that makes ``tiberius.py`` of the checkout callable, or None
    without a checkout. ``tiberius.py`` is not executable in the Tiberius
    repository, so a wrapper is written to ``~/.cache/paludamentum/bin`` (or
    ``$XDG_CACHE_HOME/paludamentum/bin``) that runs it with python3.
    """
    checkout = tiberius_checkout(root)
    if not checkout:
        return None
    script = checkout / "tiberius.py"
    if os.access(script, os.X_OK):
        return checkout
    cache = Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache") / "paludamentum" / "bin"
    wrapper = cache / "tiberius.py"
    body = f'#!/bin/sh\nexec python3 "{script}" "$@"\n'
    try:
        cache.mkdir(parents=True, exist_ok=True)
        if not wrapper.is_file() or wrapper.read_text(encoding="utf-8") != body:
            wrapper.write_text(body, encoding="utf-8")
        wrapper.chmod(0o755)
    except OSError as exc:
        print(f"[WARN] Cannot write the tiberius.py wrapper {wrapper}: {exc}")
        return None
    return cache


def nextflow_env(root: str | Path | None = None) -> Dict[str, str]:
    """
    Environment for the Nextflow process. The directory with the callable
    tiberius.py of the checkout is appended to PATH, so that runs without
    containers find it; an installed Tiberius earlier on PATH wins.
    """
    env = os.environ.copy()
    bin_dir = tiberius_bin_dir(root)
    if bin_dir:
        env["PATH"] = os.pathsep.join(part for part in (env.get("PATH", ""), str(bin_dir)) if part)
    return env


def resolve_model_cfg(value: str | Path, root: str | Path | None = None) -> Path:
    """
    Resolve a Tiberius model configuration to a file, because the Nextflow
    process stages it. Accepts a path, or a name such as ``diatoms`` or
    ``diatoms.yaml`` that is looked up in ``model_cfg/`` of the Tiberius submodule.
    """
    candidate = Path(os.path.expandvars(str(value))).expanduser()
    if candidate.is_file():
        return candidate.resolve()
    cfg_dir = submodule_root("tiberius", root) / "model_cfg"
    stem = candidate.name
    for suffix in (".yaml", ".yml"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
    searched = [candidate]
    for directory in (cfg_dir, cfg_dir / "superseded"):
        for name in dict.fromkeys((candidate.name, f"{stem}.yaml", f"{stem}.yml")):
            path = directory / name
            searched.append(path)
            if path.is_file():
                if directory.name == "superseded":
                    print(f"[WARN] The Tiberius model configuration '{value}' is superseded; "
                          "a newer model may be available.")
                return path.resolve()
    if not cfg_dir.is_dir():
        raise SystemExit(
            f"Tiberius model configuration '{value}' is not a file, and the Tiberius submodule "
            f"is not checked out to look it up by name ({cfg_dir}). Run\n"
            "    git submodule update --init tiberius\n"
            "or give the path of a model configuration file."
        )
    listing = ", ".join(sorted(p.stem for p in cfg_dir.glob("*.y*ml")))
    raise SystemExit(
        f"Tiberius model configuration not found: {value} (searched: "
        f"{', '.join(str(p) for p in searched)}). Available names: {listing}"
    )


def default_params(genefinder: str = "tiberius", root: str | Path | None = None) -> Dict:
    """
    Launcher-side defaults of the params file for a gene finder command line.
    The selected gene finder is switched on; conf/base.config holds all other defaults.
    """
    if genefinder not in GENEFINDER_CLI:
        raise SystemExit(f"Unknown gene finder '{genefinder}'. Supported: {', '.join(GENEFINDER_CLI)}.")
    repo_root, _, _ = pipeline_paths(root)
    params: Dict = {
        "threads": 48,
        "outdir": f"{genefinder}_results",
        "genome": None,
        "proteins": None,
        "odb12Partitions": [],
        "rnaseq_sra_single": [],
        "rnaseq_sra_paired": [],
        "isoseq_sra": [],
        "rnaseq_single": [],
        "rnaseq_paired": [],
        "isoseq": [],
        "genefinder": genefinder,
        "mode": None,
        "scoring_matrix": str(repo_root / "conf" / "blosum62.csv"),
    }
    for name in GENEFINDER_CLI:
        params[name] = {"run": name == genefinder}
    return params


def merge_params(base: Dict, overrides: Dict) -> Dict:
    """Recursively merge ``overrides`` into ``base`` (in place). None values do not override."""
    for key, value in overrides.items():
        if value is None:
            continue
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            merge_params(base[key], value)
        else:
            base[key] = value
    return base


def resolve_outdir(params: Dict, outdir: str | Path | None = None) -> Path:
    """Absolute output directory of a run; nothing is created."""
    out = Path(outdir or params.get("outdir") or "results").expanduser()
    if not out.is_absolute():
        out = (Path.cwd() / out).resolve()
    return out


def write_params_yaml(params: Dict, outdir: str | Path | None = None) -> Path:
    """Write the merged params to ``<outdir>/params.yaml`` and return the path."""
    out = resolve_outdir(params, outdir)
    out.mkdir(parents=True, exist_ok=True)
    params["outdir"] = str(out)
    params_path = out / "params.yaml"
    with params_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(params, handle, sort_keys=False)
    return params_path


def infer_genefinder(args, loaded: Dict) -> str:
    """
    Gene finder of a run: ``--genefinder``, else ``genefinder`` in the params
    file, else the single block with ``run: true``, else the gene finder whose
    model option was given on the command line, else Tiberius. Two blocks
    with ``run: true`` and no explicit choice are an error.
    """
    explicit = getattr(args, "genefinder", None)
    if explicit:
        return str(explicit).lower()
    if loaded.get("genefinder"):
        return str(loaded["genefinder"]).lower()
    enabled = [
        name for name in GENEFINDER_CLI
        if isinstance(loaded.get(name), dict) and loaded[name].get("run")
    ]
    if len(enabled) > 1:
        raise SystemExit(
            f"More than one gene finder has run: true in the params file ({', '.join(enabled)}). "
            "Set 'genefinder' in the file or use --genefinder."
        )
    if len(enabled) == 1:
        return enabled[0]
    if getattr(args, "model", None) and not getattr(args, "model_cfg", None):
        return "vipsania"
    return "tiberius"


def cli_overrides(args, genefinder: str) -> Dict:
    """Params values given on the command line; unset options do not override."""
    overrides: Dict = {}

    def value_of(key: str):
        value = getattr(args, key, None)
        # unset: None, an empty string or list, or a flag that was not given.
        # 0 is a value (--threads 0 is rejected by the validation later).
        if value is None or value is False or (isinstance(value, (str, list)) and not value):
            return None
        return value

    for key in TOP_LEVEL_CLI_KEYS:
        value = value_of(key)
        if value is not None:
            overrides[key] = value
    # A single command line value for paired reads is a glob that covers all
    # libraries; the pipeline expects a string then, not a list.
    paired = overrides.get("rnaseq_paired")
    if isinstance(paired, list) and len(paired) == 1:
        overrides["rnaseq_paired"] = paired[0]
    elif isinstance(paired, list) and len(paired) > 2:
        raise SystemExit(
            "--rnaseq_paired takes one glob for all libraries (\"RNA/*_{1,2}.fastq.gz\") or the two "
            f"files of one library, not {len(paired)} files. Put several libraries as [r1, r2] pairs "
            "into the params file (rnaseq_paired)."
        )

    block: Dict = {}
    for key in GENEFINDER_CLI_KEYS.get(genefinder, ()):
        value = value_of(key)
        if value is not None:
            block[key] = value
    if block:
        overrides[genefinder] = block
    return overrides


def merge_run_params(args, genefinder: str | None = None, root: str | Path | None = None) -> Tuple[str, Dict]:
    """
    The params of a run: launcher defaults, then the params file
    (``args.params_yaml``), then command line values. Nothing is written.
    Returns the gene finder and the merged params.
    """
    loaded: Dict = {}
    params_yaml = getattr(args, "params_yaml", None)
    if params_yaml:
        params_path = Path(params_yaml).expanduser().resolve()
        if not params_path.exists():
            raise SystemExit(f"Params YAML not found: {params_path}")
        loaded = load_params(params_path)

    explicit = bool(getattr(args, "genefinder", None))
    genefinder = genefinder or infer_genefinder(args, loaded)
    if genefinder not in GENEFINDER_CLI:
        raise SystemExit(f"Unknown gene finder '{genefinder}'. Supported: {', '.join(GENEFINDER_CLI)}.")
    params = default_params(genefinder, root)
    merge_params(params, loaded)
    merge_params(params, cli_overrides(args, genefinder))

    # The launcher's choice is the pipeline's choice: one gene finder runs,
    # the others are off, whatever the params file says (--genefinder wins).
    params["genefinder"] = genefinder
    for name in GENEFINDER_CLI:
        if not isinstance(params.get(name), dict):
            params[name] = {}
        if name != genefinder:
            params[name]["run"] = False
    cfg = params[genefinder]
    if explicit and not cfg.get("run"):
        print(f"[INFO] --genefinder {genefinder} switches {genefinder}.run on (the params file had run: false).")
        cfg["run"] = True

    if not params.get("genome"):
        raise SystemExit("A genome is required: --genome, or 'genome' in the params file.")

    for key in PATH_KEYS:
        if params.get(key):
            params[key] = absolute_paths(params[key], key)
    for key in GENEFINDER_PATH_KEYS:
        if cfg.get(key):
            cfg[key] = absolute_paths(cfg[key], f"{genefinder}.{key}")

    if cfg.get("run") and not cfg.get("result"):
        if genefinder == "tiberius":
            if not cfg.get("model_cfg"):
                raise SystemExit("Tiberius needs a model configuration: --model_cfg, or tiberius.model_cfg in the params file.")
            cfg["model_cfg"] = str(resolve_model_cfg(cfg["model_cfg"], root))
        elif genefinder == "vipsania" and not cfg.get("model"):
            raise SystemExit("Vipsania needs a model: --model, or vipsania.model in the params file.")

    return genefinder, params


def build_params(args, genefinder: str | None = None, root: str | Path | None = None) -> Tuple[str, Path]:
    """
    Merge the params of a run (see ``merge_run_params``) and write them to
    ``<outdir>/params.yaml``. Returns the gene finder and the written path.
    """
    genefinder, params = merge_run_params(args, genefinder, root)
    return genefinder, write_params_yaml(params)


def absolute_paths(value, key: str):
    """
    ``value`` (a path string, a glob, or a list of those, nested for read
    pairs) with '~' and environment variables expanded and relative paths
    made absolute against the launch directory, as Nextflow resolves them.
    """
    if isinstance(value, list):
        return [absolute_paths(v, key) for v in value]
    if not isinstance(value, (str, Path)):
        raise SystemExit(f"'{key}' must be a path or a list of paths, got {type(value).__name__}: {value!r}")
    expanded = os.path.expandvars(os.path.expanduser(str(value).strip()))
    if not expanded:
        return expanded
    path = Path(expanded)
    if not path.is_absolute():
        path = Path.cwd() / path
    return os.path.normpath(str(path))


def load_params(params_path: Path) -> Dict:
    with params_path.open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise SystemExit(f"Expected a mapping at the top-level of {params_path}, got {type(data).__name__}.")
    return data


def has_glob_char(value: str) -> bool:
    return any(char in value for char in GLOB_CHARS)


def expand_braces(pattern: str) -> List[str]:
    start = pattern.find("{")
    if start == -1:
        return [pattern]
    end = pattern.find("}", start)
    if end == -1:
        return [pattern]
    prefix = pattern[:start]
    suffix = pattern[end + 1 :]
    choices = pattern[start + 1 : end].split(",")
    expanded = []
    for choice in choices:
        expanded.extend(expand_braces(f"{prefix}{choice}{suffix}"))
    return expanded


def iter_strings(value, key: str = "") -> Iterable[str]:
    if isinstance(value, (list, tuple, set)):
        for sub in value:
            yield from iter_strings(sub, key)
    elif isinstance(value, (str, Path)):
        yield str(value)
    elif value is None:
        return
    else:
        raise SystemExit(f"'{key}' in the params file must be a path or a list of paths, "
                         f"got {type(value).__name__}: {value!r}")


def resolve_data_entries(value, key: str = "") -> Tuple[List[Path], List[str]]:
    resolved: List[Path] = []
    errors: List[str] = []
    for raw in iter_strings(value, key):
        cleaned = raw.strip()
        if not cleaned:
            continue
        expanded = os.path.expandvars(os.path.expanduser(cleaned))
        brace_patterns = expand_braces(expanded) if "{" in cleaned else [expanded]
        for pattern in brace_patterns:
            # Nextflow resolves relative paths against the launch directory, not
            # against the location of the params file. Validate the same way.
            candidate = Path(pattern)
            if not candidate.is_absolute():
                candidate = (Path.cwd() / candidate).resolve()

            if has_glob_char(pattern):
                matches = glob.glob(str(candidate), recursive=True)
                if matches:
                    resolved.extend(Path(match) for match in matches)
                else:
                    errors.append(f"No files matched pattern '{cleaned}'.")
            else:
                resolved.append(candidate)
    return resolved, errors


def validate_input_data(params: Dict, params_path: Path) -> List[str]:
    """Problems with the input files named in ``params`` (empty list: none)."""
    errors: List[str] = []

    def check_entries(value, label: str, key: str, directory: bool = False) -> None:
        files, errs = resolve_data_entries(value, key)
        errors.extend(errs)
        for file_path in files:
            text = str(file_path)
            if any(ch.isspace() for ch in text):
                # the shell commands of the processes do not quote file names
                errors.append(f"{label} path contains whitespace, which the pipeline does not support: {text}")
            if not file_path.exists():
                errors.append(f"{label} missing: {file_path}")
            elif directory and not file_path.is_dir():
                errors.append(f"{label} is not a directory: {file_path}")
            elif not directory and not file_path.is_file():
                errors.append(f"{label} is not a file: {file_path}")

    if not params.get("genome"):
        errors.append(f"Genome FASTA ('genome') is not set in {params_path}")
    else:
        check_entries(params["genome"], "Genome FASTA", "genome")

    optional_fields = {
        "rnaseq_single": "RNA-Seq single-end FASTQ",
        "rnaseq_paired": "RNA-Seq paired-end FASTQ",
        "rnaseq_bam": "RNA-Seq BAM",
        "isoseq": "Iso-Seq FASTQ",
        "scoring_matrix": "Scoring matrix",
        "proteins": "Protein FASTA",
    }
    for key, label in optional_fields.items():
        if params.get(key):
            check_entries(params[key], label, key)

    varus_fields = {
        "rnaseq_varus": "pyVARUS directory (rnaseq_varus)",
        "isoseq_varus": "pyVARUS directory (isoseq_varus)",
        "mixed_varus": "pyVARUS directory (mixed_varus)",
    }
    for key, label in varus_fields.items():
        if params.get(key):
            check_entries(params[key], label, key, directory=True)

    for name in GENEFINDER_CLI:
        cfg = params.get(name)
        if not isinstance(cfg, dict):
            continue
        # A mistyped result must not silently start the GPU gene finder.
        if cfg.get("result"):
            check_entries(cfg["result"], f"{name.capitalize()} result", f"{name}.result")
        if cfg.get("model_dir"):
            check_entries(cfg["model_dir"], f"{name.capitalize()} model_dir", f"{name}.model_dir", directory=True)
        if name == "tiberius" and cfg.get("run") and cfg.get("model_cfg") and not cfg.get("result"):
            check_entries(cfg["model_cfg"], "Tiberius model_cfg", "tiberius.model_cfg")

    if params.get("threads") is not None and (not isinstance(params["threads"], int) or params["threads"] < 1):
        errors.append(f"'threads' must be a positive integer, got {params['threads']!r}")

    return errors


def resolve_executable(command: str, relative_to: Path) -> Path | None:
    expanded = os.path.expandvars(os.path.expanduser(command))
    has_sep = os.sep in expanded or (os.altsep and os.altsep in expanded)
    if has_sep:
        candidate = Path(expanded)
        if not candidate.is_absolute():
            candidate = (relative_to / candidate).resolve()
        if candidate.exists() and os.access(candidate, os.X_OK):
            return candidate
        return None
    found = shutil.which(expanded)
    return Path(found).resolve() if found else None


def check_java_version(java_path: Path) -> Tuple[bool, str | None]:
    try:
        proc = subprocess.run(
            [str(java_path), "-version"],
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError as exc:
        return False, str(exc)
    output = proc.stderr or proc.stdout
    version_line = output.splitlines()[0] if output else ""
    marker = '"'
    if marker in version_line:
        version = version_line.split(marker)[1]
        major = version.split(".")[0]
        try:
            if int(major) >= MIN_JAVA:
                return True, None
        except ValueError:
            pass
        return False, f"Java version {version} detected, but {MIN_JAVA}+ is required by Nextflow."
    return False, "Unable to parse Java version output."


def build_tool_command_map(params: Dict) -> Dict[str, str]:
    overrides = params.get("tools") or {}
    command_map = DEFAULT_TOOL_BINARIES.copy()
    if isinstance(overrides, dict):
        for key, value in overrides.items():
            if key in command_map and value:
                command_map[key] = str(value)
    return command_map


def validate_executables(
    params: Dict,
    nextflow_bin: str,
    check_tool_binaries: bool,
    skip_singularity_check: bool,
    repo_root: Path,
) -> Tuple[List[str], List[str]]:
    errors: List[str] = []
    warnings: List[str] = []

    def check_command(command: str, description: str, mandatory: bool = True, java: bool = False) -> None:
        path = resolve_executable(command, repo_root)
        if not path:
            msg = f"{description} not found on PATH (looked for '{command}')."
            if mandatory:
                errors.append(msg)
            else:
                warnings.append(msg)
            return
        if java:
            ok, problem = check_java_version(path)
            if not ok:
                errors.append(problem or f"Unable to validate Java executable at {path}.")

    for cmd, desc in GENERAL_COMMANDS.items():
        if cmd == "nextflow":
            check_command(nextflow_bin, desc)
        elif cmd == "singularity":
            if skip_singularity_check:
                continue
            if not any(resolve_executable(c, repo_root) for c in CONTAINER_COMMANDS):
                errors.append(f"{desc} not found on PATH (looked for {' or '.join(CONTAINER_COMMANDS)}). "
                              "Use --skip_singularity_check to run without containers.")
        elif cmd == "java":
            check_command(cmd, desc, java=True)
        else:
            check_command(cmd, desc)

    if check_tool_binaries:
        for cmd, desc in OPTIONAL_COMMANDS.items():
            check_command(cmd, desc, mandatory=False)
        if any(params.get(k) for k in SRA_INPUT_KEYS):
            check_command("fasterq-dump", "SRA Toolkit fasterq-dump utility")
        if params.get("proteins"):
            check_command("perl", "Perl interpreter (needed for aln2hints.pl)")
        tools = build_tool_command_map(params)
        for key, cmd in tools.items():
            label = TOOL_DESCRIPTIONS.get(key, f"Tool '{key}'")
            check_command(cmd, label)

        for name, cli in GENEFINDER_CLI.items():
            finder_cfg = params.get(name) or {}
            if not (isinstance(finder_cfg, dict) and finder_cfg.get("run")) or finder_cfg.get("result"):
                continue
            if name == "tiberius" and not resolve_executable(cli, repo_root):
                # tiberius.py of the checkout is appended to PATH for the Nextflow process
                bin_dir = tiberius_bin_dir(repo_root)
                if bin_dir and os.access(bin_dir / cli, os.X_OK):
                    continue
            check_command(cli, f"{name.capitalize()} CLI ({cli})")

    return errors, warnings


def run_nextflow(
    params_path: Path,
    config_path: Path,
    base_config_path: Path,
    pipeline_main_path: Path,
    profile: str | None,
    nextflow_bin: str,
    resume: bool,
    work_dir: str | None,
    extra_args: Sequence[str],
) -> int:
    launch_cwd = Path.cwd()
    cmd = [
        nextflow_bin,
        "run",
        str(pipeline_main_path),
        "-params-file",
        str(params_path),
        "-c",
        str(base_config_path),
    ]
    # The user config is layered on top of base.config. Skip it if it is base.config itself.
    if Path(config_path).resolve() != Path(base_config_path).resolve():
        cmd.extend(["-c", str(config_path)])
    if profile:
        cmd.extend(["-profile", profile])
    if resume:
        cmd.append("-resume")
    if work_dir:
        cmd.extend(["-work-dir", work_dir])
    if extra_args:
        cmd.extend(extra_args)

    print("[INFO] Launching Nextflow with command:")
    print("       " + " ".join(shlex.quote(part) for part in cmd))

    completed = subprocess.run(cmd, cwd=launch_cwd, env=nextflow_env(Path(pipeline_main_path).parent))
    return completed.returncode


def run_nextflow_pipeline(
    args,
    genefinder: str = "tiberius",
    pipeline_root: str | Path | None = None,
    params: Dict | None = None,
) -> None:
    """
    Validate and launch the pipeline.

    ``args`` is an argparse-like namespace. Required attributes: ``params_yaml``
    (a complete params file, see ``build_params``), ``nf_config``. Optional:
    ``nextflow_args``, ``work_dir``, ``profile``, ``nextflow_bin``, ``resume``,
    ``check_tools``, ``skip_singularity_check``, ``dry_run``. ``genefinder`` is
    the gene finder of the run; the pipeline itself selects it from the params.

    ``params`` are merged params (see ``merge_run_params``) instead of
    ``args.params_yaml``: they are validated first and written to
    ``<outdir>/params.yaml`` only when the validation passes, so a failed
    run leaves no params file behind. ``args.params_yaml`` is set to it.
    """
    if genefinder not in GENEFINDER_CLI:
        raise SystemExit(f"Unknown gene finder '{genefinder}'. Supported: {', '.join(GENEFINDER_CLI)}.")

    params_yaml = getattr(args, "params_yaml", None)
    config = getattr(args, "nf_config", None)
    if (params is None and not params_yaml) or not config:
        raise SystemExit("Launching Nextflow requires --params_yaml and --nf_config.")

    extra_args = list(getattr(args, "nextflow_args", None) or [])
    if extra_args and extra_args[0] == "--":
        extra_args = extra_args[1:]

    if params is None:
        params_path = Path(params_yaml).expanduser().resolve()
    else:
        params_path = resolve_outdir(params) / "params.yaml"
    repo_root, pipeline_main, base_config = pipeline_paths(pipeline_root)
    config_path = resolve_nf_config(config, repo_root)
    work_dir: Path | None = None
    if getattr(args, "work_dir", None):
        work_dir = Path(args.work_dir).expanduser().resolve()

    if params is None and not params_path.exists():
        raise SystemExit(f"Params YAML not found: {params_path}")
    if not pipeline_main.exists():
        raise SystemExit(f"Pipeline entry point missing: {pipeline_main}")
    if not base_config.exists():
        raise SystemExit(f"Base config not found: {base_config}")

    from_file = params is None
    if from_file:
        params = load_params(params_path)
    nextflow_bin = getattr(args, "nextflow_bin", None) or "nextflow"

    print("[INFO] Validating input files...")
    data_errors = validate_input_data(params, params_path)
    if data_errors:
        for error in data_errors:
            print(f"[ERROR] {error}")
        raise SystemExit("Input validation failed.")

    print("[INFO] Validating required executables...")
    exec_errors, exec_warnings = validate_executables(
        params=params,
        nextflow_bin=nextflow_bin,
        check_tool_binaries=bool(getattr(args, "check_tools", False)),
        skip_singularity_check=bool(getattr(args, "skip_singularity_check", False)),
        repo_root=repo_root,
    )
    if exec_warnings:
        for warning in exec_warnings:
            print(f"[WARN] {warning}")
    if exec_errors:
        for error in exec_errors:
            print(f"[ERROR] {error}")
        raise SystemExit("Executable validation failed.")

    for warning in version_mismatches(repo_root):
        print(f"[WARN] {warning}")

    if not from_file:
        params_path = write_params_yaml(params)
        args.params_yaml = str(params_path)
        print(f"[INFO] Params written to {params_path}")

    if getattr(args, "dry_run", False):
        print("[INFO] Dry run requested; skipping Nextflow execution.")
        return

    returncode = run_nextflow(
        params_path=params_path,
        config_path=config_path,
        base_config_path=base_config,
        pipeline_main_path=pipeline_main,
        profile=getattr(args, "profile", None),
        nextflow_bin=nextflow_bin,
        resume=bool(getattr(args, "resume", False)),
        work_dir=str(work_dir) if work_dir else None,
        extra_args=extra_args,
    )
    if returncode != 0:
        raise SystemExit(returncode)
