"""Static checks of the Nextflow modules (no Nextflow needed)."""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MODULES = sorted((ROOT / "modules").glob("*.nf"))


def script_lines():
    for path in MODULES:
        for number, line in enumerate(path.read_text().splitlines(), 1):
            if not line.strip().startswith("//"):
                yield path.name, number, line


def test_tools_use_the_reserved_cpus():
    """Thread options must follow task.cpus, not params.threads: otherwise a task
    reserves N CPUs on the cluster and runs a different number of threads."""
    thread_option = re.compile(r"(--threads|-t|-p|-@)\s+\$\{params\.threads\}")
    offenders = [f"{name}:{n}" for name, n, line in script_lines() if thread_option.search(line)]
    assert not offenders, offenders


def test_every_diamond_call_sets_threads():
    """Without --threads DIAMOND uses every core of the node."""
    missing = []
    for path in MODULES:
        text = path.read_text()
        for match in re.finditer(r"diamond\}?\s+(makedb|blastp)", text):
            # the command may continue over several lines ending in a backslash
            end = match.end()
            while True:
                nl = text.find("\n", end)
                if nl == -1 or not text[end:nl].rstrip().endswith("\\"):
                    break
                end = nl + 1
            command = text[match.start(): nl if nl != -1 else len(text)]
            if "--threads ${task.cpus}" not in command:
                missing.append(f"{path.name}: {command.splitlines()[0].strip()}")
    assert not missing, missing


def test_base_config_reserves_params_threads_by_default():
    text = (ROOT / "conf" / "base.config").read_text()
    assert re.search(r"^\s*cpus\s*=\s*\{\s*params\.threads", text, re.MULTILINE)
