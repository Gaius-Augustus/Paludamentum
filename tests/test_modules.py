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


POSTPROCESS_MODULES = ("finalize.nf", "qc.nf", "completeness.nf", "ncrna.nf", "fantasia.nf")


def base_config_cpus_by_name():
    """Process names of the withName selectors of base.config that set cpus."""
    text = (ROOT / "conf" / "base.config").read_text()
    names = set()
    for match in re.finditer(r"withName:\s*'([^']+)'\s*\{([^}]*)\}", text):
        if re.search(r"^\s*cpus\s*=", match.group(2), re.MULTILINE):
            names.update(match.group(1).split("|"))
    return names


def test_postprocess_tasks_do_not_reserve_the_default_cpus():
    """A post-processing task whose tools do not run with task.cpus threads
    sets its cpus in base.config; otherwise it reserves the 48 of the default
    (a single-threaded report on 48 CPUs, FANTASIA with 48 CPUs on a GPU node)."""
    with_cpus = base_config_cpus_by_name()
    missing = []
    for name in POSTPROCESS_MODULES:
        text = (ROOT / "modules" / name).read_text()
        for match in re.finditer(r"^process\s+(\w+)\s*\{(.*?)^\}", text, re.MULTILINE | re.DOTALL):
            process, body = match.groups()
            if "label 'download'" in body or "${task.cpus}" in body:
                continue
            if process not in with_cpus:
                missing.append(f"{name}: {process}")
    assert not missing, missing


def base_config_tool_defaults():
    """{params.tools key: default executable name} from the tools block of base.config."""
    text = (ROOT / "conf" / "base.config").read_text()
    block = re.search(r"^  tools \{(.*?)^  \}", text, re.MULTILINE | re.DOTALL)
    assert block, "conf/base.config has no tools block"
    return dict(re.findall(r"^\s*(\w+)\s*=\s*\"([^\"]+)\"", block.group(1), re.MULTILINE))


def test_every_configured_tool_is_called_through_its_override():
    """A tool with a params.tools key must be called as ${params.tools.<key>}:
    a bare name in a script block silently ignores the override (which
    --check_tools does check and docs/parameters.md does document)."""
    defaults = base_config_tool_defaults()
    assert "gffread" in defaults and "omark" in defaults, sorted(defaults)
    # the executable at command position: start of line, or after a pipe, a
    # command separator, a $( or the ) of a case branch
    patterns = {
        key: re.compile(r"(?:^|[|;&)]|\$\()\s*" + re.escape(name) + r"(?=\s)")
        for key, name in defaults.items()
    }
    offenders = []
    for name, number, line in script_lines():
        for key, pattern in patterns.items():
            if pattern.search(line):
                offenders.append(f"{name}:{number}: {defaults[key]} (params.tools.{key})")
    assert not offenders, offenders
