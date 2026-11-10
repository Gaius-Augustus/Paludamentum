"""PREPROCESS_PROTEINDB must hand DIAMOND and miniprot a clean FASTA.

A user's Viridiplantae_train.fa started with a blank line and `diamond makedb`
failed ("Error detecting input file format. First line seems to be blank").
The process now sanitizes the FASTA: it strips CR, drops blank and
whitespace-only lines and keeps only the first token of each header.

The unit tests render the process script the way Groovy renders a triple-quoted
GString, cut out the `sanitize_fasta` shell function and the species-selection
awk (the > 1M OrthoDB path) and run them through bash on toy FASTAs. They need
bash and awk only. The Nextflow test runs the process for real on the
<= 1M proteins path (no diamond); it needs a ``nextflow`` executable (or
NEXTFLOW_BIN), like tests/test_no_writes_to_inputs.py.
"""
from __future__ import annotations

import hashlib
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
PROTEINS_NF = ROOT / "modules" / "proteins.nf"
NEXTFLOW = os.environ.get("NEXTFLOW_BIN") or shutil.which("nextflow")
needs_nextflow = pytest.mark.skipif(NEXTFLOW is None, reason="nextflow executable not found")
needs_bash_awk = pytest.mark.skipif(
    shutil.which("bash") is None or shutil.which("awk") is None, reason="bash and awk are required")

GROOVY_ESCAPES = {"\\": "\\", "$": "$", "t": "\t", "n": "\n", "r": "\r", '"': '"', "'": "'"}


def script_block(text: str) -> str:
    """The raw text between the triple quotes of PREPROCESS_PROTEINDB's script:."""
    block = text[text.index("process PREPROCESS_PROTEINDB"):]
    after = block[block.index("script:") + len("script:"):]   # Groovy statements (def top ...) precede the GString
    body = after[after.index('"""') + 3:]
    return body[:body.index('"""')]


def groovy_unescape(s: str) -> str:
    out, i = [], 0
    while i < len(s):
        if s[i] == "\\" and i + 1 < len(s):
            out.append(GROOVY_ESCAPES.get(s[i + 1], s[i:i + 2]))
            i += 2
        else:
            out.append(s[i])
            i += 1
    return "".join(out)


def render(text: str, top: int = 13) -> str:
    """The shell script Nextflow would write to .command.sh."""
    rendered = groovy_unescape(script_block(text))
    for placeholder, value in (("${proteinDB}", "db/in.fa"), ("${task.cpus}", "1"), ("${top}", str(top)),
                               ("${params.tools.diamond}", "diamond")):
        rendered = rendered.replace(placeholder, value)
    return rendered


def selection_condition(rendered: str) -> str:
    """The `[[ ... ]]` test of the `if [[ top -eq 0 || ... ]]; then` line that picks the full DB."""
    m = re.search(r"^\s*if (\[\[ \d+ -eq 0 \|\|.*\]\]); then\s*$", rendered, re.M)
    assert m, "selection condition line not found"
    return m.group(1)


def sanitize_function(rendered: str) -> str:
    start = rendered.index("sanitize_fasta() {")
    return rendered[start:rendered.index("\n    }\n", start) + len("\n    }\n")]


def species_awk(rendered: str) -> str:
    """The `awk '...' db/in.fa > protein_preprocessed.fa` command of the > 1M OrthoDB path."""
    start = re.search(r"awk '\n\s*BEGIN", rendered).start()
    end = rendered.index("> protein_preprocessed.fa", start) + len("> protein_preprocessed.fa")
    return rendered[start:end]


def run_bash(cwd: Path, code: str) -> str:
    proc = subprocess.run(["bash", "-c", code], cwd=cwd, capture_output=True, timeout=60)
    assert proc.returncode == 0, proc.stderr.decode()
    return proc.stdout.decode()


def stage(tmp_path: Path, fasta: str, top_species: str = "101020\n9606\n") -> Path:
    (tmp_path / "db").mkdir()
    (tmp_path / "db" / "in.fa").write_bytes(fasta.encode())
    (tmp_path / "top_species.txt").write_bytes(top_species.encode())
    return tmp_path


@pytest.fixture(scope="module")
def rendered() -> str:
    return render(PROTEINS_NF.read_text())


# (input FASTA, expected output); ids are OrthoDB-style so both awks apply
CASES = {
    "leading_blank": (
        "\n>101020_0:000001 desc\nMKV\nLLA\n>9606_0:000002\nMAA\n",
        ">101020_0:000001\nMKV\nLLA\n>9606_0:000002\nMAA\n"),
    "blank_and_whitespace_lines": (
        ">101020_0:000001\nMKV\n\n\n>9606_0:000002\nMAA\n   \n \t \nLLA\n\n",
        ">101020_0:000001\nMKV\n>9606_0:000002\nMAA\nLLA\n"),
    "crlf": (
        ">101020_0:000001\r\nMKV\r\n>9606_0:000002\r\nMAA\r\n",
        ">101020_0:000001\nMKV\n>9606_0:000002\nMAA\n"),
    "crlf_with_blank_first_line": (
        "\r\n>101020_0:000001 x\r\nMKV\r\n\r\n>9606_0:000002\r\nMAA\r\n",
        ">101020_0:000001\nMKV\n>9606_0:000002\nMAA\n"),
    "header_description_after_space": (
        ">101020_0:000001 Some protein OS=Arabidopsis\nMKV\n",
        ">101020_0:000001\nMKV\n"),
    "header_description_after_tab": (
        ">101020_0:000001\tOrthoDB description\nMKV\n>9606_0:000002\ttwo\tTabs\nMAA\n",
        ">101020_0:000001\n" "MKV\n>9606_0:000002\nMAA\n"),
    "clean_input_is_unchanged": (
        ">101020_0:000001\nMKV\n>9606_0:000002\nMAA\n",
        ">101020_0:000001\nMKV\n>9606_0:000002\nMAA\n"),
}
CASE_IDS = list(CASES)


@needs_bash_awk
@pytest.mark.parametrize("name", CASE_IDS)
def test_sanitize_fasta(tmp_path: Path, rendered: str, name: str) -> None:
    fasta, expected = CASES[name]
    cwd = stage(tmp_path, fasta)
    got = run_bash(cwd, sanitize_function(rendered) + "\nsanitize_fasta db/in.fa")
    assert got == expected
    assert got.startswith(">") and "\r" not in got and "\n\n" not in got
    assert got.count(">") == fasta.count(">")


@needs_bash_awk
@pytest.mark.parametrize("name", CASE_IDS)
def test_species_selection_awk_keeps_all_species(tmp_path: Path, rendered: str, name: str) -> None:
    fasta, expected = CASES[name]
    cwd = stage(tmp_path, fasta, top_species="101020\n9606\n")
    got = run_bash(cwd, species_awk(rendered) + " && cat protein_preprocessed.fa")
    assert got == expected


@needs_bash_awk
def test_species_selection_awk_drops_unwanted_species(tmp_path: Path, rendered: str) -> None:
    fasta = "\r\n>101020_0:000001 d\r\nMKV\r\n\r\n>9606_0:000002\tx\r\nMAA\r\n>101020_1:000003\r\nLLA\r\n"
    cwd = stage(tmp_path, fasta, top_species="101020\n")
    got = run_bash(cwd, species_awk(rendered) + " && cat protein_preprocessed.fa")
    assert got == ">101020_0:000001\nMKV\n>101020_1:000003\nLLA\n"


GS_FASTA = (">Acorus_americanus_Acora_01G000100_1_v1_1 desc\nMKV\n"
            ">Homo_sapiens_ENSP000001\nMAA\nLLA\n"
            ">Nannochloris_sp._g12\nGGG\n"
            ">Acorus_americanus_Acora_02G000200_1\nPPP\n")


@needs_bash_awk
def test_species_selection_awk_genus_species_database(tmp_path: Path, rendered: str) -> None:
    cwd = stage(tmp_path, GS_FASTA, top_species="Acorus_americanus\nNannochloris_sp.\n")
    got = run_bash(cwd, species_awk(rendered) + " && cat protein_preprocessed.fa")
    assert got == (">Acorus_americanus_Acora_01G000100_1_v1_1\nMKV\n"
                   ">Nannochloris_sp._g12\nGGG\n"
                   ">Acorus_americanus_Acora_02G000200_1\nPPP\n")


@needs_bash_awk
def test_species_selection_awk_orthodb_database_unchanged(tmp_path: Path, rendered: str) -> None:
    fasta = ">101020_0:000001 a\nMKV\n>9606_0:000002\nMAA\n>10090_1:000003\nGGG\n"
    cwd = stage(tmp_path, fasta, top_species="9606\n10090\n")
    got = run_bash(cwd, species_awk(rendered) + " && cat protein_preprocessed.fa")
    assert got == ">9606_0:000002\nMAA\n>10090_1:000003\nGGG\n"


SHARED_IDS = [
    "101020_0:000003", "9606_1:0a1b",
    "Acorus_americanus_Acora_01G000100_1_v1_1", "Nannochloris_sp._g12",
    "Canis_lupus_familiaris_X1", "Homo_sapiens_ENSP000001",
    "sp|P12345|NAME_HUMAN", "XP_0123.1", "plainid", "Arabidopsis_thaliana",
]


def python_get_species_id():
    """get_species_id of bin/rank_species_from_diamond.py; the script runs at import, so exec only its definitions."""
    src = (ROOT / "bin" / "rank_species_from_diamond.py").read_text()
    start = src.index("ORTHODB_ID =")
    end = src.index("# best hit")
    ns = {"re": re}
    exec(src[start:end], ns)
    return ns["get_species_id"]


def test_python_species_keys() -> None:
    f = python_get_species_id()
    assert f("101020_0:000003") == "101020"
    assert f("Acorus_americanus_Acora_01G000100_1_v1_1") == "Acorus_americanus"
    assert f("Nannochloris_sp._g12") == "Nannochloris_sp."
    assert f("Canis_lupus_familiaris_X1") == "Canis_lupus"
    for other in ("sp|P12345|NAME_HUMAN", "XP_0123.1"):
        assert f(other) == other


@needs_bash_awk
def test_awk_species_key_agrees_with_python(tmp_path: Path, rendered: str) -> None:
    f = python_get_species_id()
    ids = SHARED_IDS
    # wanted = the python keys of every second id; the awk must keep exactly the records whose key is wanted
    wanted = {f(i) for i in ids[::2]}
    fasta = "".join(f">{i} desc\nM\n" for i in ids)
    cwd = stage(tmp_path, fasta, top_species="".join(w + "\n" for w in sorted(wanted)))
    got = run_bash(cwd, species_awk(rendered) + " && cat protein_preprocessed.fa")
    expected = "".join(f">{i}\nM\n" for i in ids if f(i) in wanted)
    assert got == expected
    assert set(ids[::2]) <= {line[1:] for line in got.splitlines() if line.startswith(">")}


# (top, N_PROT, N_ODB, N_GS, whole DB kept?)
CONDITION_CASES = {
    "top_zero": (0, 5_000_000, 5_000_000, 0, True),
    "top_zero_small_db": (0, 10, 10, 0, True),
    "small_db": (13, 1_000_000, 1_000_000, 0, True),
    "mixed_ids": (13, 2_000_000, 1_000_000, 1_000_000, True),
    "partly_orthodb": (13, 2_000_000, 1_999_999, 0, True),
    "no_species_ids": (13, 2_000_000, 0, 0, True),
    "all_orthodb_large": (13, 1_000_001, 1_000_001, 0, False),
    "all_genus_species_large": (13, 2_000_000, 0, 2_000_000, False),
}


@needs_bash_awk
@pytest.mark.parametrize("name", list(CONDITION_CASES))
def test_selection_condition(tmp_path: Path, name: str) -> None:
    top, n_prot, n_odb, n_gs, full = CONDITION_CASES[name]
    cond = selection_condition(render(PROTEINS_NF.read_text(), top=top))
    out = run_bash(tmp_path, f"N_PROT={n_prot}; N_ODB={n_odb}; N_GS={n_gs}; "
                             f"if {cond}; then echo FULL; else echo SELECT; fi")
    assert out.strip() == ("FULL" if full else "SELECT")


@needs_bash_awk
def test_sanitize_fasta_empty_input_gives_empty_output(tmp_path: Path, rendered: str) -> None:
    cwd = stage(tmp_path, "")
    assert run_bash(cwd, sanitize_function(rendered) + "\nsanitize_fasta db/in.fa") == ""


# --- real Nextflow run ------------------------------------------------------

CONFIG = """process {
  executor = 'local'
  scratch  = true
  cpus     = 1
  memory   = '1 GB'
}
singularity.enabled = false
docker.enabled      = false
params.outdir = '%s'
// the GString interpolates params.tools.diamond even on the full-DB path
params.tools = [diamond: 'diamond']
"""

MAIN_NF = """nextflow.enable.dsl=2
include { PREPROCESS_PROTEINDB } from '%s/modules/proteins.nf'
workflow {
    PREPROCESS_PROTEINDB(file(params.db), file(params.tib))
    PREPROCESS_PROTEINDB.out.db.view { f -> "OUT_FA ${f}" }
}
"""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.mark.nextflow
@needs_nextflow
def test_preprocess_proteindb_cleans_blank_first_line_and_crlf(tmp_path: Path) -> None:
    original = tmp_path / "user_data" / "Viridiplantae_train.fa"
    original.parent.mkdir()
    original.write_bytes(
        b"\r\n>sp1 desc\r\nMKV\r\nLLA\r\n\r\n>sp2\tOrthoDB x\r\nMAA\r\n   \r\n>sp3\r\nGGG\r\n")
    before = sha256(original)
    tib = tmp_path / "user_data" / "tiberius.faa"
    tib.write_bytes(b"")
    (tmp_path / "main.nf").write_text(MAIN_NF % ROOT)
    (tmp_path / "test.config").write_text(CONFIG % (tmp_path / "results"))
    env = {**os.environ, "NXF_ANSI_LOG": "false"}
    proc = subprocess.run(
        [NEXTFLOW, "run", str(tmp_path / "main.nf"), "-c", str(tmp_path / "test.config"),
         "--db", str(original), "--tib", str(tib), "-work-dir", str(tmp_path / "work")],
        cwd=tmp_path, env=env, capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, proc.stdout[-3000:] + proc.stderr[-3000:]
    assert sha256(original) == before, "PREPROCESS_PROTEINDB changed the user's original FASTA"
    out = [line.split(" ", 1)[1] for line in proc.stdout.splitlines() if line.startswith("OUT_FA ")]
    assert len(out) == 1, proc.stdout
    data = Path(out[0]).read_bytes().decode()
    assert data.startswith(">") and "\r" not in data
    assert all(line.strip() for line in data.splitlines()), data
    assert data == ">sp1\nMKV\nLLA\n>sp2\nMAA\n>sp3\nGGG\n"
