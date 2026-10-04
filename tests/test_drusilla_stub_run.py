"""Nextflow tests of the Drusilla flow beyond test_stub_run.py.

Stub runs of the model names that select the flow and of a reused
tiberius.result. The options shards, fix_stop/fix_start and lgb_keep change
only the script blocks, which a stub run does not render: those processes run
their real scripts here, with fake drusilla and bin/ scripts.

Needs a ``nextflow`` executable (or NEXTFLOW_BIN), like test_stub_run.py.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest
import yaml

from test_stub_run import (DATA, EVIDENCE, NEXTFLOW, ROOT, assert_ok, drusilla_params,
                           run_pipeline)

pytestmark = [
    pytest.mark.nextflow,
    pytest.mark.skipif(NEXTFLOW is None, reason="nextflow executable not found"),
]


# ---- stub runs -----------------------------------------------------------------

def model_params(tmp_path: Path, tool: str, model: str) -> dict:
    """drusilla_params with another gene finder model."""
    params = drusilla_params(tmp_path, tool)
    if tool == "tiberius":
        # the flow reads target_species of the model, as Tiberius' model_cfg/ has it
        stem = model.split(".")[0]
        species = {"vertebrates": "Vertebrata"}.get(stem, "Mammalia" if stem.startswith("mammalia") else stem.title())
        cfg = tmp_path / model
        cfg.write_text(f'target_species: "{species}"\n')
        params["tiberius"]["model_cfg"] = str(cfg)
    else:
        params["vipsania"]["model"] = model
    return params


@pytest.mark.parametrize("tool,model,method", [
    ("tiberius", "mammalia_softmasking_v2.yaml", "drusilla"),
    ("tiberius", "vertebrates.yml", "drusilla"),
    ("tiberius", "eukaryota.yaml", "transdecoder"),
    ("vipsania", "Vertebrata", "drusilla"),
    ("vipsania", " ETB1GO6Q ", "drusilla"),
])
def test_models_of_the_drusilla_flow(tool: str, model: str, method: str, tmp_path: Path) -> None:
    proc, published = run_pipeline(tmp_path, {**model_params(tmp_path, tool, model), **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    assert f"HC genes    : {method}" in proc.stdout
    drusilla_files = {f"intermediate/{tool}_lgb_filtered.gtf", "intermediate/drusilla_orfs.gtf",
                      "intermediate/hint_rescue.gtf"}
    if method == "drusilla":
        assert drusilla_files | {f"{tool}_evidence.gff3"} <= published, sorted(published)
        assert "intermediate/hc.gff3" not in published
    else:
        assert not drusilla_files & published
        assert {f"{tool}_evidence.gff3", "intermediate/hc.gff3"} <= published
    # citations.md: Drusilla and LightGBM for the flow, TD2 (the ORF finder of
    # hc_table) else; the hint rescue runs Tiberius, also for Vipsania; the
    # mammalian Tiberius models have no reference to the clade models
    text = (tmp_path / "out" / "citations.md").read_text()
    for name in ("Drusilla", "LightGBM"):
        assert (f"**{name}**" in text) == (method == "drusilla"), text
    assert ("**TD2**" in text) == (method == "transdecoder"), text
    assert "**TransDecoder**" not in text, text
    assert ("**Tiberius**" in text) == (tool == "tiberius" or method == "drusilla"), text
    assert ("multiple clades" in text) == (tool == "tiberius" and not model.startswith("mammalia")), text


def test_tiberius_result_with_drusilla(tmp_path: Path) -> None:
    """A reused prediction of a vertebrate model: the LightGBM filter and the rescue use it."""
    result = tmp_path / "previous.gtf"
    result.write_text("")
    params = model_params(tmp_path, "tiberius", "mammalia_softmasking_v2.yaml")
    params["tiberius"]["result"] = str(result)
    proc, published = run_pipeline(tmp_path, {**params, **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    assert "HC genes    : drusilla" in proc.stdout
    assert "RUN_TIBERIUS" not in proc.stdout and "SPLIT_GENOME" not in proc.stdout
    assert "GENEFINDER_LGB_FILTER" in proc.stdout and "HINT_RESCUE_TIBERIUS" in proc.stdout
    assert {"tiberius_evidence.gff3", "intermediate/tiberius_lgb_filtered.gtf",
            "intermediate/hint_rescue.gtf"} <= published, sorted(published)


def test_tiberius_result_without_model_is_not_drusilla(tmp_path: Path) -> None:
    result = tmp_path / "previous.gtf"
    result.write_text("")
    params = drusilla_params(tmp_path, "tiberius")
    params["tiberius"] = {"run": True, "result": str(result)}
    proc, published = run_pipeline(tmp_path, {**params, **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    assert "cannot tell which clade the prediction comes from" in proc.stdout + proc.stderr
    assert "HC genes    : transdecoder" in proc.stdout
    assert "intermediate/hc.gff3" in published


# ---- processes with their real scripts -------------------------------------------

PROCESS_TEST = """
include { %s } from '%s/modules/drusilla.nf'
workflow {
%s
}
"""

# Nextflow puts the bin/ directory next to the test script on the PATH.
FAKE_DRUSILLA = r'''#!/usr/bin/env python3
"""annotate: the ORFs are the input lines; orfs.partial.gtf records the input file,
its genes and --threads. Fails on a gene "fail". filter-subseq: copies."""
import os, re, shutil, sys
a = sys.argv[1:]
opt = lambda name: a[a.index(name) + 1]
if a[0] == "annotate":
    gtf, out = opt("--stringtie-gtf"), opt("--out-dir")
    os.makedirs(out, exist_ok=True)
    lines = [l for l in open(gtf) if l.strip() and not l.startswith("#")]
    genes = sorted({re.search(r'gene_id "([^"]+)"', l).group(1) for l in lines})
    if "fail" in genes:
        sys.exit("drusilla annotate failed on " + os.path.basename(gtf))
    open(os.path.join(out, "orfs.gtf"), "w").writelines(lines)
    with open(opt("--partial-out"), "w") as fh:
        fh.write(f"{os.path.basename(gtf)}\t{','.join(genes)}\t{opt('--threads')}\n")
elif a[0] == "filter-subseq":
    shutil.copy(opt("--orfs-gtf"), opt("--out-gtf"))
    open(opt("--report-tsv"), "w").close()
'''

FAKE_FIX_STOP = r'''#!/usr/bin/env python3
import shutil, sys
a = sys.argv[1:]
out = a[a.index("--out") + 1]
shutil.copy(a[a.index("--orfs") + 1], out)
open(out, "a").write("# fix_stop_by_miniprot.py " + " ".join(a) + "\n")
'''

# apply_lgb_model_gtf.py writes the transcripts that pass the threshold, with their class
FAKE_LGB = {
    "gff_to_cds_gtf.py": '#!/bin/sh\ncat "$1"\n',
    "compute_orf_features.py": r'''#!/usr/bin/env python3
import sys
open(sys.argv[sys.argv.index("--out") + 1], "w").write("transcript_id\n")
''',
    "apply_lgb_model_gtf.py": r'''#!/usr/bin/env python3
import sys
a = sys.argv[1:]
out = a[a.index("--out-gtf") + 1]
with open(out, "w") as fh:
    for tx, cls in [("t1", "correct"), ("t2", "partial"), ("t3", "correct")]:
        fh.write(f'c1\tTiberius\tCDS\t1\t300\t.\t+\t0\ttranscript_id "{tx}"; lgb_class "{cls}";\n')
open(out.replace(".gtf", ".scores.tsv"), "w").write(" ".join(a) + "\n")
''',
}


def run_process(tmp_path: Path, process: str, body: str, drusilla: dict,
                fakes: dict[str, str]) -> tuple[subprocess.CompletedProcess, Path]:
    """Runs one process of modules/drusilla.nf without -stub-run; returns the outdir."""
    (tmp_path / "bin").mkdir()
    for name, text in fakes.items():
        (tmp_path / "bin" / name).write_text(text)
        (tmp_path / "bin" / name).chmod(0o755)
    script = tmp_path / "test.nf"
    script.write_text(PROCESS_TEST % (process, ROOT, body))
    config = tmp_path / "test.config"
    config.write_text("process.executor = 'local'\nprocess.container = null\nprocess.cpus = 2\n")
    outdir = tmp_path / "out"
    params = tmp_path / "params.yaml"
    params.write_text(yaml.safe_dump({"outdir": str(outdir), "genome": str(DATA / "tiny.fa"), "drusilla": drusilla}))
    proc = subprocess.run(
        [NEXTFLOW, "run", str(script), "-c", str(config), "-params-file", str(params),
         "-work-dir", str(tmp_path / "work")],
        cwd=tmp_path, env=dict(os.environ, NXF_ANSI_LOG="false"), capture_output=True, text=True,
    )
    return proc, outdir


ANNOTATE = """
    DRUSILLA_ANNOTATE(file(params.gtf), file(params.genome), false, [])
    DRUSILLA_ANNOTATE.out.gtf.mix(DRUSILLA_ANNOTATE.out.partial, DRUSILLA_ANNOTATE.out.partial5)
        .subscribe { f -> f.copyTo("${params.outdir}/${f.name}") }
"""


def stringtie_lines(genes: list[str]) -> list[str]:
    """Two lines per gene, the genes interleaved (g1 g2 g1 g2 ...)."""
    return [f'c1\tStringTie\t{feature}\t1\t300\t1000\t+\t.\tgene_id "{g}"; transcript_id "{g}.1";\n'
            for feature in ("transcript", "exon") for g in genes]


def annotate(tmp_path: Path, shards: int, genes: list[str]) -> tuple[subprocess.CompletedProcess, Path]:
    gtf = tmp_path / "stringtie.filtered.gtf"
    gtf.write_text("# StringTie\n" + "".join(stringtie_lines(genes)))
    body = ANNOTATE.replace("params.gtf", f"'{gtf}'")
    return run_process(tmp_path, "DRUSILLA_ANNOTATE", body,
                       {"model": "vertebrates", "shards": shards, "min_coding_length": 200},
                       {"drusilla": FAKE_DRUSILLA})


@pytest.mark.parametrize("shards", [1, 2])
def test_drusilla_shards(shards: int, tmp_path: Path) -> None:
    """With shards > 1 the transcripts are split by gene, annotated in parallel and joined."""
    genes = ["g1", "g2", "g3", "g4", "g5"]
    proc, out = annotate(tmp_path, shards, genes)
    assert_ok(proc)
    assert sorted((out / "orfs.gtf").read_text().splitlines(keepends=True)) == sorted(stringtie_lines(genes))
    calls = sorted(line.split("\t") for line in (out / "orfs.partial.gtf").read_text().splitlines())
    if shards == 1:
        # the 2 CPUs of the task in one process
        assert calls == [["stringtie.filtered.gtf", "g1,g2,g3,g4,g5", "2"]]
    else:
        # genes in turn by their first line; the CPUs shared
        assert calls == [["shard_0.gtf", "g1,g3,g5", "1"], ["shard_1.gtf", "g2,g4", "1"]]
    assert (out / "orfs.partial5.gtf").read_text() == ""


def test_a_failing_shard_fails_the_task(tmp_path: Path) -> None:
    proc, _ = annotate(tmp_path, 2, ["g1", "fail", "g3"])
    assert proc.returncode != 0
    assert "drusilla annotate failed on shard_1.gtf" in proc.stdout + proc.stderr


FIX_ORFS = """
    FIX_ORFS(file("${projectDir}/orfs.gtf"), file("${projectDir}/orfs.partial.gtf"),
             file("${projectDir}/orfs.partial5.gtf"), file("${projectDir}/miniprot.gff"),
             file("${projectDir}/hc.gff"), file(params.genome))
"""

RAW_ORFS = 'c1\tDrusilla\tCDS\t1\t300\t.\t+\t0\ttranscript_id "o1"; gene_id "o1.g";\n'


@pytest.mark.parametrize("options,fix,fix_starts", [
    ({}, True, True),
    ({"fix_stop": False}, False, False),
    ({"fix_stop": "false", "fix_start": True}, False, False),   # the start fix needs the stop fix
    ({"fix_start": False}, True, False),
])
def test_fix_stop_and_fix_start(options: dict, fix: bool, fix_starts: bool, tmp_path: Path) -> None:
    for name in ("orfs.gtf", "orfs.partial.gtf", "orfs.partial5.gtf", "miniprot.gff", "hc.gff"):
        (tmp_path / name).write_text(RAW_ORFS if name == "orfs.gtf" else "")
    proc, out = run_process(tmp_path, "FIX_ORFS", FIX_ORFS, {"fix_stop": True, "fix_start": True, **options},
                            {"drusilla": FAKE_DRUSILLA, "fix_stop_by_miniprot.py": FAKE_FIX_STOP})
    assert_ok(proc)
    orfs = (out / "intermediate" / "drusilla_orfs.gtf").read_text()
    if not fix:
        assert orfs == RAW_ORFS
        return
    raw, call = orfs.splitlines()
    assert raw + "\n" == RAW_ORFS
    assert "--partial orfs.partial.gtf --partial5 orfs.partial5.gtf --miniprot miniprot.gff --hints hc.gff" in call
    assert ("--fix-starts --fix-starts-classes LORF_NOUPSTOP upLORF" in call) == fix_starts


LGB_FILTER = """
    GENEFINDER_LGB_FILTER('tiberius', file("${projectDir}/ab_initio.gff3"), file("${projectDir}/miniprot.gff"),
                          file("${projectDir}/hc.gff"), file(params.genome), file("${projectDir}/lgb.tar.gz"))
    GENEFINDER_LGB_FILTER.out.partial.subscribe { f -> f.copyTo("${params.outdir}/${f.name}") }
"""


def lgb_filter(tmp_path: Path, drusilla: dict) -> tuple[subprocess.CompletedProcess, Path]:
    for name in ("ab_initio.gff3", "miniprot.gff", "hc.gff", "lgb.tar.gz"):
        (tmp_path / name).write_text("")
    return run_process(tmp_path, "GENEFINDER_LGB_FILTER", LGB_FILTER, drusilla, FAKE_LGB)


@pytest.mark.parametrize("keep,kept", [
    (None, ["t1", "t3"]),                 # default: correct
    ("correct", ["t1", "t3"]),
    ("correct,partial", ["t1", "t2", "t3"]),
    ("partial correct", ["t1", "t2", "t3"]),
    ("wrong", []),                        # never in the scored GTF
])
def test_lgb_keep(keep: str | None, kept: list[str], tmp_path: Path) -> None:
    proc, out = lgb_filter(tmp_path, {"lgb_keep": keep, "lgb_threshold": 0.4, "lgb_model_sha256": "abc"})
    assert_ok(proc)
    tx = lambda path: [line.split('"')[1] for line in path.read_text().splitlines()]
    assert tx(out / "intermediate" / "tiberius_lgb_filtered.gtf") == kept
    # the candidates of the hint rescue do not depend on lgb_keep
    assert tx(out / "tiberius_lgb_partial.gtf") == ["t2"]
    call = (out / "intermediate" / "tiberius_lgb_scores.tsv").read_text()
    assert "--sha256 abc" in call and "--threshold 0.4" in call


def test_lgb_filter_without_checksum(tmp_path: Path) -> None:
    proc, out = lgb_filter(tmp_path, {"lgb_keep": "correct", "lgb_threshold": 0.5, "lgb_model_sha256": None})
    assert_ok(proc)
    assert "--sha256" not in (out / "intermediate" / "tiberius_lgb_scores.tsv").read_text()
