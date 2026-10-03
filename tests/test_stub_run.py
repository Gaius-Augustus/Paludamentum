"""Nextflow stub runs: check the wiring and the published file names.

No tools, containers or GPU are needed, only a ``nextflow`` executable (or the
environment variable ``NEXTFLOW_BIN``). The tests are skipped without it.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "tests" / "data"
NEXTFLOW = os.environ.get("NEXTFLOW_BIN") or shutil.which("nextflow")

pytestmark = [
    pytest.mark.nextflow,
    pytest.mark.skipif(NEXTFLOW is None, reason="nextflow executable not found"),
]

PROTEINS = {"proteins": str(DATA / "tiny_proteins.faa")}
SHORT_READS = {"rnaseq_paired": str(DATA / "reads_{1,2}.fastq")}
ISOSEQ = {"isoseq": [str(DATA / "isoseq.fastq")]}

EVIDENCE = {
    "abinitio": {},
    "proteins": PROTEINS,
    "rnaseq": {**PROTEINS, **SHORT_READS},
    "isoseq": {**PROTEINS, **ISOSEQ},
    "mixed": {**PROTEINS, **SHORT_READS, **ISOSEQ},
}

GENEFINDER = {
    # any existing file works as Tiberius model config in a stub run
    "tiberius": {"tiberius": {"run": True, "model_cfg": str(DATA / "tiny.fa")}},
    "vipsania": {"vipsania": {"run": True, "model": "Fungi"}},
}


def expected_outputs(tool: str, mode: str) -> set[str]:
    if mode == "abinitio":
        return {f"{tool}_ab_initio.gff3"}
    return {
        f"{tool}_evidence.gff3",
        f"{tool}_evidence_proteins.fa",
        f"intermediate/{tool}_ab_initio.gff3",
        "intermediate/hc.gff3",
        "hintsfile.gff",
    }


def run_pipeline(tmp_path: Path, params: dict) -> tuple[subprocess.CompletedProcess, set[str]]:
    outdir = tmp_path / "out"
    params = {"genome": str(DATA / "tiny.fa"), "outdir": str(outdir), **params}
    params_file = tmp_path / "params.yaml"
    params_file.write_text(yaml.safe_dump(params))
    proc = subprocess.run(
        [
            NEXTFLOW, "run", str(ROOT / "main.nf"), "-stub-run",
            "-params-file", str(params_file),
            "-c", str(ROOT / "tests" / "stub.config"),
            "-work-dir", str(tmp_path / "work"),
        ],
        cwd=tmp_path, env=dict(os.environ, NXF_ANSI_LOG="false"),
        capture_output=True, text=True,
    )
    published = set()
    if outdir.exists():
        published = {str(p.relative_to(outdir)) for p in outdir.rglob("*") if p.is_file()}
    return proc, published


def assert_ok(proc: subprocess.CompletedProcess) -> None:
    assert proc.returncode == 0, proc.stdout[-4000:] + proc.stderr[-2000:]


@pytest.mark.parametrize("mode", sorted(EVIDENCE))
@pytest.mark.parametrize("tool", sorted(GENEFINDER))
def test_published_files(tool: str, mode: str, tmp_path: Path) -> None:
    proc, published = run_pipeline(tmp_path, {**GENEFINDER[tool], **EVIDENCE[mode]})
    assert_ok(proc)
    missing = expected_outputs(tool, mode) - published
    assert not missing, f"missing {missing}; published {sorted(published)}"
    other = "vipsania" if tool == "tiberius" else "tiberius"
    assert not [f for f in published if Path(f).name.startswith(other)], sorted(published)
    if mode == "abinitio":
        assert f"{tool}_evidence.gff3" not in published


def test_tiberius_file_names_are_unchanged(tmp_path: Path) -> None:
    """The names that Tiberius users rely on (before the gene finder abstraction)."""
    proc, published = run_pipeline(tmp_path, {**GENEFINDER["tiberius"], **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    assert published == {
        "tiberius_evidence.gff3",
        "tiberius_evidence_proteins.fa",
        "intermediate/tiberius_ab_initio.gff3",
        "intermediate/hc.gff3",
        "hintsfile.gff",
    }


def test_legacy_mode_name_tiberius(tmp_path: Path) -> None:
    proc, published = run_pipeline(tmp_path, {**GENEFINDER["tiberius"], "mode": "tiberius"})
    assert_ok(proc)
    assert "Running mode: abinitio" in proc.stdout
    assert "tiberius_ab_initio.gff3" in published


def test_vipsania_is_chunked_without_finetuning(tmp_path: Path) -> None:
    proc, published = run_pipeline(tmp_path, GENEFINDER["vipsania"])
    assert_ok(proc)
    assert "Gene finder : vipsania" in proc.stdout
    logs = {f for f in published if f.startswith("intermediate/vipsania/")}
    # the stub of SPLIT_GENOME makes two chunks
    assert logs == {
        "intermediate/vipsania/vipsania.genome_chunk_1.fa.log",
        "intermediate/vipsania/vipsania.genome_chunk_2.fa.log",
    }


def test_vipsania_finetuning_uses_the_whole_genome(tmp_path: Path) -> None:
    params = {"vipsania": {"run": True, "model": "Fungi", "finetune": True, "finetune_epochs": 1}}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "SPLIT_GENOME" not in proc.stdout
    assert "intermediate/vipsania/vipsania.tiny.fa.log" in published
    assert "vipsania_ab_initio.gff3" in published


def test_vipsania_local_model_dir_skips_the_download(tmp_path: Path) -> None:
    model_dir = tmp_path / "models"
    (model_dir / "abc123").mkdir(parents=True)
    (model_dir / "versions.json").write_text("{}")
    params = {"vipsania": {"run": True, "model": "abc123", "model_dir": str(model_dir)}}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "DOWNLOAD_VIPSANIA_MODEL" not in proc.stdout
    assert "vipsania_ab_initio.gff3" in published


def test_existing_result_is_reused(tmp_path: Path) -> None:
    result = tmp_path / "previous.gtf"
    result.write_text("")
    params = {"vipsania": {"run": True, "result": str(result)}, **PROTEINS}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "RUN_VIPSANIA" not in proc.stdout and "SPLIT_GENOME" not in proc.stdout
    assert "vipsania_evidence.gff3" in published


def test_evidence_without_gene_finder(tmp_path: Path) -> None:
    proc, published = run_pipeline(tmp_path, PROTEINS)
    assert_ok(proc)
    assert "Gene finder : none" in proc.stdout
    assert {"hintsfile.gff", "intermediate/hc.gff3"} <= published
    assert not [f for f in published if "evidence.gff3" in f]


def test_explicit_genefinder_wins(tmp_path: Path) -> None:
    params = {**GENEFINDER["tiberius"], **GENEFINDER["vipsania"], "genefinder": "vipsania"}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "vipsania_ab_initio.gff3" in published


def test_two_gene_finders_are_rejected(tmp_path: Path) -> None:
    proc, _ = run_pipeline(tmp_path, {**GENEFINDER["tiberius"], **GENEFINDER["vipsania"]})
    assert proc.returncode != 0
    assert "More than one gene finder" in proc.stdout + proc.stderr


def test_vipsania_needs_a_model(tmp_path: Path) -> None:
    proc, _ = run_pipeline(tmp_path, {"vipsania": {"run": True}})
    assert proc.returncode != 0
    assert "params.vipsania.model is required" in proc.stdout + proc.stderr


def test_cli_end_to_end_stub_run(tmp_path: Path) -> None:
    """paludamentum builds the params file and starts Nextflow with the config."""
    proc = subprocess.run(
        [
            sys.executable, "-m", "paludamentum",
            "--nf_config", str(ROOT / "tests" / "stub.config"),
            "--genome", str(DATA / "tiny.fa"), "--model_cfg", str(DATA / "tiny.fa"),
            "--proteins", str(DATA / "tiny_proteins.faa"),
            "--outdir", "out", "--work_dir", str(tmp_path / "work"),
            "--nextflow_bin", NEXTFLOW, "--skip_singularity_check",
            "--", "-stub-run",
        ],
        cwd=tmp_path, env=dict(os.environ, NXF_ANSI_LOG="false", PYTHONPATH=str(ROOT)),
        capture_output=True, text=True,
    )
    assert_ok(proc)
    assert "Gene finder: tiberius" in proc.stdout
    outdir = tmp_path / "out"
    written = yaml.safe_load((outdir / "params.yaml").read_text())
    assert written["tiberius"]["run"] is True and written["outdir"] == str(outdir)
    published = {str(p.relative_to(outdir)) for p in outdir.rglob("*") if p.is_file()}
    assert {"tiberius_evidence.gff3", "tiberius_evidence_proteins.fa", "params.yaml"} <= published


R1, R2 = str(DATA / "reads_1.fastq"), str(DATA / "reads_2.fastq")

INPUT_FORMS = {
    # any file stands in for a BAM in a stub run
    "bam_list":          ({"rnaseq_bam": [R1], "mode": "rnaseq"}, []),
    "bam_string":        ({"rnaseq_bam": R1, "mode": "rnaseq"}, []),
    "paired_pairs_list": ({"rnaseq_paired": [[R1, R2]]}, []),
    "paired_flat_two":   ({"rnaseq_paired": [R1, R2]}, []),
    "single_list":       ({"rnaseq_single": [R1]}, []),
    "sra_paired_list":   ({"rnaseq_sra_paired": ["SRR0000001"]},
                          ["sra_downloads/rnaseq_sra_paired/SRR0000001_1.fastq.gz"]),
    "sra_single_string": ({"rnaseq_sra_single": "SRR0000002"},
                          ["sra_downloads/rnaseq_sra_single/SRR0000002.fastq.gz"]),
    "isoseq_sra_list":   ({"isoseq_sra": ["DRR0000003"]},
                          ["sra_downloads/isoseq_sra/DRR0000003.fastq.gz"]),
    "two_protein_files": ({"proteins": [str(DATA / "tiny_proteins.faa"), str(DATA / "tiny_proteins.faa")],
                           "rnaseq_paired": str(DATA / "reads_{1,2}.fastq")}, []),
}


@pytest.mark.parametrize("form", sorted(INPUT_FORMS))
def test_input_forms(form: str, tmp_path: Path) -> None:
    """Every accepted form of the evidence inputs is wired through to the final annotation."""
    extra, downloads = INPUT_FORMS[form]
    params = {**GENEFINDER["tiberius"], **PROTEINS, **extra}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert {"tiberius_evidence.gff3", "intermediate/hc.gff3", "hintsfile.gff"} <= published, sorted(published)
    for f in downloads:
        assert f in published, (f, sorted(published))


@pytest.mark.parametrize("bam", [[R1], R1], ids=["list", "string"])
def test_bam_alone_is_inferred_as_rnaseq(bam, tmp_path: Path) -> None:
    proc, _ = run_pipeline(tmp_path, {**GENEFINDER["tiberius"], **PROTEINS, "rnaseq_bam": bam})
    assert_ok(proc)
    assert "Running mode: rnaseq" in proc.stdout


# ---------------------------------------------------------------- input handling and guards

def test_gzipped_genome_and_proteins(tmp_path: Path) -> None:
    """Both are decompressed under their own names, so they do not collide in MINIPROT_ALIGN."""
    import gzip
    genome_gz = tmp_path / "genome.fa.gz"
    proteins_gz = tmp_path / "proteins.faa.gz"
    genome_gz.write_bytes(gzip.compress((DATA / "tiny.fa").read_bytes()))
    proteins_gz.write_bytes(gzip.compress((DATA / "tiny_proteins.faa").read_bytes()))
    params = {**GENEFINDER["tiberius"], "genome": str(genome_gz), "proteins": str(proteins_gz)}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "tiberius_evidence.gff3" in published


def test_two_bams_with_the_same_file_name(tmp_path: Path) -> None:
    """Typical STAR output: <sample>/Aligned.out.bam for every sample."""
    bams = []
    for sample in ("s1", "s2"):
        (tmp_path / sample).mkdir()
        bam = tmp_path / sample / "Aligned.out.bam"
        bam.write_text("")
        bams.append(str(bam))
    proc, published = run_pipeline(tmp_path, {**GENEFINDER["tiberius"], **PROTEINS, "rnaseq_bam": bams})
    assert_ok(proc)
    assert {"tiberius_evidence.gff3", "hintsfile.gff"} <= published


def test_tiberius_model_dir_is_staged(tmp_path: Path) -> None:
    model_dir = tmp_path / "weights"
    (model_dir / "tiny_weights").mkdir(parents=True)
    (model_dir / "tiny_weights" / "model.h5").write_text("")
    params = {"tiberius": {"run": True, "model_cfg": str(DATA / "tiny.fa"), "model_dir": str(model_dir)}}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "tiberius_ab_initio.gff3" in published


@pytest.mark.parametrize("params,message", [
    ({"mode": "rnaseqq", **PROTEINS}, "Unknown params.mode"),
    ({"mode": "rnaseq", **PROTEINS}, "needs short reads"),
    ({"mode": "isoseq", **PROTEINS}, "needs Iso-Seq reads"),
    ({"mode": "proteins"}, "needs protein evidence"),
    (ISOSEQ, "need protein evidence"),
    ({"tiberius": {"run": True, "result": "/nonexistent/old.gtf"}, **PROTEINS}, "does not exist"),
    ({"tiberius": {"run": True, "model_cfg": "/nonexistent/m.yaml"}}, "not a file"),
], ids=["mode-typo", "rnaseq-without-reads", "isoseq-without-reads", "proteins-without-proteins",
        "isoseq-without-proteins", "missing-result", "missing-model-cfg"])
def test_bad_inputs_stop_the_run_early(params, message, tmp_path: Path) -> None:
    """Instead of a run that ends with exit 0 and without its main outputs, or a GPU task for nothing."""
    base = GENEFINDER["tiberius"] if "tiberius" not in params else {}
    proc, _ = run_pipeline(tmp_path, {**base, **params})
    assert proc.returncode != 0
    assert message in proc.stdout + proc.stderr, proc.stdout[-2000:] + proc.stderr[-1000:]


# ---- Drusilla flow -----------------------------------------------------------

def drusilla_params(tmp_path: Path, tool: str, lgb: bool = True) -> dict:
    """A vertebrate gene finder model, and a LightGBM model file if lgb."""
    if tool == "tiberius":
        cfg = tmp_path / "vertebrates.yaml"
        cfg.write_text("")
        params = {"tiberius": {"run": True, "model_cfg": str(cfg)}}
    else:
        params = {"vipsania": {"run": True, "model": "etb1go6q"}}
    if lgb:
        model = tmp_path / "lgb.tar.gz"
        model.write_text("")
        params["drusilla"] = {"lgb_model": str(model)}
    else:
        params["drusilla"] = {"lgb_model": None}
    return params


@pytest.mark.parametrize("mode", ["rnaseq", "isoseq", "mixed"])
@pytest.mark.parametrize("tool", sorted(GENEFINDER))
def test_drusilla_flow_for_vertebrate_models(tool: str, mode: str, tmp_path: Path) -> None:
    proc, published = run_pipeline(tmp_path, {**drusilla_params(tmp_path, tool), **EVIDENCE[mode]})
    assert_ok(proc)
    assert "HC genes    : drusilla" in proc.stdout
    assert "TD_ALL" not in proc.stdout
    published = {f for f in published if not f.startswith("intermediate/vipsania/")}
    assert published == {
        f"{tool}_evidence.gff3",
        f"{tool}_evidence_proteins.fa",
        f"intermediate/{tool}_ab_initio.gff3",
        f"intermediate/{tool}_lgb_filtered.gtf",
        f"intermediate/{tool}_lgb_scores.tsv",
        "intermediate/drusilla_orfs.gtf",
        "intermediate/hint_rescue.gtf",
        "hintsfile.gff",
    }
    # one StringTie assembly of all reads
    assert proc.stdout.count("STRINGTIE_ASSEMBLE") == 1, proc.stdout


def test_drusilla_flow_not_used_without_transcripts(tmp_path: Path) -> None:
    proc, published = run_pipeline(tmp_path, {**drusilla_params(tmp_path, "tiberius"), **EVIDENCE["proteins"]})
    assert_ok(proc)
    assert "DRUSILLA" not in proc.stdout
    assert "intermediate/hc.gff3" in published


def test_drusilla_flow_not_used_for_other_models(tmp_path: Path) -> None:
    params = {**GENEFINDER["vipsania"], "drusilla": {"lgb_model": str(DATA / "tiny.fa")}, **EVIDENCE["rnaseq"]}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "HC genes    : transdecoder" in proc.stdout
    assert "intermediate/hc.gff3" in published


def test_drusilla_flow_needs_the_lgb_model(tmp_path: Path) -> None:
    proc, published = run_pipeline(tmp_path, {**drusilla_params(tmp_path, "tiberius", lgb=False), **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    assert "params.drusilla.lgb_model is not set" in proc.stdout
    assert "intermediate/hc.gff3" in published


def test_drusilla_flow_can_be_switched_off(tmp_path: Path) -> None:
    params = drusilla_params(tmp_path, "tiberius")
    params["drusilla"]["run"] = False
    proc, published = run_pipeline(tmp_path, {**params, **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    assert "HC genes    : transdecoder" in proc.stdout
    assert "intermediate/drusilla_orfs.gtf" not in published


def test_drusilla_flow_forced_needs_transcripts(tmp_path: Path) -> None:
    params = drusilla_params(tmp_path, "tiberius")
    params["drusilla"]["run"] = True
    proc, _ = run_pipeline(tmp_path, {**params, **EVIDENCE["proteins"]})
    assert proc.returncode != 0
    assert "needs transcripts" in proc.stdout + proc.stderr


def test_drusilla_flow_without_rescue(tmp_path: Path) -> None:
    params = drusilla_params(tmp_path, "tiberius")
    params["drusilla"]["rescue"] = False
    proc, published = run_pipeline(tmp_path, {**params, **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    assert "HINT_RESCUE" not in proc.stdout
    assert "intermediate/drusilla_orfs.gtf" in published
    assert "intermediate/hint_rescue.gtf" not in published
