"""Nextflow stub runs: check the wiring and the published file names.

No tools, containers or GPU are needed, only a ``nextflow`` executable (or the
environment variable ``NEXTFLOW_BIN``). The tests are skipped without it.
"""
from __future__ import annotations

import hashlib
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


STATISTICS = {
    "qc/gene_set_statistics.txt", "qc/isoform_and_exon_structure.png", "qc/transcript_lengths.png",
    "qc/introns_per_gene.png",
}


def final_files(stem: str) -> set[str]:
    """The published annotation, its post-processing and QC files (default params)."""
    return {
        f"{stem}.gff3", f"{stem}.gtf", f"{stem}_proteins.fa", f"{stem}_cds.fa",
        f"intermediate/{stem}_sanity_filtered.gff3", "qc/sanity_filter.tsv",
        "qc/software_versions.tsv", "report.html", *STATISTICS,
    }


def expected_outputs(tool: str, mode: str) -> set[str]:
    if mode == "abinitio":
        return final_files(f"{tool}_ab_initio") | {f"intermediate/{tool}_ab_initio.gff3"}
    files = final_files(f"{tool}_evidence") | {
        f"intermediate/{tool}_ab_initio.gff3",
        f"intermediate/{tool}_merged.gff3",
        "intermediate/hc.gff3",
        "hintsfile.gff",
        "qc/gene_support.tsv",
        "qc/evidence_support.png",
    }
    if mode in ("rnaseq", "isoseq", "mixed"):
        files.add("qc/utr_report.tsv")
    return files


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
            "-with-trace", str(tmp_path / "trace.txt"),
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
        "citations.md",
        "methods.md",
        # post-processing (docs/postprocessing.md)
        "tiberius_evidence.gtf",
        "tiberius_evidence_cds.fa",
        "intermediate/tiberius_merged.gff3",
        "intermediate/tiberius_evidence_sanity_filtered.gff3",
        "qc/sanity_filter.tsv",
        "qc/utr_report.tsv",
        "qc/gene_support.tsv",
        "qc/gene_set_statistics.txt",
        "qc/isoform_and_exon_structure.png",
        "qc/transcript_lengths.png",
        "qc/introns_per_gene.png",
        "qc/evidence_support.png",
        "qc/software_versions.tsv",
        "report.html",
    }


def citations(tmp_path: Path) -> str:
    return (tmp_path / "out" / "citations.md").read_text()


def test_citations_follow_the_run(tmp_path: Path) -> None:
    """citations.md lists the references of the tools of this run only."""
    proc, _ = run_pipeline(tmp_path, {**GENEFINDER["tiberius"], **EVIDENCE["rnaseq"], "transdecoder": "td1"})
    assert_ok(proc)
    text = citations(tmp_path)
    for tool in ("Tiberius", "miniprot", "HISAT2", "StringTie", "TransDecoder", "DIAMOND"):
        assert f"**{tool}**" in text, text
    for tool in ("Vipsania", "minimap2", "TD2", "Drusilla", "LightGBM", "OrthoDB v12"):
        assert f"**{tool}**" not in text, text


def test_citations_ab_initio_vipsania(tmp_path: Path) -> None:
    proc, _ = run_pipeline(tmp_path, GENEFINDER["vipsania"])
    assert_ok(proc)
    text = citations(tmp_path)
    assert "**Vipsania**" in text and "**Nextflow**" in text
    assert "**Tiberius**" not in text and "**miniprot**" not in text


def test_citations_td2_and_isoseq(tmp_path: Path) -> None:
    """TD2 is the ORF finder of hc_table for clades without Drusilla."""
    proc, _ = run_pipeline(tmp_path, {**GENEFINDER["tiberius"], **EVIDENCE["isoseq"]})
    assert_ok(proc)
    text = citations(tmp_path)
    for tool in ("minimap2", "TD2", "PSAURON"):
        assert f"**{tool}**" in text, text
    assert "**TransDecoder**" not in text and "**HISAT2**" not in text


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


@pytest.mark.parametrize("tool", ["tiberius", "vipsania"])
def test_existing_gff3_result_without_model(tool: str, tmp_path: Path) -> None:
    """A GFF3 of the gene finder is enough: no model, no gene finder task, in every mode."""
    result = tmp_path / "previous.gff3"
    result.write_text("##gff-version 3\n")
    finder = {tool: {"run": True, "result": str(result)}}
    proc, published = run_pipeline(tmp_path, finder)
    assert_ok(proc)
    assert f"RUN_{tool.upper()}" not in proc.stdout and "SPLIT_GENOME" not in proc.stdout
    assert expected_outputs(tool, "abinitio") <= published
    proc, published = run_pipeline(tmp_path, {**finder, **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    assert f"RUN_{tool.upper()}" not in proc.stdout and "DOWNLOAD_" not in proc.stdout
    assert expected_outputs(tool, "rnaseq") <= published


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
    # SRA reads are published with keep_downloads only
    "sra_paired_list":   ({"rnaseq_sra_paired": ["SRR0000001"], "keep_downloads": True},
                          ["sra_downloads/rnaseq_sra_paired/SRR0000001_1.fastq.gz"]),
    "sra_single_string": ({"rnaseq_sra_single": "SRR0000002"}, []),
    "isoseq_sra_list":   ({"isoseq_sra": ["DRR0000003"], "keep_downloads": True},
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
    if not extra.get("keep_downloads"):
        assert not any(f.startswith("sra_downloads/") for f in published), sorted(published)


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


@pytest.mark.parametrize("flow", ["transdecoder", "drusilla"])
def test_file_names_with_whitespace(flow: str, tmp_path: Path) -> None:
    """A params file used without the launcher (which rejects such names)."""
    inputs = tmp_path / "my inputs"
    inputs.mkdir()
    names = {"tiny.fa": "tiny genome.fa", "tiny_proteins.faa": "tiny proteins.faa",
             "reads_1.fastq": "reads 1.fastq", "isoseq.fastq": "iso seq.fastq"}
    for src, dst in names.items():
        shutil.copy(DATA / src, inputs / dst)
    params = {
        "genome": str(inputs / "tiny genome.fa"),
        "proteins": str(inputs / "tiny proteins.faa"),
        "rnaseq_single": [str(inputs / "reads 1.fastq")],
        "isoseq": [str(inputs / "iso seq.fastq")],
    }
    params.update(drusilla_params(tmp_path, "tiberius") if flow == "drusilla" else GENEFINDER["tiberius"])
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert {"tiberius_evidence.gff3", "hintsfile.gff"} <= published


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
    assert "DOWNLOAD_TIBERIUS_WEIGHTS" not in proc.stdout


@pytest.mark.parametrize("model_dir", [False, True], ids=["download", "model-dir"])
def test_tiberius_weights_are_downloaded_once_without_model_dir(model_dir: bool, tmp_path: Path) -> None:
    cfg = tmp_path / "tiny.yaml"
    cfg.write_text('weights_url: "https://example.org/models/tiny_weights.tar.gz"\n')
    params = {"tiberius": {"run": True, "model_cfg": str(cfg)}}
    if model_dir:
        (tmp_path / "weights" / "tiny_weights").mkdir(parents=True)
        params["tiberius"]["model_dir"] = str(tmp_path / "weights")
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "tiberius_ab_initio.gff3" in published
    assert ("DOWNLOAD_TIBERIUS_WEIGHTS" in proc.stdout) != model_dir


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
        cfg.write_text('target_species: "Vertebrata"\n')
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
        f"intermediate/{tool}_ab_initio.gff3",
        f"intermediate/{tool}_lgb_filtered.gtf",
        f"intermediate/{tool}_lgb_scores.tsv",
        "intermediate/drusilla_orfs.gtf",
        "intermediate/hint_rescue.gtf",
        "citations.md",
        "methods.md",
    } | expected_outputs(tool, mode) - {"intermediate/hc.gff3"}
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


def hc_table(tmp_path: Path, clades: dict, orf_finder: str = "td2") -> str:
    """An HC table (conf/hc_genes.yaml format) for params.hc_table."""
    table = tmp_path / "hc_genes.yaml"
    table.write_text(yaml.safe_dump({"orf_finder": orf_finder, "clades": clades}))
    return str(table)


def test_drusilla_flow_needs_the_lgb_model(tmp_path: Path) -> None:
    """hc: drusilla without a LightGBM model in the table or the params: TransDecoder flow."""
    table = hc_table(tmp_path, {"Vertebrata": {"hc": "drusilla", "drusilla_model": "vertebrates"}})
    params = {**drusilla_params(tmp_path, "tiberius", lgb=False), "hc_table": table}
    proc, published = run_pipeline(tmp_path, {**params, **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    assert "no LightGBM model (params.drusilla.lgb_model)" in proc.stdout
    assert "HC genes    : transdecoder (clade Vertebrata, ORF finder td2)" in proc.stdout
    assert "intermediate/hc.gff3" in published


def test_hc_table_released_models_for_vertebrates(tmp_path: Path) -> None:
    """Without drusilla params, conf/hc_genes.yaml gives the Drusilla and LightGBM models."""
    params = {**drusilla_params(tmp_path, "tiberius"), "drusilla": {}}
    proc, published = run_pipeline(tmp_path, {**params, **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    assert "HC genes    : drusilla (clade Vertebrata, Drusilla model vertebrates)" in proc.stdout
    assert "intermediate/drusilla_orfs.gtf" in published


@pytest.mark.parametrize("tool,model,clade", [
    ("tiberius", "Insecta", "Insecta"),
    ("vipsania", "Fungi", "Fungi"),
    ("vipsania", "fh1kg88z", "fh1kg88z"),
], ids=["tiberius-insecta", "vipsania-fungi", "vipsania-model-id"])
def test_other_clades_get_td2(tool: str, model: str, clade: str, tmp_path: Path) -> None:
    if tool == "tiberius":
        cfg = tmp_path / "insecta.yaml"
        cfg.write_text(f"target_species: {model}\n")
        params = {"tiberius": {"run": True, "model_cfg": str(cfg)}}
    else:
        params = {"vipsania": {"run": True, "model": model}}
    proc, published = run_pipeline(tmp_path, {**params, **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    assert f"HC genes    : transdecoder (clade {clade}, ORF finder td2)" in proc.stdout
    assert "TD_ALL" in proc.stdout and "DRUSILLA" not in proc.stdout
    assert "intermediate/hc.gff3" in published


def test_transdecoder_param_overrides_the_table(tmp_path: Path) -> None:
    params = {**GENEFINDER["vipsania"], "transdecoder": "td1", **EVIDENCE["rnaseq"]}
    proc, _ = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "HC genes    : transdecoder (clade Fungi, ORF finder td1)" in proc.stdout


def test_transdecoder_param_is_checked(tmp_path: Path) -> None:
    proc, _ = run_pipeline(tmp_path, {**GENEFINDER["vipsania"], "transdecoder": "td3", **EVIDENCE["rnaseq"]})
    assert proc.returncode != 0
    assert "must be 'td1' or 'td2'" in proc.stdout + proc.stderr


def test_hc_table_adds_a_drusilla_clade(tmp_path: Path) -> None:
    """A clade added to the table with hc: drusilla gets the Drusilla flow and its models."""
    lgb = tmp_path / "insect_lgb.tar.gz"
    lgb.write_text("")
    table = hc_table(tmp_path, {
        "Insecta": {"hc": "drusilla", "drusilla_model": "insects", "lgb_model": str(lgb), "lgb_model_sha256": None},
        "Fungi": {"hc": "transdecoder", "orf_finder": "td1"},
    })
    cfg = tmp_path / "insecta.yaml"
    cfg.write_text("target_species: Insecta\n")
    params = {"tiberius": {"run": True, "model_cfg": str(cfg)}, "hc_table": table}
    proc, published = run_pipeline(tmp_path, {**params, **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    assert "HC genes    : drusilla (clade Insecta, Drusilla model insects)" in proc.stdout
    assert "intermediate/drusilla_orfs.gtf" in published
    # a clade's own orf_finder, and the vertebrates of conf/hc_genes.yaml are gone from this table
    (tmp_path / "fungi").mkdir()
    (tmp_path / "vert").mkdir()
    proc, _ = run_pipeline(tmp_path / "fungi", {**GENEFINDER["vipsania"], "hc_table": table, **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    assert "HC genes    : transdecoder (clade Fungi, ORF finder td1)" in proc.stdout
    proc, _ = run_pipeline(tmp_path / "vert", {**drusilla_params(tmp_path, "vipsania", lgb=False), "hc_table": table,
                                               **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    assert "HC genes    : transdecoder (clade etb1go6q, ORF finder td2)" in proc.stdout


def test_drusilla_forced_uses_the_drusilla_forced_clade(tmp_path: Path) -> None:
    """drusilla.run = true on a clade without Drusilla models: those of drusilla_forced (Vertebrata)."""
    params = {**GENEFINDER["vipsania"], "drusilla": {"run": True}, **EVIDENCE["rnaseq"]}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "HC genes    : drusilla (clade Fungi, Drusilla model vertebrates)" in proc.stdout
    assert "intermediate/drusilla_orfs.gtf" in published


def test_drusilla_forced_needs_models_for_the_clade(tmp_path: Path) -> None:
    table = hc_table(tmp_path, {"Vertebrata": {"hc": "drusilla", "drusilla_model": "vertebrates"}})
    params = {**GENEFINDER["vipsania"], "drusilla": {"run": True}, "hc_table": table, **EVIDENCE["rnaseq"]}
    proc, _ = run_pipeline(tmp_path, params)
    assert proc.returncode != 0
    assert "the HC table has none for clade Fungi" in proc.stdout + proc.stderr


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


def rescue_model(tmp_path: Path) -> str:
    """The model line that the HINT_RESCUE_TIBERIUS stub writes."""
    return (tmp_path / "out" / "intermediate" / "hint_rescue.gtf").read_text().strip()


@pytest.mark.parametrize("name,species,flow", [
    ("my_model.yaml", "Vertebrata", "drusilla"),
    ("mammals_v3.yaml", "mammalia", "drusilla"),
    ("vertebrates.yaml", "Insecta", "transdecoder"),
], ids=["vertebrata-any-name", "mammalia-any-case", "insecta-named-vertebrates"])
def test_drusilla_flow_reads_target_species(name: str, species: str, flow: str, tmp_path: Path) -> None:
    """The model's target_species decides, not its file name."""
    params = drusilla_params(tmp_path, "tiberius")
    cfg = tmp_path / "models" / name
    cfg.parent.mkdir()
    cfg.write_text(f"target_species: {species}\nncbi_tax_id: 1\n")
    params["tiberius"]["model_cfg"] = str(cfg)
    proc, published = run_pipeline(tmp_path, {**params, **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    assert f"HC genes    : {flow}" in proc.stdout
    if flow == "drusilla":
        # the rescue gets the model file itself, not a model of the image with the same stem
        assert rescue_model(tmp_path) == f"# model: {name}, staged weights: false"


@pytest.mark.parametrize("drusilla,download", [
    ({}, True),
    ({"cache_dir": "/no/such/cache"}, False),
    ({"weights": "/no/such/model.weights.h5", "config": "/no/such/arch.yaml"}, False),
    ({"model": "no_such_model"}, False),
], ids=["released-model", "cache-dir", "local-weights", "name-not-in-submodule"])
def test_drusilla_model_is_downloaded_once(drusilla: dict, download: bool, tmp_path: Path) -> None:
    """A model of the submodule's model_cfg/ is downloaded once by Nextflow, not by the shards."""
    params = drusilla_params(tmp_path, "tiberius")
    params["drusilla"].update(drusilla, shards=4)
    proc, _ = run_pipeline(tmp_path, {**params, **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    assert "HC genes    : drusilla" in proc.stdout
    assert ("DOWNLOAD_DRUSILLA_MODEL" in proc.stdout) == download


def test_drusilla_flow_needs_target_species(tmp_path: Path) -> None:
    params = drusilla_params(tmp_path, "tiberius")
    (tmp_path / "vertebrates.yaml").write_text("default_seq_len: 400050\n")
    proc, published = run_pipeline(tmp_path, {**params, **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    assert "has no target_species" in proc.stdout
    assert "HC genes    : transdecoder" in proc.stdout
    assert "intermediate/hc.gff3" in published


def test_drusilla_flow_forced_with_result_and_no_model(tmp_path: Path) -> None:
    """The rescue falls back to vertebrates instead of a model named 'null'."""
    result = tmp_path / "previous.gtf"
    result.write_text("")
    params = drusilla_params(tmp_path, "tiberius")
    params["tiberius"] = {"run": True, "result": str(result)}
    params["drusilla"]["run"] = True
    proc, published = run_pipeline(tmp_path, {**params, **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    assert "HC genes    : drusilla" in proc.stdout
    assert "the hint rescue uses the Tiberius model vertebrates" in proc.stdout
    # the weights of vertebrates (weights_url in the submodule's model_cfg/) are downloaded once
    assert "DOWNLOAD_RESCUE_WEIGHTS" in proc.stdout
    assert rescue_model(tmp_path) == "# model: vertebrates, staged weights: true"


@pytest.mark.parametrize("rescue_cfg,expected,weights", [
    ("mammalia_softmasking_v2", "mammalia_softmasking_v2", "true"),
    ("rescue.yaml", "rescue.yaml", "false"),
    ("no_such_model", "no_such_model", "false"),
], ids=["image-model-name", "model-file-without-weights-url", "name-not-in-submodule"])
def test_drusilla_rescue_model_cfg(rescue_cfg: str, expected: str, weights: str, tmp_path: Path) -> None:
    params = drusilla_params(tmp_path, "tiberius")
    if rescue_cfg.endswith(".yaml"):
        (tmp_path / rescue_cfg).write_text('target_species: "Vertebrata"\n')
        rescue_cfg = str(tmp_path / rescue_cfg)
    params["drusilla"]["rescue_model_cfg"] = rescue_cfg
    proc, _ = run_pipeline(tmp_path, {**params, **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    assert rescue_model(tmp_path) == f"# model: {expected}, staged weights: {weights}"
    assert ("DOWNLOAD_RESCUE_WEIGHTS" in proc.stdout) == (weights == "true")


def test_drusilla_rescue_downloads_the_weights_of_its_model_file(tmp_path: Path) -> None:
    params = drusilla_params(tmp_path, "tiberius")
    Path(params["tiberius"]["model_cfg"]).write_text(
        'target_species: "Vertebrata"\nweights_url: "https://example.org/models/tiny_weights.tar.gz"\n')
    proc, _ = run_pipeline(tmp_path, {**params, **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    # one download for the gene finder, one for the rescue
    assert "DOWNLOAD_TIBERIUS_WEIGHTS" in proc.stdout and "DOWNLOAD_RESCUE_WEIGHTS" in proc.stdout
    assert rescue_model(tmp_path) == "# model: vertebrates.yaml, staged weights: true"


def test_drusilla_rescue_stages_the_tiberius_weights(tmp_path: Path) -> None:
    model_dir = tmp_path / "weights"
    (model_dir / "tiny_weights").mkdir(parents=True)
    (model_dir / "tiny_weights" / "model.h5").write_text("")
    params = drusilla_params(tmp_path, "tiberius")
    params["tiberius"]["model_dir"] = str(model_dir)
    proc, _ = run_pipeline(tmp_path, {**params, **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    assert rescue_model(tmp_path) == "# model: vertebrates.yaml, staged weights: true"


def test_drusilla_rescue_of_vipsania_runs_uses_vertebrates(tmp_path: Path) -> None:
    proc, _ = run_pipeline(tmp_path, {**drusilla_params(tmp_path, "vipsania"), **EVIDENCE["rnaseq"]})
    assert_ok(proc)
    assert rescue_model(tmp_path) == "# model: vertebrates, staged weights: true"


# ---- pyVARUS output directories ------------------------------------------------

GENOME_MD5 = hashlib.md5((DATA / "tiny.fa").read_bytes()).hexdigest()


def varus_dir(tmp_path: Path, name: str, mode: str, gtf: bool = True, md5: str = GENOME_MD5) -> str:
    """A pyVARUS output directory: varus run (shortreads, longreads) or varus assemble --short --long (mixed)."""
    d = tmp_path / name
    d.mkdir()
    manifest = "VARUS.assembly.tsv" if mode == "mixed" else "VARUS.manifest.tsv"
    (d / manifest).write_text(
        f"#varus_manifest=1\n# a comment line\n#genome_md5={md5}\n#mode={mode}\n#assembly=stringtie.gtf\n"
        + ("" if mode == "mixed" else "#hints=hints.gff\n") + "batch\taccession\tn\n"
    )
    if gtf:
        (d / "stringtie.gtf").write_text(
            'seq1\tStringTie\ttranscript\t1\t60\t1000\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.1"; '
            'cov "5.0"; FPKM "1.0"; TPM "2.0";\n'
        )
    if mode != "mixed":
        (d / "hints.gff").write_text("seq1\tb2h\tintron\t11\t50\t3\t+\t.\tmult=3;pri=4;src=E\n")
    return str(d)


def flow_params(tmp_path: Path, flow: str) -> dict:
    """Tiberius with the TransDecoder or the Drusilla flow."""
    return drusilla_params(tmp_path, "tiberius") if flow == "drusilla" else GENEFINDER["tiberius"]


MAPPING_AND_ASSEMBLY = ("HISAT2", "MINIMAP2", "SAMTOOLS_MERGE", "BAM2HINTS", "STRINGTIE_ASSEMBLE")


@pytest.mark.parametrize("flow", ["transdecoder", "drusilla"])
@pytest.mark.parametrize("mode", ["rnaseq", "isoseq"])
def test_one_varus_directory_replaces_the_reads(mode: str, flow: str, tmp_path: Path) -> None:
    key, varus_mode = ("rnaseq_varus", "shortreads") if mode == "rnaseq" else ("isoseq_varus", "longreads")
    params = {**flow_params(tmp_path, flow), **PROTEINS, key: [varus_dir(tmp_path, "varus", varus_mode)]}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert f"Running mode: {mode}" in proc.stdout
    assert f"HC genes    : {flow}" in proc.stdout
    assert "VARUS_INPUT" in proc.stdout
    assert not [p for p in MAPPING_AND_ASSEMBLY if p in proc.stdout], proc.stdout
    # one hint source: used as it is
    assert "MERGE_INTRON_HINTS" not in proc.stdout
    assert {"tiberius_evidence.gff3", "hintsfile.gff"} <= published
    assert ("intermediate/drusilla_orfs.gtf" in published) == (flow == "drusilla")


def test_varus_directories_and_reads_are_combined(tmp_path: Path) -> None:
    """TransDecoder flow: the assemblies of all sources, the hints summed."""
    params = {**GENEFINDER["tiberius"], **PROTEINS, **SHORT_READS,
              "rnaseq_varus": [varus_dir(tmp_path, "a", "shortreads"), varus_dir(tmp_path, "b", "shortreads")]}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert proc.stdout.count("VARUS_INPUT_RNA") == 2
    assert "STRINGTIE_ASSEMBLE_RNA" in proc.stdout and "BAM2HINTS_RNA" in proc.stdout
    assert "MERGE_INTRON_HINTS_RNA" in proc.stdout
    assert "tiberius_evidence.gff3" in published


@pytest.mark.parametrize("tool", sorted(GENEFINDER))
def test_mixed_varus_is_the_drusilla_assembly(tool: str, tmp_path: Path) -> None:
    params = {
        **drusilla_params(tmp_path, tool), **PROTEINS,
        "rnaseq_varus": varus_dir(tmp_path, "short", "shortreads"),
        "isoseq_varus": [varus_dir(tmp_path, "long", "longreads")],
        "mixed_varus": varus_dir(tmp_path, "mixed", "mixed"),
    }
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "Running mode: mixed" in proc.stdout and "HC genes    : drusilla" in proc.stdout
    assert all(p in proc.stdout for p in ("VARUS_INPUT_RNA", "VARUS_INPUT_ISO", "VARUS_INPUT_MIX"))
    assert not [p for p in MAPPING_AND_ASSEMBLY if p in proc.stdout], proc.stdout
    assert {f"{tool}_evidence.gff3", "intermediate/drusilla_orfs.gtf", "hintsfile.gff"} <= published


def test_mixed_varus_is_not_used_by_transdecoder(tmp_path: Path) -> None:
    params = {
        **GENEFINDER["tiberius"], **PROTEINS,
        "rnaseq_varus": varus_dir(tmp_path, "short", "shortreads"),
        "isoseq_varus": varus_dir(tmp_path, "long", "longreads"),
        "mixed_varus": varus_dir(tmp_path, "mixed", "mixed"),
    }
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "mixed_varus is not used" in proc.stdout + proc.stderr
    assert "VARUS_INPUT_MIX" not in proc.stdout
    assert "intermediate/hc.gff3" in published


@pytest.mark.parametrize("case,message", [
    ("drusilla-two-dirs", "needs one StringTie assembly of all short reads"),
    ("drusilla-dir-and-reads", "needs one StringTie assembly of all Iso-Seq reads"),
    ("drusilla-mixed-without-mixed-varus", "varus assemble GENOME --short A/VARUS.bam --long B/VARUS.bam"),
    ("drusilla-mixed-varus-and-reads", "no other short-read or Iso-Seq input"),
    ("wrong-mode", "is a pyVARUS run in mode 'longreads'"),
    ("mixed-varus-wrong-mode", "is a pyVARUS run in mode 'shortreads'"),
    ("no-stringtie-gtf", "varus assemble GENOME --long"),
    ("not-a-directory", "is not a directory"),
])
def test_bad_varus_inputs_stop_the_run_early(case: str, message: str, tmp_path: Path) -> None:
    drusilla = {**drusilla_params(tmp_path, "tiberius"), **PROTEINS}
    transdecoder = {**GENEFINDER["tiberius"], **PROTEINS}
    short = lambda name="short": varus_dir(tmp_path, name, "shortreads")
    long = lambda name="long": varus_dir(tmp_path, name, "longreads")
    params = {
        "drusilla-two-dirs": lambda: {**drusilla, "rnaseq_varus": [short("a"), short("b")]},
        "drusilla-dir-and-reads": lambda: {**drusilla, **ISOSEQ, "isoseq_varus": long()},
        "drusilla-mixed-without-mixed-varus": lambda: {**drusilla, "rnaseq_varus": short(), "isoseq_varus": long()},
        "drusilla-mixed-varus-and-reads": lambda: {**drusilla, **SHORT_READS, "rnaseq_varus": short(),
                                                   "isoseq_varus": long(),
                                                   "mixed_varus": varus_dir(tmp_path, "m", "mixed")},
        "wrong-mode": lambda: {**transdecoder, "rnaseq_varus": long()},
        "mixed-varus-wrong-mode": lambda: {**drusilla, "rnaseq_varus": short(), "isoseq_varus": long(),
                                           "mixed_varus": short("m")},
        "no-stringtie-gtf": lambda: {**transdecoder, "isoseq_varus": varus_dir(tmp_path, "old", "longreads", gtf=False)},
        "not-a-directory": lambda: {**transdecoder, "rnaseq_varus": str(DATA / "tiny.fa")},
    }[case]()
    proc, _ = run_pipeline(tmp_path, params)
    assert proc.returncode != 0
    out = proc.stdout + proc.stderr
    assert message in out, out[-3000:]
    # checked when the workflow is built, before any task
    assert "Submitted process" not in out


# ---- VARUS_INPUT without -stub-run: the genome check -----------------------------

VARUS_INPUT_TEST = """
include { VARUS_INPUT } from '%s/modules/varus.nf'
include { varusInputs } from '%s/lib_nf/functions.nf'
workflow {
    VARUS_INPUT(channel.fromList(varusInputs(params.dir, 'rnaseq_varus', 'shortreads')), file(params.genome))
    VARUS_INPUT.out.gtf.view { f -> "GTF " + f.name }
    VARUS_INPUT.out.hints.view { f -> "HINTS " + f.text.trim() }
}
"""


def run_varus_input(tmp_path: Path, genome: Path, md5: str = GENOME_MD5) -> subprocess.CompletedProcess:
    script = tmp_path / "test.nf"
    script.write_text(VARUS_INPUT_TEST % (ROOT, ROOT))
    config = tmp_path / "test.config"
    config.write_text("process.executor = 'local'\nprocess.container = null\n")
    return subprocess.run(
        [NEXTFLOW, "run", str(script), "-c", str(config), "-work-dir", str(tmp_path / "work"),
         "--dir", varus_dir(tmp_path, "varus", "shortreads", md5=md5), "--genome", str(genome)],
        cwd=tmp_path, env=dict(os.environ, NXF_ANSI_LOG="false"), capture_output=True, text=True,
    )


@pytest.mark.parametrize("gzipped", [False, True], ids=["plain", "gzip"])
def test_varus_input_passes_the_files_on(gzipped: bool, tmp_path: Path) -> None:
    import gzip
    genome = DATA / "tiny.fa"
    if gzipped:
        genome = tmp_path / "tiny.fa.gz"
        genome.write_bytes(gzip.compress((DATA / "tiny.fa").read_bytes()))
    proc = run_varus_input(tmp_path, genome)
    assert_ok(proc)
    assert "GTF varus.stringtie.gtf" in proc.stdout
    assert "HINTS seq1\tb2h\tintron\t11\t50\t3\t+\t.\tmult=3;pri=4;src=E" in proc.stdout


def test_varus_input_of_another_genome_fails(tmp_path: Path) -> None:
    proc = run_varus_input(tmp_path, DATA / "tiny.fa", md5="0" * 32)
    assert proc.returncode != 0
    assert "the VARUS run used another genome" in proc.stdout + proc.stderr


# ---- StringTie assemblies made outside the pipeline (stringtie) ---------------

def stringtie_gtf(tmp_path: Path, name: str = "assembly.gtf") -> str:
    gtf = tmp_path / name
    gtf.write_text(
        'seq1\tStringTie\ttranscript\t1\t60\t1000\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.1"; '
        'cov "5.0"; FPKM "1.0"; TPM "2.0";\n'
    )
    return str(gtf)


@pytest.mark.parametrize("flow", ["transdecoder", "drusilla"])
@pytest.mark.parametrize("form", ["string", "list"])
def test_stringtie_assembly_replaces_the_reads(form: str, flow: str, tmp_path: Path) -> None:
    """Proteins and one StringTie GTF are a complete evidence run, in both flows."""
    gtf = stringtie_gtf(tmp_path)
    params = {**flow_params(tmp_path, flow), **PROTEINS, "stringtie": gtf if form == "string" else [gtf]}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "Running mode: rnaseq" in proc.stdout
    assert f"HC genes    : {flow}" in proc.stdout
    assert not [p for p in MAPPING_AND_ASSEMBLY if p in proc.stdout], proc.stdout
    assert {"tiberius_evidence.gff3", "tiberius_evidence_proteins.fa", "hintsfile.gff"} <= published
    assert ("intermediate/drusilla_orfs.gtf" in published) == (flow == "drusilla")
    assert ("intermediate/hc.gff3" in published) == (flow == "transdecoder")


def test_stringtie_assemblies_and_reads_are_combined(tmp_path: Path) -> None:
    """TransDecoder flow: several given assemblies (here a glob) are merged with those of the reads."""
    stringtie_gtf(tmp_path, "a.gtf")
    stringtie_gtf(tmp_path, "b.gtf")
    params = {**GENEFINDER["tiberius"], **PROTEINS, **ISOSEQ, "stringtie": str(tmp_path / "*.gtf")}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "Running mode: mixed" in proc.stdout
    assert "STRINGTIE_ASSEMBLE_ISO" in proc.stdout and "STRINGTIE_MERGE" in proc.stdout
    assert "STRINGTIE_ASSEMBLE_RNA" not in proc.stdout and "HISAT2" not in proc.stdout
    assert "tiberius_evidence.gff3" in published


def test_stringtie_assembly_with_an_existing_prediction(tmp_path: Path) -> None:
    """Nothing is predicted, mapped or assembled: an existing prediction and an existing assembly."""
    result = tmp_path / "previous.gff3"
    result.write_text("##gff-version 3\n")
    params = {"tiberius": {"run": True, "result": str(result)}, **PROTEINS, "stringtie": stringtie_gtf(tmp_path)}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert not [p for p in (*MAPPING_AND_ASSEMBLY, "RUN_TIBERIUS", "SPLIT_GENOME") if p in proc.stdout], proc.stdout
    assert "tiberius_evidence.gff3" in published


@pytest.mark.parametrize("case,message", [
    ("missing-file", "stringtie: no such file"),
    ("without-proteins", "need protein evidence"),
    ("drusilla-two-assemblies", "needs one StringTie assembly of all reads"),
    ("drusilla-assembly-and-reads", "needs one StringTie assembly of all reads"),
])
def test_bad_stringtie_inputs_stop_the_run_early(case: str, message: str, tmp_path: Path) -> None:
    drusilla = {**drusilla_params(tmp_path, "tiberius"), **PROTEINS}
    gtf = stringtie_gtf(tmp_path)
    params = {
        "missing-file": lambda: {**GENEFINDER["tiberius"], **PROTEINS, "stringtie": str(tmp_path / "no.gtf")},
        "without-proteins": lambda: {**GENEFINDER["tiberius"], "stringtie": gtf},
        "drusilla-two-assemblies": lambda: {**drusilla, "stringtie": [gtf, stringtie_gtf(tmp_path, "b.gtf")]},
        "drusilla-assembly-and-reads": lambda: {**drusilla, **SHORT_READS, "stringtie": gtf},
    }[case]()
    proc, _ = run_pipeline(tmp_path, params)
    assert proc.returncode != 0
    out = proc.stdout + proc.stderr
    assert message in out, out[-3000:]
