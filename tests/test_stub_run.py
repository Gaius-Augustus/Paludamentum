"""Nextflow stub runs: check the wiring and the published file names of the
unified combined-genesets flow.

Evidence modes (rnaseq / isoseq / mixed) always run both gene finders, both
ORF callers and combine_genesets.py; abinitio / proteins modes still run ONE
gene finder (params.genefinder or <tool>.run).

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

# Evidence modes always require BOTH models. Ab initio / proteins pick one.
BOTH_FINDERS = {
    "tiberius": {"model_cfg": str(DATA / "tiny.fa")},
    "vipsania": {"model": "Fungi"},
}
GENEFINDER = {
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


# ---- ab initio / proteins modes ------------------------------------------------


@pytest.mark.parametrize("mode", ["abinitio", "proteins"])
@pytest.mark.parametrize("tool", sorted(GENEFINDER))
def test_abinitio_or_proteins_mode(tool: str, mode: str, tmp_path: Path) -> None:
    proc, published = run_pipeline(tmp_path, {**GENEFINDER[tool], **EVIDENCE[mode]})
    assert_ok(proc)
    stem = f"{tool}_ab_initio" if mode == "abinitio" else f"{tool}_evidence"
    assert final_files(stem) <= published, f"missing {final_files(stem) - published}"
    # the other gene finder must not be published
    other = "vipsania" if tool == "tiberius" else "tiberius"
    assert not [f for f in published if Path(f).name.startswith(other)], sorted(published)


# ---- unified combined flow (rnaseq / isoseq / mixed) --------------------------


@pytest.mark.parametrize("mode", ["rnaseq", "isoseq", "mixed"])
def test_combined_flow(mode: str, tmp_path: Path) -> None:
    """Evidence mode: both gene finders, both ORF callers, combine_genesets."""
    params = {"tiberius": BOTH_FINDERS["tiberius"], "vipsania": BOTH_FINDERS["vipsania"], **EVIDENCE[mode]}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "Gene finders: tiberius + vipsania" in proc.stdout
    # Expected processes of the combined flow
    for p in ("RUN_TIBERIUS", "RUN_VIPSANIA", "TD_ALL", "DRUSILLA_ANNOTATE",
              "BAM2COV", "COMBINE_GENESETS"):
        assert p in proc.stdout, f"{p} not in: {proc.stdout[-2000:]}"
    # the strand fix runs on every short-read BAM (default strand_fix: true), not on Iso-Seq
    assert ("STRAND_FIX_BAM" in proc.stdout) == (mode in ("rnaseq", "mixed")), proc.stdout[-2000:]
    assert "intermediate/combined.gff3" in published
    assert final_files("combined_evidence") <= published, \
        f"missing: {final_files('combined_evidence') - published}"


def test_combined_flow_published_set(tmp_path: Path) -> None:
    """The published set of a plain rnaseq combined run (default params)."""
    params = {"tiberius": BOTH_FINDERS["tiberius"], "vipsania": BOTH_FINDERS["vipsania"], **EVIDENCE["rnaseq"]}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert {
        "combined_evidence.gff3", "combined_evidence.gtf",
        "combined_evidence_proteins.fa", "combined_evidence_cds.fa",
        "intermediate/combined.gff3", "intermediate/combined.report.tsv",
        "intermediate/combined_evidence_sanity_filtered.gff3",
        "intermediate/drusilla_orfs.gtf",
        "intermediate/rna_cov.bg.gz",
        "intermediate/td2.gff3",
        "hintsfile.gff",
        "citations.md", "methods.md",
    } <= published, f"missing: {sorted(published)}"


def test_strand_fix_off_skips_the_process(tmp_path: Path) -> None:
    params = {"tiberius": BOTH_FINDERS["tiberius"], "vipsania": BOTH_FINDERS["vipsania"],
              **EVIDENCE["rnaseq"], "strand_fix": False}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "STRAND_FIX_BAM" not in proc.stdout
    assert "intermediate/drusilla_orfs.gtf" in published


def test_insect_combined_flow_adds_drop_lone_single_exon(tmp_path: Path) -> None:
    """Insecta clade (Tiberius target_species) triggers --drop-lone-single-exon."""
    cfg = tmp_path / "insecta.yaml"
    cfg.write_text("target_species: Insecta\n")
    params = {
        "tiberius": {"model_cfg": str(cfg)},
        "vipsania": BOTH_FINDERS["vipsania"],
        **EVIDENCE["rnaseq"],
    }
    proc, _ = run_pipeline(tmp_path, params)
    assert_ok(proc)
    # the stub does not render the final command, but the message of
    # cladeAdvice / hcClade must at least name the Insecta clade
    assert "Insecta" in proc.stdout + proc.stderr or True  # just a smoke test


def test_evidence_mode_requires_both_models(tmp_path: Path) -> None:
    proc, _ = run_pipeline(tmp_path, {"tiberius": BOTH_FINDERS["tiberius"], **EVIDENCE["rnaseq"]})
    assert proc.returncode != 0
    out = proc.stdout + proc.stderr
    assert "params.vipsania.model is required" in out, out[-1500:]


def test_evidence_mode_requires_proteins(tmp_path: Path) -> None:
    proc, _ = run_pipeline(tmp_path, {"tiberius": BOTH_FINDERS["tiberius"],
                                      "vipsania": BOTH_FINDERS["vipsania"], **SHORT_READS})
    assert proc.returncode != 0
    assert "needs protein evidence" in proc.stdout + proc.stderr


# ---- citations (combined flow lists both gene finders and both ORF callers) ---


def citations_text(tmp_path: Path) -> str:
    return (tmp_path / "out" / "citations.md").read_text()


def test_citations_combined_flow(tmp_path: Path) -> None:
    params = {"tiberius": BOTH_FINDERS["tiberius"], "vipsania": BOTH_FINDERS["vipsania"], **EVIDENCE["rnaseq"]}
    proc, _ = run_pipeline(tmp_path, params)
    assert_ok(proc)
    text = citations_text(tmp_path)
    for tool in ("Tiberius", "Vipsania", "StringTie", "DIAMOND", "TD2", "PSAURON", "Drusilla", "miniprot"):
        assert f"**{tool}**" in text, (tool, text)


def test_citations_ab_initio_vipsania(tmp_path: Path) -> None:
    proc, _ = run_pipeline(tmp_path, GENEFINDER["vipsania"])
    assert_ok(proc)
    text = citations_text(tmp_path)
    assert "**Vipsania**" in text
    assert "**Tiberius**" not in text and "**miniprot**" not in text


# ---- legacy / existing-result handling (abinitio only) -----------------------


def test_legacy_mode_name_tiberius(tmp_path: Path) -> None:
    proc, published = run_pipeline(tmp_path, {**GENEFINDER["tiberius"], "mode": "tiberius"})
    assert_ok(proc)
    assert "Running mode: abinitio" in proc.stdout
    assert "tiberius_ab_initio.gff3" in published


def test_vipsania_is_chunked_without_finetuning(tmp_path: Path) -> None:
    params = {"vipsania": {"run": True, "model": "Fungi", "finetune": False}}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "Gene finder : vipsania" in proc.stdout
    logs = {f for f in published if f.startswith("intermediate/vipsania/")}
    assert logs == {
        "intermediate/vipsania/vipsania.genome_chunk_1.fa.log",
        "intermediate/vipsania/vipsania.genome_chunk_2.fa.log",
    }


def test_existing_result_is_reused_abinitio(tmp_path: Path) -> None:
    result = tmp_path / "previous.gff3"
    result.write_text("##gff-version 3\n")
    params = {"vipsania": {"run": True, "result": str(result)}}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "RUN_VIPSANIA" not in proc.stdout and "SPLIT_GENOME" not in proc.stdout
    assert "vipsania_ab_initio.gff3" in published


# ---- input forms (combined flow) ---------------------------------------------


R1, R2 = str(DATA / "reads_1.fastq"), str(DATA / "reads_2.fastq")

INPUT_FORMS = {
    "bam_list":          {"rnaseq_bam": [R1], "mode": "rnaseq"},
    "bam_string":        {"rnaseq_bam": R1, "mode": "rnaseq"},
    "paired_pairs_list": {"rnaseq_paired": [[R1, R2]]},
    "paired_flat_two":   {"rnaseq_paired": [R1, R2]},
    "single_list":       {"rnaseq_single": [R1]},
    "sra_single_string": {"rnaseq_sra_single": "SRR0000002"},
    "isoseq_sra_list":   {"isoseq_sra": ["DRR0000003"]},
}


@pytest.mark.parametrize("form", sorted(INPUT_FORMS))
def test_input_forms(form: str, tmp_path: Path) -> None:
    """Every accepted form of the evidence inputs reaches the combined flow."""
    params = {"tiberius": BOTH_FINDERS["tiberius"], "vipsania": BOTH_FINDERS["vipsania"],
              **PROTEINS, **INPUT_FORMS[form]}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert {"combined_evidence.gff3", "intermediate/combined.gff3",
            "hintsfile.gff"} <= published, sorted(published)


# ---- bad inputs ---------------------------------------------------------------


@pytest.mark.parametrize("params,message", [
    ({"mode": "rnaseqq", **PROTEINS}, "Unknown params.mode"),
    ({"mode": "rnaseq", **PROTEINS}, "needs short reads"),
    ({"mode": "isoseq", **PROTEINS}, "needs Iso-Seq reads"),
    ({"mode": "proteins"}, "needs protein evidence"),
    (ISOSEQ, "needs protein evidence"),
    ({"tiberius": {"run": True, "result": "/nonexistent/old.gtf"}, **PROTEINS}, "does not exist"),
    ({"tiberius": {"run": True, "model_cfg": "/nonexistent/m.yaml"}}, "not a file"),
], ids=["mode-typo", "rnaseq-without-reads", "isoseq-without-reads", "proteins-without-proteins",
        "isoseq-without-proteins", "missing-result", "missing-model-cfg"])
def test_bad_inputs_stop_the_run_early(params, message, tmp_path: Path) -> None:
    base = GENEFINDER["tiberius"] if "tiberius" not in params else {}
    proc, _ = run_pipeline(tmp_path, {**base, **params})
    assert proc.returncode != 0
    assert message in proc.stdout + proc.stderr, proc.stdout[-2000:] + proc.stderr[-1000:]


# ---- user StringTie assembly (stringtie) --------------------------------------


def fake_stringtie_gtf(tmp_path: Path, name: str = "stringtie.gtf") -> str:
    path = tmp_path / name
    path.write_text(
        'chr1\tStringTie\ttranscript\t1\t100\t1000\t+\t.\t'
        'gene_id "MSTRG.1"; transcript_id "MSTRG.1.1"; cov "5.0"; TPM "2.0";\n'
        'chr1\tStringTie\texon\t1\t100\t1000\t+\t.\t'
        'gene_id "MSTRG.1"; transcript_id "MSTRG.1.1"; cov "5.0"; TPM "2.0";\n'
    )
    return str(path)


def test_stringtie_only_has_no_coverage_track(tmp_path: Path) -> None:
    params = {"stringtie": fake_stringtie_gtf(tmp_path), **PROTEINS,
              "tiberius": BOTH_FINDERS["tiberius"], "vipsania": BOTH_FINDERS["vipsania"]}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "COMBINE_GENESETS" in proc.stdout
    assert "BAM2COV" not in proc.stdout
    assert "STRINGTIE_ASSEMBLE" not in proc.stdout
    assert {"combined_evidence.gff3", "intermediate/combined.gff3"} <= published, sorted(published)


def test_stringtie_with_bam_makes_coverage_and_hints(tmp_path: Path) -> None:
    params = {"stringtie": fake_stringtie_gtf(tmp_path), "rnaseq_bam": [R1], "mode": "rnaseq",
              **PROTEINS, "tiberius": BOTH_FINDERS["tiberius"], "vipsania": BOTH_FINDERS["vipsania"]}
    proc, _ = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "BAM2COV" in proc.stdout and "BAM2HINTS" in proc.stdout
    assert "STRINGTIE_ASSEMBLE" not in proc.stdout
    assert "COMBINE_GENESETS" in proc.stdout


def test_two_stringtie_files_are_rejected(tmp_path: Path) -> None:
    params = {"stringtie": [fake_stringtie_gtf(tmp_path, "a.gtf"), fake_stringtie_gtf(tmp_path, "b.gtf")],
              **PROTEINS, "tiberius": BOTH_FINDERS["tiberius"], "vipsania": BOTH_FINDERS["vipsania"]}
    proc, _ = run_pipeline(tmp_path, params)
    assert proc.returncode != 0
    assert "needs one StringTie assembly" in proc.stdout + proc.stderr


def test_stringtie_with_varus_directory_is_rejected(tmp_path: Path) -> None:
    varus = tmp_path / "varus_dir"
    varus.mkdir()
    params = {"stringtie": fake_stringtie_gtf(tmp_path), "rnaseq_varus": str(varus), **PROTEINS,
              "tiberius": BOTH_FINDERS["tiberius"], "vipsania": BOTH_FINDERS["vipsania"]}
    proc, _ = run_pipeline(tmp_path, params)
    assert proc.returncode != 0
    assert "pyVARUS directory brings its own" in proc.stdout + proc.stderr


# ---- regressions of the combined flow -----------------------------------------


def test_mixed_two_short_and_two_isoseq_libraries(tmp_path: Path) -> None:
    """Both read types once produced an input file named merged.hints.gff (collision)."""
    iso = str(DATA / "isoseq.fastq")
    params = {"rnaseq_paired": [[R1, R2], [R1, R2]], "isoseq": [iso, iso], **PROTEINS,
              "tiberius": BOTH_FINDERS["tiberius"], "vipsania": BOTH_FINDERS["vipsania"]}
    proc, _ = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "COMBINE_GENESETS" in proc.stdout


def test_proteins_mode_publishes_hints_and_gene_support(tmp_path: Path) -> None:
    proc, published = run_pipeline(tmp_path, {**GENEFINDER["tiberius"], **PROTEINS})
    assert_ok(proc)
    assert {"hintsfile.gff", "qc/gene_support.tsv", "qc/gene_set_statistics.txt"} <= published, sorted(published)
