"""Stub runs of the post-processing (subworkflows/postprocess.nf): which
processes run and which files are published for each switch.

Needs a ``nextflow`` executable (or NEXTFLOW_BIN), like tests/test_stub_run.py.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from test_stub_run import (DATA, EVIDENCE, GENEFINDER, NEXTFLOW, assert_ok, citations, expected_outputs,
                           run_pipeline)

pytestmark = [
    pytest.mark.nextflow,
    pytest.mark.skipif(NEXTFLOW is None, reason="nextflow executable not found"),
]

TIBERIUS_RNASEQ = {**GENEFINDER["tiberius"], **EVIDENCE["rnaseq"]}


def busco_cache(tmp_path: Path, lineage: str = "eukaryota_odb12", ready: bool = True) -> str:
    """A lineage cache as DOWNLOAD_BUSCO_LINEAGE leaves it (or an empty one)."""
    cache = tmp_path / "busco"
    lineages = cache / "lineages"
    (lineages / lineage).mkdir(parents=True)
    if ready:
        odb = lineage.split("_")[1]
        (lineages / lineage / "dataset.cfg").write_text("name=x\n")
        for marker in (f"{lineage}.done", f"eukaryota_{odb}.done", "placement_files.done", "file_versions.tsv.done"):
            (lineages / marker).write_text("")
    return str(cache)


def test_postprocess_abinitio_publishes_proteins(tmp_path: Path) -> None:
    """The ab initio mode gets the same final files as the evidence modes, with its own stem."""
    proc, published = run_pipeline(tmp_path, GENEFINDER["tiberius"])
    assert_ok(proc)
    assert {"tiberius_ab_initio.gff3", "tiberius_ab_initio_proteins.fa", "tiberius_ab_initio.gtf",
            "tiberius_ab_initio_cds.fa", "intermediate/tiberius_ab_initio.gff3"} <= published
    # no hints, no transcripts: neither gene support nor UTRs
    assert "GENE_SUPPORT" not in proc.stdout and "ADD_UTRS" not in proc.stdout
    # without fantasia.run there is no GPU probe
    assert "FANTASIA_GPU_CHECK" not in proc.stdout
    assert "Completeness assessment off (no qc.busco_lineage)" in proc.stdout


@pytest.mark.parametrize("mode", ["proteins", "rnaseq"])
def test_utrs_only_with_transcripts(mode: str, tmp_path: Path) -> None:
    proc, published = run_pipeline(tmp_path, {**GENEFINDER["tiberius"], **EVIDENCE[mode]})
    assert_ok(proc)
    assert ("ADD_UTRS" in proc.stdout) == (mode == "rnaseq")
    assert ("qc/utr_report.tsv" in published) == (mode == "rnaseq")


def test_utrs_can_be_switched_off(tmp_path: Path) -> None:
    proc, published = run_pipeline(tmp_path, {**TIBERIUS_RNASEQ, "postprocess": {"utr": False}})
    assert_ok(proc)
    assert "ADD_UTRS" not in proc.stdout
    assert published == expected_outputs("tiberius", "rnaseq") - {"qc/utr_report.tsv"} | {"citations.md"}


@pytest.mark.parametrize("ready", [True, False], ids=["cached", "download"])
def test_completeness_with_a_lineage(ready: bool, tmp_path: Path) -> None:
    qc = {"busco_lineage": "eukaryota", "busco_download_path": busco_cache(tmp_path, ready=ready)}
    proc, published = run_pipeline(tmp_path, {**TIBERIUS_RNASEQ, "qc": qc})
    assert_ok(proc)
    assert ("DOWNLOAD_BUSCO_LINEAGE" in proc.stdout) != ready
    assert {
        "qc/compleasm_genome/summary.txt", "qc/compleasm_proteins/summary.txt",
        "qc/busco_genome_short_summary.txt", "qc/busco_proteins_short_summary.txt", "qc/completeness.tsv",
        "qc/tiberius_evidence_longest_isoform.gff3", "qc/tiberius_evidence_longest_isoform_proteins.fa",
    } <= published
    text = citations(tmp_path)
    assert "**BUSCO**" in text and "**compleasm**" in text


def test_completeness_compleasm_only(tmp_path: Path) -> None:
    qc = {"busco_lineage": "eukaryota_odb12", "busco": False, "busco_download_path": busco_cache(tmp_path)}
    proc, published = run_pipeline(tmp_path, {**GENEFINDER["tiberius"], "qc": qc})
    assert_ok(proc)
    assert "BUSCO_GENOME" not in proc.stdout and "COMPLEASM_GENOME" in proc.stdout
    assert {"qc/completeness.tsv", "qc/compleasm_genome/summary.txt"} <= published
    assert "qc/busco_genome_short_summary.txt" not in published
    assert "**BUSCO**" not in citations(tmp_path)


def rfam_dir(tmp_path: Path) -> str:
    rfam = tmp_path / "rfam"
    rfam.mkdir()
    (rfam / "Rfam.cm").write_text("")
    (rfam / "Rfam.clanin").write_text("")
    return str(rfam)


NCRNA_FILES = {"ncrna/rRNA.gff3", "ncrna/tRNAs.gff3", "ncrna/tRNAs.txt", "ncrna/infernal.tblout",
               "ncrna/ncRNAs_infernal.gff3"}


@pytest.mark.parametrize("mode", ["abinitio", "rnaseq"])
def test_ncrna(mode: str, tmp_path: Path) -> None:
    params = {**GENEFINDER["tiberius"], **EVIDENCE[mode], "ncrna": {"run": True, "rfam_dir": rfam_dir(tmp_path)}}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    stem = "tiberius_ab_initio" if mode == "abinitio" else "tiberius_evidence"
    assert NCRNA_FILES | {f"{stem}_with_ncRNA.gff3"} <= published
    assert "DOWNLOAD_RFAM" not in proc.stdout
    # one cmscan task per genome chunk (the stub of SPLIT_GENOME makes two),
    # on the chunks of the gene finder: the genome is split once
    assert proc.stdout.count("CMSCAN") >= 2
    assert proc.stdout.count("SPLIT_GENOME") == 1, proc.stdout
    lnc = {"ncrna/lncRNAs.gff3", "ncrna/feelnc_classifier.txt"}
    assert (lnc <= published) == (mode == "rnaseq")
    text = citations(tmp_path)
    for tool in ("tRNAscan-SE", "Infernal", "Rfam", "barrnap"):
        assert f"**{tool}" in text, text
    assert ("**FEELnc**" in text) == (mode == "rnaseq")


def test_ncrna_with_an_existing_result_splits_the_genome(tmp_path: Path) -> None:
    """A reused prediction reads no chunks; the genome is split for cmscan alone."""
    result = tmp_path / "previous.gff3"
    result.write_text("##gff-version 3\n")
    params = {"tiberius": {"run": True, "result": str(result)}, "ncrna": {"run": True, "rfam_dir": rfam_dir(tmp_path)}}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "RUN_TIBERIUS" not in proc.stdout
    assert proc.stdout.count("SPLIT_GENOME") == 1, proc.stdout
    assert proc.stdout.count("CMSCAN") >= 2
    assert NCRNA_FILES | {"tiberius_ab_initio_with_ncRNA.gff3"} <= published


def test_ncrna_rfam_dir_without_rfam_files(tmp_path: Path) -> None:
    (tmp_path / "empty").mkdir()
    params = {**GENEFINDER["tiberius"], "ncrna": {"run": True, "rfam_dir": str(tmp_path / "empty")}}
    proc, _ = run_pipeline(tmp_path, params)
    assert proc.returncode != 0
    assert "no Rfam.cm and Rfam.clanin" in proc.stdout + proc.stderr


def fantasia_dirs(tmp_path: Path) -> dict:
    for name in ("hf", "lookup"):
        (tmp_path / name).mkdir()
    return {"run": True, "hf_cache_dir": str(tmp_path / "hf"), "lookup_dir": str(tmp_path / "lookup")}


@pytest.mark.parametrize("ncrna", [False, True], ids=["coding", "with-ncrna"])
def test_fantasia(ncrna: bool, tmp_path: Path) -> None:
    params = {**TIBERIUS_RNASEQ, "fantasia": fantasia_dirs(tmp_path)}
    if ncrna:
        params["ncrna"] = {"run": True, "rfam_dir": rfam_dir(tmp_path)}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert {"tiberius_evidence_go.gff3", "qc/fantasia/results.csv", "qc/fantasia/fantasia_summary.txt",
            "qc/fantasia/fantasia_go_terms.tsv", "qc/fantasia/fantasia_go_categories.png"} <= published
    assert ("tiberius_evidence_with_ncRNA_go.gff3" in published) == ncrna
    # the GPU probe that fails a run without a usable GPU at the start
    assert "FANTASIA_GPU_CHECK" in proc.stdout
    assert "**FANTASIA" in citations(tmp_path)


def test_fantasia_needs_its_directories(tmp_path: Path) -> None:
    proc, _ = run_pipeline(tmp_path, {**GENEFINDER["tiberius"], "fantasia": {"run": True}})
    assert proc.returncode != 0
    assert "fantasia.run = true needs fantasia.hf_cache_dir" in proc.stdout + proc.stderr


@pytest.mark.parametrize("ready", [True, False], ids=["cached", "download"])
def test_omark_and_gffcompare(ready: bool, tmp_path: Path) -> None:
    """The NCBI taxonomy (ete3's taxa.sqlite) is built once on the submitting host."""
    db = tmp_path / "LUCA.h5"
    db.write_text("")
    reference = tmp_path / "reference.gff3"
    reference.write_text("##gff-version 3\n")
    taxa = tmp_path / "ncbi_taxonomy"
    taxa.mkdir()
    if ready:
        (taxa / "taxa.sqlite").write_text("")
    qc = {"omark": True, "omamer_db": str(db), "reference_annotation": str(reference), "ete_taxa_path": str(taxa)}
    proc, published = run_pipeline(tmp_path, {**TIBERIUS_RNASEQ, "qc": qc})
    assert_ok(proc)
    assert ("DOWNLOAD_NCBI_TAXONOMY" in proc.stdout) != ready
    assert {"qc/omark_summary.txt", "qc/gffcompare.stats"} <= published
    text = citations(tmp_path)
    assert "**OMArk**" in text and "GffCompare" in text


def test_omark_needs_the_database(tmp_path: Path) -> None:
    proc, _ = run_pipeline(tmp_path, {**GENEFINDER["tiberius"], "qc": {"omark": True}})
    assert proc.returncode != 0
    assert "qc.omark = true needs the OMAmer database" in proc.stdout + proc.stderr


def test_report_and_statistics_can_be_switched_off(tmp_path: Path) -> None:
    qc = {"report": False, "statistics": False, "gene_support": False}
    proc, published = run_pipeline(tmp_path, {**TIBERIUS_RNASEQ, "qc": qc})
    assert_ok(proc)
    assert "report.html" not in published
    assert not [f for f in published if f.endswith(".png")]
    assert "qc/gene_support.tsv" not in published
    assert {"tiberius_evidence.gff3", "tiberius_evidence.gtf"} <= published
