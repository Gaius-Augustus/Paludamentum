"""methods.md: the methods text of a run (lib_nf/methods.nf).

Stub runs check that main.nf writes methods.md and that its citations are
those of the citations.md next to it. The exhaustive run-map tests of the
previous flow (methods_check.nf over thousands of maps) are gone with the
unified HC flow; a tighter per-run check lives here instead.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

from test_stub_run import (BOTH_FINDERS, EVIDENCE, GENEFINDER, NEXTFLOW, PROTEINS, assert_ok, fake_stringtie_gtf,
                           run_pipeline)

pytestmark = [
    pytest.mark.nextflow,
    pytest.mark.skipif(NEXTFLOW is None, reason="nextflow executable not found"),
]

TITLE = "# Methods of this Paludamentum run"
CITATION = re.compile(r"\(([^()]*?(?:\b(?:19|20)\d\d[a-z]?|n\.d\.))\)")


def citations_of(text: str) -> list[str]:
    return [c.strip() for group in CITATION.findall(text) for c in group.split(";")]


def short_cite(ref: str) -> str:
    end = ref.find(". ")
    authors = [a.strip() for a in (ref[:end] if end > 0 else ref).split(",") if a.strip()]
    et_al = authors and authors[-1] == "et al"
    if et_al:
        authors = authors[:-1]
    def surname(author: str) -> str:
        parts = author.split()
        return " ".join(parts[:-1]) if len(parts) > 1 else author
    if et_al or len(authors) > 2:
        names = f"{surname(authors[0])} et al."
    elif len(authors) == 2:
        names = f"{surname(authors[0])} & {surname(authors[1])}"
    else:
        names = surname(authors[0])
    years = re.findall(r"\. ((?:19|20)\d\d)[;.]", ref)
    return f"{names}, {years[-1] if years else 'n.d.'}"


def citations_md_shorts(text: str) -> set[str]:
    shorts = set()
    for line in text.splitlines():
        if not line.startswith("- **"):
            continue
        m = re.search(r" \(cited as ([^()]+)\)$", line)
        if m:
            shorts.add(m.group(1))
            continue
        ref = line[line.index("). ") + 3:].rsplit(" ", 1)[0]
        shorts.add(short_cite(ref))
    return shorts


def check_format(text: str) -> None:
    assert text.startswith(TITLE + "\n\n"), text[:200]
    assert text.endswith("\n") and not text.endswith("\n\n"), repr(text[-50:])
    paragraphs = text[len(TITLE) + 2:-1].split("\n\n")
    assert paragraphs and all(p and "\n" not in p for p in paragraphs), text
    for p in paragraphs:
        assert not p.startswith(("#", "- ", "* ", "1. ")), p
    for bad in ("null", "None", "[]", "${", "http", "](", "  ", " .", " ,", "..", "unknown"):
        assert bad not in text, (bad, text)
    assert all(m == "*ab initio*" for m in re.findall(r"\*?ab initio\*?", text)), text
    assert text.count("`") % 2 == 0, text


def methods_and_citations(tmp_path: Path) -> tuple[str, str]:
    out = tmp_path / "out"
    return (out / "methods.md").read_text(), (out / "citations.md").read_text()


def check_stub_run(tmp_path: Path) -> str:
    methods, citations = methods_and_citations(tmp_path)
    check_format(methods)
    cited = set(citations_of(methods))
    assert cited and cited <= citations_md_shorts(citations), cited - citations_md_shorts(citations)
    return methods


def test_stub_run_combined_rnaseq(tmp_path: Path) -> None:
    """Evidence mode: both gene finders, both ORF callers, combine_genesets."""
    params = {"tiberius": BOTH_FINDERS["tiberius"], "vipsania": BOTH_FINDERS["vipsania"], **EVIDENCE["rnaseq"]}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "methods.md" in published
    methods = check_stub_run(tmp_path)
    assert "in mode `rnaseq`" in methods
    for phrase in ("Tiberius", "Vipsania", "StringTie", "miniprot", "TD2", "Drusilla"):
        assert phrase in methods, phrase
    # default selection: the classifier; the stub models' clade (Vipsania model Fungi) has classifier_clade
    # fungi in the HC table, so the t_locus of the classifier model (0.25) applies, not the unknown-clade one
    assert "gradient-boosted tree classifier" in methods and "reached a tier" not in methods
    assert "if its probability reached 0.25," in methods
    assert "outside the training clades" not in methods
    assert "further models of the locus if theirs reached 0.3" in methods


def test_stub_run_combined_clade_outside_the_classifier(tmp_path: Path) -> None:
    """A clade not in the HC table (Stramenopiles) gets the conservative locus threshold 0.3 and says so."""
    params = {"tiberius": BOTH_FINDERS["tiberius"], "vipsania": {"model": "Stramenopiles"}, **EVIDENCE["rnaseq"]}
    proc, _ = run_pipeline(tmp_path, params)
    assert_ok(proc)
    methods = check_stub_run(tmp_path)
    assert "if its probability reached 0.3 (the threshold for clades outside the training clades" in methods
    assert "reached 0.25" not in methods


def test_stub_run_combined_tiberius_mesangiospermae(tmp_path: Path) -> None:
    """Tiberius target_species Mesangiospermae (classifier_clade embryophyta) is a known clade; it wins
    over the Vipsania model (Stramenopiles, not in the table)."""
    cfg = tmp_path / "angiosperms.yaml"
    cfg.write_text('target_species: "Mesangiospermae"\n')
    params = {"tiberius": {"model_cfg": str(cfg)}, "vipsania": {"model": "Stramenopiles"}, **EVIDENCE["rnaseq"]}
    proc, _ = run_pipeline(tmp_path, params)
    assert_ok(proc)
    methods = check_stub_run(tmp_path)
    assert "if its probability reached 0.25," in methods
    assert "outside the training clades" not in methods


def test_stub_run_combined_standard_thresholds_clade(tmp_path: Path) -> None:
    """Chlorophyta is in the HC table with standard_thresholds (no classifier_clade): not a new clade, so the
    t_locus of the classifier model (0.25) applies and the log names standard_thresholds."""
    params = {"tiberius": BOTH_FINDERS["tiberius"], "vipsania": {"model": "Chlorophyta"}, **EVIDENCE["rnaseq"]}
    proc, _ = run_pipeline(tmp_path, params)
    assert_ok(proc)
    methods = check_stub_run(tmp_path)
    assert "if its probability reached 0.25," in methods
    assert "outside the training clades" not in methods
    assert "standard_thresholds" in proc.stdout


def test_stub_run_combined_tiberius_bacillariophyta(tmp_path: Path) -> None:
    """Tiberius target_species Bacillariophyta (standard_thresholds) wins over the Vipsania model
    (Stramenopiles, not in the table): standard thresholds, no unknown-clade note."""
    cfg = tmp_path / "diatoms.yaml"
    cfg.write_text('target_species: "Bacillariophyta"\n')
    params = {"tiberius": {"model_cfg": str(cfg)}, "vipsania": {"model": "Stramenopiles"}, **EVIDENCE["rnaseq"]}
    proc, _ = run_pipeline(tmp_path, params)
    assert_ok(proc)
    methods = check_stub_run(tmp_path)
    assert "if its probability reached 0.25," in methods
    assert "outside the training clades" not in methods


def test_stub_run_combined_rules_and_thresholds(tmp_path: Path) -> None:
    """combine.classifier = false: the tier rules are described; combine.t_locus overrides the model's."""
    base = {"tiberius": BOTH_FINDERS["tiberius"], "vipsania": BOTH_FINDERS["vipsania"], **EVIDENCE["rnaseq"]}
    proc, _ = run_pipeline(tmp_path, {**base, "combine": {"classifier": False}})
    assert_ok(proc)
    methods = check_stub_run(tmp_path)
    assert "reached a tier of 4" in methods
    assert "classifier" not in methods
    other = tmp_path / "t"
    other.mkdir()
    proc, _ = run_pipeline(other, {**base, "combine": {"t_locus": 0.4}})
    assert_ok(proc)
    methods = check_stub_run(other)
    assert "if its probability reached 0.4," in methods


def test_stub_run_rnaseq_without_tiberius(tmp_path: Path) -> None:
    """Vipsania as the only ab initio gene set: three gene sets, the classifier for runs without Tiberius,
    no Tiberius anywhere in the methods or the references."""
    params = {"vipsania": {"model": "Nematoda"}, **EVIDENCE["rnaseq"]}
    proc, _ = run_pipeline(tmp_path, params)
    assert_ok(proc)
    methods = check_stub_run(tmp_path)
    assert "Genes were predicted *ab initio* with Vipsania" in methods
    assert "The three gene sets (Vipsania, revised TD2, Drusilla) were integrated" in methods
    assert "gradient-boosted tree classifier for runs without Tiberius" in methods
    assert "RNA and protein features" in methods
    # thresholds of the model, which also apply to clades outside the HC table (Nematoda)
    assert "if its probability reached 0.225" in methods
    assert "further models of the locus if theirs reached 0.3" in methods
    # Tiberius is named only as the gene set the classifier was trained without, never as a run tool
    for phrase in ("predicted *ab initio* with Tiberius", "(Tiberius, Vipsania", "Tiberius-only", "drop-lone-single-exon"):
        assert phrase not in methods, phrase
    citations = (tmp_path / "out" / "citations.md").read_text()
    assert "**Vipsania**" in citations and "**Tiberius**" not in citations and "gene finder: vipsania" in citations


def test_stub_run_stringtie_without_reads(tmp_path: Path) -> None:
    """A given assembly without reads: the methods say that the intron hints and the coverage track were
    derived from the assembly, and the integration and the gene support count RNA evidence."""
    params = {"tiberius": BOTH_FINDERS["tiberius"], "vipsania": BOTH_FINDERS["vipsania"], **PROTEINS,
              "stringtie": fake_stringtie_gtf(tmp_path)}
    proc, _ = run_pipeline(tmp_path, params)
    assert_ok(proc)
    methods = check_stub_run(tmp_path)
    assert "in mode `rnaseq` (with protein evidence and transcript assemblies)" in methods
    assert "was used as the transcript assembly; no reads were aligned." in methods
    assert "Intron hints were derived from the assembly" in methods
    assert "a coverage track from the coverage values of its exons" in methods
    assert "the intron hints and the coverage track derived from the assembly were used" in methods
    assert "RNA and protein features" in methods
    assert "the hints of the run were the protein hints and the intron hints of the assembly" in methods
    assert "supported by RNA-Seq and protein hints were counted" in methods
    for phrase in ("bam2hints", "genomecov", "HISAT2", "no coverage track and no intron hints"):
        assert phrase not in methods, phrase


def test_stub_run_proteins_combined(tmp_path: Path) -> None:
    """Protein-only mode: both gene finders and combine_genesets.py with the protein-only classifier,
    no transcript steps in the text."""
    cfg = tmp_path / "vertebrates.yaml"
    cfg.write_text('target_species: "Vertebrata"\n')
    params = {"tiberius": {"model_cfg": str(cfg)}, "vipsania": BOTH_FINDERS["vipsania"], **EVIDENCE["proteins"]}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "methods.md" in published
    methods = check_stub_run(tmp_path)
    assert "in mode `proteins` (with protein evidence)" in methods
    for phrase in ("Tiberius", "Vipsania", "miniprot", "DIAMOND", "two *ab initio* gene sets (Tiberius, Vipsania)",
                   "gradient-boosted tree classifier for runs without transcript evidence",
                   "without transcript evidence no ORF caller ran"):
        assert phrase in methods, phrase
    for phrase in ("StringTie", "TD2", "Drusilla", "bedtools", "bam2hints", "RNA and protein features", "UTR"):
        assert phrase not in methods, phrase
    # the protein-only model's thresholds (Vertebrata has classifier_clade vertebrates: no unknown-clade threshold)
    assert "if its probability reached 0.2," in methods
    assert "further models of the locus if theirs reached 0.3" in methods
    assert "outside the training clades" not in methods
    assert "reduced to the proteins of 13 species" in methods
    assert "the genus and species in front of the identifiers of reference proteomes" in methods
    citations = (tmp_path / "out" / "citations.md").read_text()
    for tool in ("Tiberius", "Vipsania", "miniprot", "DIAMOND"):
        assert f"**{tool}**" in citations, tool
    for tool in ("StringTie", "TD2", "Drusilla", "PSAURON"):
        assert f"**{tool}**" not in citations, tool


def test_stub_run_proteins_no_species_reduction(tmp_path: Path) -> None:
    cfg = tmp_path / "vertebrates.yaml"
    cfg.write_text('target_species: "Vertebrata"\n')
    params = {"tiberius": {"model_cfg": str(cfg)}, "vipsania": BOTH_FINDERS["vipsania"],
              "protein_top_species": 0, **EVIDENCE["proteins"]}
    proc, _ = run_pipeline(tmp_path, params)
    assert_ok(proc)
    methods = check_stub_run(tmp_path)
    assert "reduced to the proteins of" not in methods


def test_stub_run_proteins_tier_rules(tmp_path: Path) -> None:
    """combine.classifier_proteins = false: the tier rules are described, without RNA rules."""
    params = {"tiberius": BOTH_FINDERS["tiberius"], "vipsania": BOTH_FINDERS["vipsania"],
              "combine": {"classifier_proteins": False}, **EVIDENCE["proteins"]}
    proc, _ = run_pipeline(tmp_path, params)
    assert_ok(proc)
    methods = check_stub_run(tmp_path)
    assert "reached a tier of 4" in methods and "classifier" not in methods
    assert "adjusted the tier with protein and DIAMOND features" in methods


def test_stub_run_abinitio_both_models(tmp_path: Path) -> None:
    """The ab initio combination: both gene finders and combine_genesets.py without evidence; the methods
    describe the evidence-free classifier at 0.30 / 0.25 and cite no evidence tool."""
    params = {"tiberius": BOTH_FINDERS["tiberius"], "vipsania": BOTH_FINDERS["vipsania"]}
    proc, _ = run_pipeline(tmp_path, params)
    assert_ok(proc)
    methods = check_stub_run(tmp_path)
    assert "in mode `abinitio` (without extrinsic evidence)" in methods
    assert "A second set of *ab initio* predictions was made with Vipsania" in methods
    assert "integrated by `combine_genesets.py` without extrinsic evidence" in methods
    assert "gradient-boosted tree classifier for runs without evidence" in methods
    assert "if its probability reached 0.3, and further models of the locus if theirs reached 0.25" in methods
    for phrase in ("miniprot", "DIAMOND", "StringTie", "TD2", "Drusilla", "hints", "UTR"):
        assert phrase not in methods, phrase
    citations = (tmp_path / "out" / "citations.md").read_text()
    for tool in ("Tiberius", "Vipsania", "GffRead", "GenomeTools"):
        assert f"**{tool}" in citations, tool
    for tool in ("miniprot", "DIAMOND", "StringTie", "TD2", "Drusilla"):
        assert f"**{tool}**" not in citations, tool


def test_stub_run_ab_initio_vipsania(tmp_path: Path) -> None:
    proc, _ = run_pipeline(tmp_path, GENEFINDER["vipsania"])
    assert_ok(proc)
    methods = check_stub_run(tmp_path)
    assert "in mode `abinitio`" in methods
    assert "Vipsania" in methods
    assert "miniprot" not in methods and "Tiberius" not in methods
