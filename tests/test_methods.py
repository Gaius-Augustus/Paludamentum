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

from test_stub_run import (BOTH_FINDERS, EVIDENCE, GENEFINDER, NEXTFLOW, assert_ok, run_pipeline)

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


def test_stub_run_proteins_tiberius(tmp_path: Path) -> None:
    cfg = tmp_path / "vertebrates.yaml"
    cfg.write_text('target_species: "Vertebrata"\n')
    params = {"tiberius": {"run": True, "model_cfg": str(cfg)}, **EVIDENCE["proteins"]}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "methods.md" in published
    methods = check_stub_run(tmp_path)
    assert "in mode `proteins`" in methods and "Vertebrata" in methods
    assert "StringTie" not in methods


def test_stub_run_ab_initio_vipsania(tmp_path: Path) -> None:
    proc, _ = run_pipeline(tmp_path, GENEFINDER["vipsania"])
    assert_ok(proc)
    methods = check_stub_run(tmp_path)
    assert "in mode `abinitio`" in methods
    assert "Vipsania" in methods
    assert "miniprot" not in methods and "Tiberius" not in methods
