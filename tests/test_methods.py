"""methods.md: the methods text of a run (lib_nf/methods.nf).

Two kinds of tests, both need a ``nextflow`` executable (or NEXTFLOW_BIN):

- tests/methods_check.nf renders the text of some thousand run maps (every
  mode, gene finder, HC flow, input form and post-processing switch) in one
  Nextflow run without processes. Every citation of a text must be the short
  form of a reference that citations.md of that run lists, and the text must
  keep the format that bin/paludamentum_report.py reads.
- Stub runs check that main.nf writes methods.md, and that its citations are
  those of the citations.md next to it.
"""
from __future__ import annotations

import itertools
import json
import os
import re
import subprocess
from pathlib import Path

import pytest

from test_stub_run import (EVIDENCE, GENEFINDER, NEXTFLOW, ROOT, assert_ok, drusilla_params,
                           run_pipeline)

pytestmark = [
    pytest.mark.nextflow,
    pytest.mark.skipif(NEXTFLOW is None, reason="nextflow executable not found"),
]

TITLE = "# Methods of this Paludamentum run"


# ---------------------------------------------------------------- short citations

def short_cite(ref: str) -> str:
    """The short citation of a reference string, written independently of shortCite in lib_nf/methods.nf."""
    end = ref.find(". ")
    authors = [a.strip() for a in (ref[:end] if end > 0 else ref).split(",") if a.strip()]
    et_al = authors[-1] == "et al"
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


CITATION = re.compile(r"\(([^()]*?(?:\b(?:19|20)\d\d|n\.d\.))\)")


def citations_of(text: str) -> list[str]:
    """The single citations of a methods text: '(A, 2020; B et al., 2021)' -> two."""
    return [c.strip() for group in CITATION.findall(text) for c in group.split(";")]


def citations_md_shorts(text: str) -> set[str]:
    """Short citations of the references that a citations.md lists ('- **Tool** (use). REF LINK')."""
    shorts = set()
    for line in text.splitlines():
        if not line.startswith("- **"):
            continue
        ref = line[line.index("). ") + 3:].rsplit(" ", 1)[0]
        shorts.add(short_cite(ref))
    return shorts


def check_format(text: str) -> None:
    """The format that bin/paludamentum_report.py reads, and the style rules of the text."""
    assert text.startswith(TITLE + "\n\n"), text[:200]
    assert text.endswith("\n") and not text.endswith("\n\n"), repr(text[-50:])
    paragraphs = text[len(TITLE) + 2:-1].split("\n\n")
    assert paragraphs and all(p and "\n" not in p for p in paragraphs), text
    for p in paragraphs:
        assert not p.startswith(("#", "- ", "* ", "1. ")), p
    for bad in ("null", "None", "[]", "${", "http", "](", "  ", " .", " ,", "..", "unknown"):
        assert bad not in text, (bad, text)
    # "ab initio" is always italic
    assert all(m == "*ab initio*" for m in re.findall(r"\*?ab initio\*?", text)), text
    assert text.count("`") % 2 == 0, text


# ---------------------------------------------------------------- run maps

BASE = {
    "mode": "proteins", "genefinder": "tiberius", "tiberiusModel": "/m/vertebrates.yaml",
    "hc": "transdecoder", "orfFinder": "td2", "rescue": True, "odb12": False,
    "shortFastq": False, "shortBam": False, "shortVarus": False, "isoFastq": False, "isoVarus": False,
    "busco": False, "compleasm": False, "omark": False, "gffcompare": False,
    "ncrna": False, "lncrna": True, "fantasia": False,
    "model": "vertebrates", "clade": "Vertebrata", "result": None, "finetune": False,
    "splitMinSize": 20000000, "splitMaxFiles": 20,
    "proteinFiles": 1, "odb12Partitions": [], "scoringMatrix": "blosum62.csv",
    "shortPaired": False, "shortSingle": False, "shortSra": False, "isoSra": False,
    "minAlignmentRate": 80, "stringtieFiles": 0, "mixVarus": False,
    "td2PredictArgs": None, "drusillaSettings": None,
    "sanityFilter": True, "utr": True, "maxUtrExtension": 5000, "buscoLineage": None,
    "geneSupport": True, "statistics": True, "trnascanHighConfidence": False, "fantasiaMinScore": 0.5,
}

DRUSILLA = {
    "model": "vertebrates", "weights": None, "lgbModel": "drusilla_lgb_3class_v1",
    "minLength": 300, "minCov": 3, "minTpm": 1, "longLength": 3000, "minTpmLong": 0.5,
    "minCodingLength": 200, "fixStop": True, "fixStart": True, "lgbThreshold": 0.5, "lgbKeep": ["correct"],
    "rescueModel": "vertebrates", "rescueFlank": 25000, "rescueHintWeight": 2.5, "rescueOrfFilter": False,
}

INPUTS = {
    "abinitio": {"none": {}},
    "proteins": {"local": {}, "odb12": {"odb12": True, "odb12Partitions": ["Vertebrata", "Metazoa"]},
                 "both": {"proteinFiles": 2, "odb12": True, "odb12Partitions": ["Fungi"]}},
    "rnaseq": {
        "fastq": {"shortFastq": True, "shortPaired": True},
        "sra": {"shortFastq": True, "shortSingle": True, "shortSra": True},
        "bam": {"shortBam": True},
        "varus": {"shortVarus": True},
        "stringtie": {"stringtieFiles": 1},
        "fastq+stringtie": {"shortFastq": True, "shortPaired": True, "stringtieFiles": 2},
    },
    "isoseq": {
        "fastq": {"isoFastq": True},
        "sra": {"isoFastq": True, "isoSra": True},
        "varus": {"isoVarus": True},
        "stringtie": {"stringtieFiles": 1},
    },
    "mixed": {
        "fastq": {"shortFastq": True, "shortPaired": True, "shortSingle": True, "isoFastq": True},
        "bam": {"shortBam": True, "isoFastq": True},
        "varus": {"shortVarus": True, "isoVarus": True},
        "mixvarus": {"shortVarus": True, "isoVarus": True, "mixVarus": True},
        "stringtie": {"stringtieFiles": 1},
    },
}

GENEFINDERS = {
    "tiberius": {},
    "tiberius-mammalia": {"tiberiusModel": "/m/mammalia_softmasking_v2.yaml", "model": "mammalia_softmasking_v2",
                          "clade": "Mammalia"},
    "tiberius-result": {"tiberiusModel": None, "model": None, "clade": None, "result": "previous.gtf"},
    "vipsania": {"genefinder": "vipsania", "tiberiusModel": None, "model": "Fungi", "clade": "Fungi"},
    "vipsania-vertebrata": {"genefinder": "vipsania", "tiberiusModel": None, "model": "etb1go6q",
                            "clade": "Vertebrata", "finetune": True},
    "none": {"genefinder": None},
}

POSTPROCESSING = {
    "defaults": {},
    "all": {"busco": True, "compleasm": True, "buscoLineage": "vertebrata_odb12", "omark": True,
            "gffcompare": True, "ncrna": True, "trnascanHighConfidence": True, "fantasia": True},
    "off": {"sanityFilter": False, "utr": False, "geneSupport": False, "statistics": False, "ncrna": True,
            "lncrna": False, "compleasm": True, "buscoLineage": "eukaryota_odb12"},
}


def cases() -> list[tuple[str, dict]]:
    out = []
    for mode, inputs in INPUTS.items():
        transcripts = mode in ("rnaseq", "isoseq", "mixed")
        for (inp, extra), (gf, finder), (pp, post) in itertools.product(
                inputs.items(), GENEFINDERS.items(), POSTPROCESSING.items()):
            if mode == "abinitio" and gf == "none":
                continue   # main.nf stops: the ab initio mode needs a gene finder
            flows = [("transdecoder", "td1"), ("transdecoder", "td2")]
            if transcripts and gf != "none":
                flows += [("drusilla", None)]
            for (hc, orf), rescue in itertools.product(flows, [True, False]):
                if inp == "mixvarus" and hc != "drusilla":
                    continue
                run = {**BASE, "mode": mode, **extra, **finder, **post, "hc": hc,
                       "orfFinder": orf, "rescue": rescue}
                if hc == "drusilla":
                    run["drusillaSettings"] = {**DRUSILLA, "rescueOrfFilter": not rescue,
                                               "lgbKeep": ["correct"] if rescue else ["correct", "partial"]}
                if hc == "transdecoder" and orf == "td2" and rescue:
                    run["td2PredictArgs"] = "--precise"
                out.append((f"{mode}-{inp}-{gf}-{hc}-{orf}-{pp}-rescue{rescue}", run))
    return out


@pytest.fixture(scope="module")
def rendered(tmp_path_factory) -> dict:
    tmp = tmp_path_factory.mktemp("methods")
    named = cases()
    (tmp / "cases.json").write_text(json.dumps([run for _name, run in named]))
    proc = subprocess.run(
        [NEXTFLOW, "run", str(ROOT / "tests" / "methods_check.nf"),
         "--cases", str(tmp / "cases.json"), "--out", str(tmp / "result.json")],
        cwd=tmp, env=dict(os.environ, NXF_ANSI_LOG="false"), capture_output=True, text=True,
    )
    assert_ok(proc)
    result = json.loads((tmp / "result.json").read_text())
    assert len(result["cases"]) == len(named)
    return {"cases": [(name, run, res) for (name, run), res in zip(named, result["cases"])],
            "references": result["references"]}


def test_short_citations_of_the_references(rendered) -> None:
    """shortCite (Groovy) agrees with short_cite (above), and no two references share one."""
    refs = rendered["references"]
    for key, r in refs.items():
        assert r["short"] == short_cite(r["ref"]), key
    shorts = [r["short"] for r in refs.values()]
    assert len(shorts) == len(set(shorts)), sorted(shorts)
    assert refs["miniprot"]["short"] == "Li, 2023"
    assert refs["infernal"]["short"] == "Nawrocki & Eddy, 2013"
    assert refs["nextflow"]["short"] == "Di Tommaso et al., 2017"
    assert refs["transdecoder"]["short"] == "Haas et al., 2013"
    assert refs["paludamentum"]["short"] == "Gabriel & Hoff, n.d."


def test_every_citation_is_a_reference_of_the_run(rendered) -> None:
    refs = rendered["references"]
    assert len(rendered["cases"]) > 1000
    for name, _run, res in rendered["cases"]:
        assert res["problems"] == [], name
        assert set(res["cited"]) <= set(res["selected"]), name
        selected = {refs[k]["short"] for k in res["selected"]}
        cited = citations_of(res["text"])
        assert cited, name
        assert set(cited) <= selected, (name, set(cited) - selected)
        # every cited key appears as its short citation
        assert {refs[k]["short"] for k in res["cited"]} == set(cited), name


def test_format_of_the_text(rendered) -> None:
    for name, _run, res in rendered["cases"]:
        try:
            check_format(res["text"])
        except AssertionError as e:
            raise AssertionError(f"{name}: {e}") from None


def test_text_follows_the_run(rendered) -> None:
    for name, run, res in rendered["cases"]:
        text = res["text"]
        transcripts = run["mode"] in ("rnaseq", "isoseq", "mixed")
        drusilla = transcripts and run["hc"] == "drusilla"
        assert f"in mode `{run['mode']}`" in text, name
        assert ("miniprot" in text) == (run["mode"] != "abinitio"), name
        assert ("Drusilla" in text) == drusilla, name
        assert ("LightGBM" in text) == drusilla, name
        assert ("hint rescue" in text) == (drusilla and run["rescue"]), name
        assert ("TD2" in text) == (transcripts and not drusilla and run["orfFinder"] == "td2"), name
        assert ("TransDecoder" in text) == (transcripts and not drusilla and run["orfFinder"] == "td1"), name
        assert ("StringTie" in text) == transcripts, name
        assert ("HISAT2" in text) == (run["mode"] in ("rnaseq", "mixed") and run["shortFastq"]), name
        assert ("minimap2" in text) == (run["mode"] in ("isoseq", "mixed") and run["isoFastq"]), name
        assert ("Sequence Read Archive" in text) == (
            (run["mode"] in ("rnaseq", "mixed") and (run["shortSra"] or run["shortVarus"]))
            or (run["mode"] in ("isoseq", "mixed") and (run["isoSra"] or run["isoVarus"]))), name
        assert ("OrthoDB v12" in text) == (run["odb12"] and run["mode"] != "abinitio"), name
        gf = run["genefinder"]
        # no gene finder: no final annotation, no post-processing
        assert ("No gene finder was run" in text) == (gf is None), name
        assert ("GFF3" in text) == (gf is not None), name
        if gf is None:
            continue
        assert ("Tiberius" if gf == "tiberius" else "Vipsania") in text, name
        if run["result"]:
            assert f"`{run['result']}`" in text and "chunks" not in text, name
        assert ("UTRs from the StringTie" in text) == (transcripts and run["utr"]), name
        assert ("sanity filter removed" in text) == run["sanityFilter"], name
        assert ("BUSCO" in text) == run["busco"] and ("compleasm" in text) == run["compleasm"], name
        if run["busco"] or run["compleasm"]:
            assert f"`{run['buscoLineage']}`" in text, name
        assert ("OMArk" in text) == run["omark"] and ("GffCompare" in text) == run["gffcompare"], name
        assert ("tRNAscan-SE" in text) == run["ncrna"], name
        assert ("FEELnc" in text) == (run["ncrna"] and run["lncrna"] and transcripts), name
        assert ("FANTASIA" in text) == run["fantasia"], name


# ---------------------------------------------------------------- stub runs

def methods_and_citations(tmp_path: Path) -> tuple[str, str]:
    out = tmp_path / "out"
    return (out / "methods.md").read_text(), (out / "citations.md").read_text()


def check_stub_run(tmp_path: Path) -> str:
    methods, citations = methods_and_citations(tmp_path)
    check_format(methods)
    cited = set(citations_of(methods))
    assert cited and cited <= citations_md_shorts(citations), cited - citations_md_shorts(citations)
    return methods


def test_stub_run_proteins_tiberius(tmp_path: Path) -> None:
    cfg = tmp_path / "vertebrates.yaml"
    cfg.write_text('target_species: "Vertebrata"\n')
    params = {"tiberius": {"run": True, "model_cfg": str(cfg)}, **EVIDENCE["proteins"],
              "qc": {"busco_lineage": "vertebrata"}}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "methods.md" in published
    methods = check_stub_run(tmp_path)
    assert "in mode `proteins`" in methods and "the model `vertebrates`" in methods
    assert "against the lineage `vertebrata_odb12`" in methods
    assert "UTRs" not in methods and "StringTie" not in methods


def test_stub_run_mixed_drusilla(tmp_path: Path) -> None:
    proc, published = run_pipeline(tmp_path, {**drusilla_params(tmp_path, "tiberius"), **EVIDENCE["mixed"],
                                               "ncrna": {"run": True}})
    assert_ok(proc)
    assert "HC genes    : drusilla" in proc.stdout
    methods = check_stub_run(tmp_path)
    for phrase in ("in mode `mixed`", "`stringtie --mix`", "Drusilla flow", "LightGBM", "hint rescue",
                   "HISAT2", "minimap2", "UTRs from the StringTie", "FEELnc", "tRNAscan-SE"):
        assert phrase in methods, (phrase, methods)
    assert "TD2" not in methods and "TransDecoder" not in methods


def test_stub_run_ab_initio_vipsania(tmp_path: Path) -> None:
    proc, published = run_pipeline(tmp_path, GENEFINDER["vipsania"])
    assert_ok(proc)
    methods = check_stub_run(tmp_path)
    assert "in mode `abinitio` (without extrinsic evidence)" in methods
    assert "Vipsania" in methods and "the pretrained model `Fungi`" in methods
    assert "miniprot" not in methods and "Tiberius" not in methods


def test_stub_run_without_gene_finder(tmp_path: Path) -> None:
    proc, published = run_pipeline(tmp_path, {**EVIDENCE["rnaseq"], "transdecoder": "td1"})
    assert_ok(proc)
    methods = check_stub_run(tmp_path)
    assert "No gene finder was run" in methods and "TransDecoder" in methods
    assert "GFF3" not in methods
