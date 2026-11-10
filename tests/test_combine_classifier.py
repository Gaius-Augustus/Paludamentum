"""bin/combine_genesets.py --classifier: learned selection with the shipped models conf/combiner_classifier_v1.json
(evidence modes) and conf/combiner_classifier_protonly_v1.json (mode proteins: no RNA-Seq features).
Unit tests of the tree predictor (against probabilities of the scikit-learn models that produced the JSON files)
and of the per-locus selection, plus runs of the script on a tiny synthetic set of four tool GTFs."""
from __future__ import annotations

import csv
import math
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "bin"
SCRIPT = BIN / "combine_genesets.py"
MODEL = ROOT / "conf" / "combiner_classifier_v1.json"
MODEL_PROTEINS = ROOT / "conf" / "combiner_classifier_protonly_v1.json"
MODEL_VIPSANIA = ROOT / "conf" / "combiner_classifier_vipsania_v1.json"
MODEL_ABINITIO = ROOT / "conf" / "combiner_classifier_abinitio_v1.json"
FIXTURES = ROOT / "tests" / "data" / "combine_classifier"
FIXTURE = FIXTURES / "fixture_Tilletiopsis_washingtonensis_300.tsv"
FIXTURE_PROTEINS = FIXTURES / "fixture_Tilletiopsis_washingtonensis_protonly_300.tsv"
FIXTURE_VIPSANIA = FIXTURES / "fixture_Tilletiopsis_washingtonensis_noT_300.tsv"
FIXTURE_ABINITIO = FIXTURES / "fixture_Tilletiopsis_washingtonensis_abinitio_300.tsv"
# model file, its feature count, its fixture (300 structures of a test species never used for training)
MODELS = {"evidence": (MODEL, 47, FIXTURE), "proteins": (MODEL_PROTEINS, 41, FIXTURE_PROTEINS),
          "vipsania": (MODEL_VIPSANIA, 47, FIXTURE_VIPSANIA), "abinitio": (MODEL_ABINITIO, 18, FIXTURE_ABINITIO)}

sys.path.insert(0, str(BIN))
import combine_genesets as cg  # noqa: E402

pytest.importorskip("numpy")


# ---------------------------------------------------------------- model files and predictor

@pytest.mark.parametrize("which", sorted(MODELS))
def test_model_file_is_complete(which):
    model, n_features, _ = MODELS[which]
    clf = cg.load_classifier(str(model))
    assert clf["format"] == cg.CLASSIFIER_FORMAT
    assert len(clf["features"]) == n_features and "clade_code" not in clf["features"]
    assert len(clf["trees"]) == 400
    assert 0 < clf["t_locus"] <= clf["t_locus_unknown_clade"] < 1 and 0 < clf["t_iso"] < 1
    assert len(clf["trained_on"]) == 15
    for t in clf["trees"]:
        n = len(t["feat"])
        assert all(len(t[k]) == n for k in ("thr", "left", "right", "leaf", "value", "nan_left"))
        assert all(0 <= c < n for c, is_leaf in zip(t["left"] + t["right"], t["leaf"] + t["leaf"]) if not is_leaf)


def test_protein_only_model_has_no_rna_features():
    """The protein-only model is trained on pools without reads: its features are a subset of the evidence
    model's without the RNA-Seq columns (intron support, multiplicity, coverage), which the script still
    computes (NaN or 0 without --introns / --coverage) but the model never reads."""
    clf, clf_p = cg.load_classifier(str(MODEL)), cg.load_classifier(str(MODEL_PROTEINS))
    assert set(clf_p["features"]) < set(clf["features"])
    assert set(clf["features"]) - set(clf_p["features"]) == \
        {"rna_frac", "log_mult", "log_mult_rel", "rna_cov", "novel_rna_frac", "rna_frac_rel"}
    assert clf_p["trained_on"] == clf["trained_on"]
    # the sensitive end of the F1 plateau of the species-out cross-validation (0.20-0.30), not its optimum 0.275 / 0.40
    assert (clf_p["t_locus"], clf_p["t_iso"], clf_p["t_locus_unknown_clade"]) == (0.2, 0.3, 0.3)


@pytest.mark.parametrize("which", sorted(MODELS))
def test_predictor_reproduces_sklearn_probabilities(which):
    """300 structures of a test species (never used for training) with the raw feature values and the
    probability of the scikit-learn model that was exported to the JSON."""
    model, _, fixture = MODELS[which]
    clf = cg.load_classifier(str(model))
    with open(fixture) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 300
    F = [{c: cg._num(r[c]) for c in clf["features"]} for r in rows]
    p = cg.classifier_predict(clf, F)
    diff = max(abs(a - float(r["p"])) for a, r in zip(p, rows))
    assert diff < 1e-6, diff
    # missing evidence (structures without a DIAMOND hit; for the evidence model also single-exon structures
    # without introns) is NaN; the ab initio model has no evidence feature at all
    if which != "abinitio":
        assert any(math.isnan(f["dmd_pident"]) for f in F)
    if which == "evidence":
        assert any(math.isnan(f["rna_frac"]) for f in F)
    assert 0.0 < sum(p) / len(p) < 1.0


def test_predictor_handles_nan_and_inf():
    clf = cg.load_classifier(str(MODEL))
    nan_row = {c: float("nan") for c in clf["features"]}
    inf_row = {c: float("inf") for c in clf["features"]}
    p = cg.classifier_predict(clf, [nan_row, inf_row])
    assert all(0.0 <= x <= 1.0 for x in p)


# ---------------------------------------------------------------- selection

def _s(ix, tools, n_exons=2):
    return {"_ix": ix, "tools": set(tools), "segs": tuple((100 * k, 100 * k + 50) for k in range(n_exons)),
            "rank": 0, "keep": False, "why": ""}


def test_selection_primary_isoforms_and_locus_threshold():
    S = [_s(0, ["tiberius"]), _s(1, ["tiberius", "drusilla"]), _s(2, ["td2"]), _s(3, ["vipsania"])]
    by_cl = {0: [S[0], S[1], S[2]], 1: [S[3]]}
    cg.classifier_select(S, by_cl, [0.4, 0.9, 0.35, 0.2], t_locus=0.25, t_iso=0.3, unselectable=set())
    assert (S[1]["keep"], S[1]["why"], S[1]["rank"]) == (True, "primary", 1)
    assert (S[0]["keep"], S[0]["why"], S[0]["rank"]) == (True, "isoform_classifier", 2)
    assert (S[2]["keep"], S[2]["why"], S[2]["rank"]) == (True, "isoform_classifier", 3)
    assert (S[3]["keep"], S[3]["why"]) == (False, "locus_below_min")
    # the isoform threshold is independent of the locus threshold
    cg.classifier_select(S, by_cl, [0.4, 0.9, 0.35, 0.2], t_locus=0.25, t_iso=0.38, unselectable=set())
    assert S[0]["keep"] and not S[2]["keep"] and S[2]["why"] == "isoform_rejected"


def test_selection_ties_keep_input_order_and_unselectable_never_wins():
    S = [_s(0, ["tiberius"], 1), _s(1, ["drusilla", "td2"]), _s(2, ["td2"])]
    by_cl = {0: S[:]}
    cg.classifier_select(S, by_cl, [0.6, 0.6, 0.1], 0.25, 0.3, unselectable=set())
    assert S[0]["why"] == "primary" and S[1]["why"] == "isoform_classifier"
    cg.classifier_select(S, by_cl, [0.6, 0.6, 0.1], 0.25, 0.3, unselectable={0})
    assert (S[0]["keep"], S[0]["why"]) == (False, "dropped_lone_single_exon")
    assert S[1]["why"] == "primary"
    # a locus whose only candidate is unselectable is dropped
    cg.classifier_select(S[:1], {0: S[:1]}, [0.99], 0.25, 0.3, unselectable={0})
    assert not S[0]["keep"]


def test_chain_adjust_moves_log_odds():
    assert (cg.CHAIN_LOGIT_FULL, cg.CHAIN_LOGIT_PARTIAL) == (1.5, 0.0)      # calibrated defaults
    assert cg.chain_adjust(0.5, "full", 1.0, 0.0) == pytest.approx(1 / (1 + math.exp(-1.0)))
    assert cg.chain_adjust(0.5, "partial", 1.0, -1.0) == pytest.approx(1 / (1 + math.exp(1.0)))
    assert cg.chain_adjust(0.3, "none", 2.0, -1.0) == 0.3 and cg.chain_adjust(0.3, "", 2.0, -1.0) == 0.3
    assert cg.chain_adjust(0.3, "full", 0.0, 0.0) == 0.3
    assert 0.0 < cg.chain_adjust(0.0, "full", 3.0, 0.0) < 0.01 and cg.chain_adjust(1.0, "full", -3.0, 0.0) > 0.99


# ---------------------------------------------------------------- script runs on synthetic input

def write_gtf(path: Path, src: str, recs: list) -> Path:
    """recs: (transcript id, [(start, end), ...], strand)."""
    with open(path, "w") as fh:
        for tid, segs, strand in recs:
            for a, b in segs:
                fh.write(f'chr1\t{src}\tCDS\t{a}\t{b}\t.\t{strand}\t0\tgene_id "{tid}"; transcript_id "{tid}";\n')
    return path


@pytest.fixture
def synth(tmp_path: Path) -> Path:
    A = [(1001, 1300), (1401, 1600)]      # four tools agree
    B = [(1001, 1300), (1401, 1700)]      # Tiberius alone, alternative end
    C = [(3001, 3600)]                    # lone single-exon Tiberius model
    D = [(5001, 5300), (5401, 5900)]      # Vipsania alone
    write_gtf(tmp_path / "tiberius.gtf", "Tiberius", [("tA", A, "+"), ("tB", B, "+"), ("tC", C, "+")])
    write_gtf(tmp_path / "vipsania.gtf", "Vipsania", [("vA", A, "+"), ("vD", D, "+")])
    write_gtf(tmp_path / "drusilla.gtf", "Drusilla", [("dA", A, "+")])
    write_gtf(tmp_path / "td2.gtf", "TD2", [("xA", A, "+")])
    return tmp_path


def combine(synth: Path, name: str, *extra: str, check: bool = True):
    out, rep = synth / f"{name}.gtf", synth / f"{name}.report.tsv"
    cmd = [sys.executable, str(SCRIPT), "--tiberius", str(synth / "tiberius.gtf"),
           "--vipsania", str(synth / "vipsania.gtf"), "--drusilla", str(synth / "drusilla.gtf"),
           "--td2", str(synth / "td2.gtf"), "-o", str(out), "--report", str(rep), *extra]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if not check:
        return r
    assert r.returncode == 0, r.stderr
    with open(rep) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    by_id = {row["member_ids"].split(",")[0].split(":")[1]: row for row in rows}
    kept = {m.group(1) for m in (cg.GTF_TID.search(line) for line in open(out)) if m}
    return by_id, kept, r.stderr


def test_classifier_mode_writes_probabilities_and_selects_by_thresholds(synth):
    rows, kept, err = combine(synth, "clf", "--classifier", str(MODEL))
    assert set(rows) == {"tA", "tB", "tC", "vD"}
    assert all(0.0 <= float(r["p_classifier"]) <= 1.0 for r in rows.values())
    assert "classifier:" in err and "t_locus=0.25" in err
    # the four-tool structure is the most probable one and is the primary of its locus
    pa = float(rows["tA"]["p_classifier"])
    assert pa > float(rows["tB"]["p_classifier"]) and pa > float(rows["vD"]["p_classifier"])
    assert rows["tA"]["reason"] == "primary" and rows["tA"]["kept"] == "1"
    for r in rows.values():
        assert r["reason"] in {"primary", "isoform_classifier", "isoform_rejected", "locus_below_min"}
        assert (r["kept"] == "1") == (r["reason"] in {"primary", "isoform_classifier"})
    assert len(kept) == sum(r["kept"] == "1" for r in rows.values())


def test_classifier_thresholds_override_and_unknown_clade(synth):
    rows, kept, _ = combine(synth, "none", "--classifier", str(MODEL), "--t-locus", "1.01")
    assert not kept and all(r["reason"] == "locus_below_min" for r in rows.values())
    rows, kept, _ = combine(synth, "all", "--classifier", str(MODEL), "--t-locus", "0", "--t-iso", "0")
    assert kept == {"g1.t1", "g1.t2", "g2.t1", "g3.t1"}
    assert rows["tB"]["reason"] == "isoform_classifier"
    _, _, err = combine(synth, "unk", "--classifier", str(MODEL), "--unknown-clade")
    assert "t_locus=0.3 " in err and "unknown_clade=True" in err


def test_abinitio_model_without_any_evidence(synth):
    """The ab initio combination: Tiberius and Vipsania alone, no evidence option at all; the evidence-free
    model (tools, structure, locus context) selects at 0.30 / 0.25 and ranks the two-tool structure first."""
    clf = cg.load_classifier(str(MODEL_ABINITIO))
    assert not any(c.startswith(("dmd_", "prot_", "rna_", "mp_")) or c in ("log_mult", "log_mult_rel", "pif_rel", "cl_has_hit")
                   for c in clf["features"])
    assert (clf["t_locus"], clf["t_iso"], clf["t_locus_unknown_clade"]) == (0.3, 0.25, 0.3)
    (synth / "none.gtf").write_text("")
    out, rep = synth / "ab.gtf", synth / "ab.report.tsv"
    cmd = [sys.executable, str(SCRIPT), "--tiberius", str(synth / "tiberius.gtf"), "--vipsania", str(synth / "vipsania.gtf"),
           "--drusilla", str(synth / "none.gtf"), "--td2", str(synth / "none.gtf"),
           "--classifier", str(MODEL_ABINITIO), "-o", str(out), "--report", str(rep)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "t_locus=0.3 t_iso=0.25" in r.stderr
    with open(rep) as fh:
        rows = {row["member_ids"].split(",")[0].split(":")[1]: row for row in csv.DictReader(fh, delimiter="\t")}
    assert set(rows) == {"tA", "tB", "tC", "vD"}
    pa = float(rows["tA"]["p_classifier"])
    assert pa > float(rows["tB"]["p_classifier"]) and pa > float(rows["vD"]["p_classifier"])
    assert rows["tA"]["reason"] == "primary"
    assert all(rows[k]["dmd_searched"] == "" for k in rows)       # no DIAMOND table was given


def test_vipsania_model_without_tiberius_set(synth):
    """Evidence modes without a Tiberius model: the Tiberius slot and its DIAMOND table are empty files; the
    model trained without Tiberius selects among Vipsania, Drusilla and TD2 at 0.225 / 0.30, also for a clade
    outside the HC table (--unknown-clade)."""
    clf = cg.load_classifier(str(MODEL_VIPSANIA))
    assert clf["features"] == cg.load_classifier(str(MODEL))["features"]
    assert (clf["t_locus"], clf["t_iso"], clf["t_locus_unknown_clade"]) == (0.225, 0.3, 0.225)
    (synth / "none.gtf").write_text("")
    (synth / "none.tsv").write_text("")
    out, rep = synth / "vip.gtf", synth / "vip.report.tsv"
    cmd = [sys.executable, str(SCRIPT), "--tiberius", str(synth / "none.gtf"), "--diamond-tiberius", str(synth / "none.tsv"),
           "--diamond-format", "full", "--vipsania", str(synth / "vipsania.gtf"), "--drusilla", str(synth / "drusilla.gtf"),
           "--td2", str(synth / "td2.gtf"), "--classifier", str(MODEL_VIPSANIA), "--unknown-clade",
           "-o", str(out), "--report", str(rep)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "input transcripts with CDS: tiberius=0, vipsania=2, drusilla=1, td2=1" in r.stderr
    assert "t_locus=0.225 t_iso=0.3 unknown_clade=True" in r.stderr
    with open(rep) as fh:
        rows = {row["member_ids"].split(",")[0].split(":")[1]: row for row in csv.DictReader(fh, delimiter="\t")}
    assert set(rows) == {"vA", "vD"}      # the Tiberius-only structures tB and tC are gone
    assert rows["vA"]["tools"] == "drusilla,td2,vipsania" and rows["vA"]["reason"] == "primary"
    assert all(0.0 <= float(row["p_classifier"]) <= 1.0 for row in rows.values())


def test_protein_only_model_on_empty_orf_caller_sets(synth):
    """Mode proteins: the Drusilla and TD2 slots and their DIAMOND tables are empty files, no --introns and
    no --coverage; the protein-only model selects among the Tiberius and Vipsania structures."""
    (synth / "none.gtf").write_text("")
    (synth / "none.tsv").write_text("")
    out, rep = synth / "prot.gtf", synth / "prot.report.tsv"
    cmd = [sys.executable, str(SCRIPT), "--tiberius", str(synth / "tiberius.gtf"),
           "--vipsania", str(synth / "vipsania.gtf"), "--drusilla", str(synth / "none.gtf"),
           "--td2", str(synth / "none.gtf"), "--diamond-drusilla", str(synth / "none.tsv"),
           "--diamond-td2", str(synth / "none.tsv"), "--diamond-format", "full",
           "--classifier", str(MODEL_PROTEINS), "-o", str(out), "--report", str(rep)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "input transcripts with CDS: tiberius=3, vipsania=2, drusilla=0, td2=0" in r.stderr
    assert "t_locus=0.2 t_iso=0.3" in r.stderr
    with open(rep) as fh:
        rows = {row["member_ids"].split(",")[0].split(":")[1]: row for row in csv.DictReader(fh, delimiter="\t")}
    assert set(rows) == {"tA", "tB", "tC", "vD"}
    assert all(0.0 <= float(row["p_classifier"]) <= 1.0 for row in rows.values())
    assert rows["tA"]["tools"] == "tiberius,vipsania" and rows["tA"]["reason"] == "primary"
    # RNA columns stay empty without reads; the model does not read them
    assert rows["tA"]["rna_min_mult"] == "" and rows["tA"]["rna_cov"] == ""


def test_classifier_mode_keeps_the_insect_flag(synth):
    rows, kept, err = combine(synth, "ins", "--classifier", str(MODEL), "--t-locus", "0", "--t-iso", "0",
                              "--drop-lone-single-exon", "tiberius")
    assert rows["tC"]["reason"] == "dropped_lone_single_exon" and rows["tC"]["kept"] == "0"
    assert rows["tA"]["kept"] == "1" and rows["tB"]["kept"] == "1"    # multi-exon Tiberius models are unaffected
    assert kept == {"g1.t1", "g1.t2", "g2.t1"}
    assert "drop_lone_single_exon=tiberius" in err


def test_classifier_chain_correction_raises_supported_structures(synth):
    """A full-length transcript holding the CDS chain of tA raises its P by --chain-logit-full; tB shares the
    intron but ends beyond the transcript (partial) and vD is untouched (none); weight 0 changes nothing."""
    with open(synth / "longreads.gtf", "w") as fh:
        for a, b in [(900, 1300), (1401, 1650)]:
            fh.write(f'chr1\tStringTie\texon\t{a}\t{b}\t.\t+\t.\tgene_id "L"; transcript_id "L.1";\n')
    base, _, _ = combine(synth, "nochain", "--classifier", str(MODEL))
    same, _, err0 = combine(synth, "w0", "--classifier", str(MODEL), "--chains", str(synth / "longreads.gtf"),
                            "--chain-logit-full", "0", "--chain-logit-partial", "0")
    assert (same["tA"]["chain_support"], same["tB"]["chain_support"], same["vD"]["chain_support"]) == \
        ("full", "partial", "none")
    assert all(same[k]["p_classifier"] == base[k]["p_classifier"] for k in base)
    assert "chain_logit_full=0.0" in err0
    up, _, _ = combine(synth, "w2", "--classifier", str(MODEL), "--chains", str(synth / "longreads.gtf"),
                       "--chain-logit-full", "2", "--chain-logit-partial", "-1")
    pa0, pa2 = float(base["tA"]["p_classifier"]), float(up["tA"]["p_classifier"])
    assert pa2 > pa0 and pa2 == pytest.approx(cg.chain_adjust(pa0, "full", 2.0, -1.0), abs=2e-4)
    pb0, pb2 = float(base["tB"]["p_classifier"]), float(up["tB"]["p_classifier"])
    assert pb2 < pb0 and pb2 == pytest.approx(cg.chain_adjust(pb0, "partial", 2.0, -1.0), abs=2e-4)
    assert up["vD"]["p_classifier"] == base["vD"]["p_classifier"]


def test_classifier_refuses_a_tier_table_and_rules_mode_is_unchanged(synth):
    table = synth / "tiers.tsv"
    table.write_text("tools\tstruct_class\ttier\ntiberius\t*\t5\n")
    r = combine(synth, "bad", "--classifier", str(MODEL), "--tier-table", str(table), check=False)
    assert r.returncode != 0 and "exclude each other" in r.stderr
    rows, kept, _ = combine(synth, "rules")
    assert "p_classifier" not in next(iter(rows.values()))
    assert rows["tA"]["reason"] == "primary" and rows["tA"]["locus_tier"] == "10"
