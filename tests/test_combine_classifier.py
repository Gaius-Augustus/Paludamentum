"""bin/combine_genesets.py --classifier: learned selection with the shipped models conf/combiner_classifier_v1.json
(evidence modes) and conf/combiner_classifier_protonly_v1.json (mode proteins: no RNA-Seq features).
Unit tests of the tree predictor (against probabilities of the scikit-learn models that produced the JSON files)
and of the per-locus selection, plus runs of the script on a tiny synthetic set of four tool GTFs."""
from __future__ import annotations

import csv
import math
import re
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


# ---------------------------------------------------------------- fragmented assemblies: L50 and the locus threshold

def write_fasta(path: Path, records: list, width: int = 60) -> Path:
    """records: (name, length); sequences wrapped at `width`, a .gz path is written compressed."""
    import gzip
    text = ""
    for name, n in records:
        seq = "ACGT" * (n // 4 + 1)
        seq = seq[:n]
        text += f">{name} toy\n" + "".join(seq[i:i + width] + "\n" for i in range(0, n, width))
    if str(path).endswith(".gz"):
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


@pytest.mark.parametrize("suffix", [".fa", ".fa.gz"])
def test_assembly_stats_n50_l50(tmp_path, suffix):
    # total 1000: 400 + 300 = 700 >= 500 -> L50 2, N50 300 (the shortest of those two)
    g = write_fasta(tmp_path / f"g{suffix}", [("s1", 100), ("s2", 400), ("s3", 50), ("s4", 300), ("s5", 150)])
    assert cg.assembly_stats(str(g)) == {"n_seqs": 5, "length": 1000, "n50": 300, "l50": 2}


def test_assembly_stats_exact_half_and_single_sequence(tmp_path):
    # the longest sequence holds exactly half the length: L50 1 (cumulative >= half)
    g = write_fasta(tmp_path / "half.fa", [("a", 500), ("b", 250), ("c", 250)])
    assert cg.assembly_stats(str(g)) == {"n_seqs": 3, "length": 1000, "n50": 500, "l50": 1}
    one = write_fasta(tmp_path / "one.fa.gz", [("only", 77)])
    assert cg.assembly_stats(str(one)) == {"n_seqs": 1, "length": 77, "n50": 77, "l50": 1}


def test_assembly_stats_empty_record_multiline_and_empty_file(tmp_path):
    g = tmp_path / "e.fa"
    g.write_text(">a\nACGT\nACGT\nAC\n>empty\n>b desc\nACGTACGTAC\n\n")
    st = cg.assembly_stats(str(g))
    assert st["n_seqs"] == 3 and st["length"] == 20            # the empty record counts as a sequence of length 0
    assert (st["n50"], st["l50"]) == (10, 1)                   # 10 of 20 is exactly half
    (tmp_path / "none.fa").write_text("")
    assert cg.assembly_stats(str(tmp_path / "none.fa")) == {"n_seqs": 0, "length": 0, "n50": 0, "l50": 0}


def test_select_t_locus_precedence():
    clf = {"t_locus": 0.25, "t_locus_unknown_clade": 0.3}
    assert cg.select_t_locus(clf) == (0.25, "model")
    assert cg.select_t_locus(clf, unknown_clade=True) == (0.3, "model_unknown_clade")
    assert cg.select_t_locus(clf, unknown_clade=True, t_locus_fragmented=0.1) == (0.1, "fragmented_assembly")
    assert cg.select_t_locus(clf, t_locus=0.4, unknown_clade=True, t_locus_fragmented=0.1) == (0.4, "t_locus")
    assert cg.select_t_locus(clf, t_locus=0.0, t_locus_fragmented=0.1) == (0.0, "t_locus")   # 0 is a value, not unset
    # a model file without the unknown-clade threshold falls back to t_locus
    assert cg.select_t_locus({"t_locus": 0.2}, unknown_clade=True) == (0.2, "model_unknown_clade")


def read_settings(path: Path) -> dict:
    lines = path.read_text().splitlines()
    assert lines[0] == "key\tvalue"
    return dict(line.split("\t", 1) for line in lines[1:])


FRAG = ("--fragmented-l50", "5", "--t-locus-fragmented", "0.1")


def test_fragmented_assembly_lowers_the_locus_threshold(synth):
    """20 sequences of equal length: L50 10 > 5. The synth loci have P 0.68 (tA), 0.20 (tC), 0.18 (vD): the
    model's 0.25 keeps one locus, 0.10 keeps all three."""
    write_fasta(synth / "frag.fa", [("chr1", 1000)] + [(f"s{i}", 1000) for i in range(19)])
    _, kept0, _ = combine(synth, "base", "--classifier", str(MODEL))
    assert len(kept0) == 1
    rows, kept, err = combine(synth, "frag", "--classifier", str(MODEL), "--genome", str(synth / "frag.fa"), *FRAG,
                              "--settings", str(synth / "frag.settings.tsv"))
    assert len(kept) == 3 and rows["tC"]["kept"] == "1" and rows["vD"]["kept"] == "1"
    assert "t_locus=0.1 " in err and "L50 10" in err and "fragmented (L50 > 5)" in err
    s = read_settings(synth / "frag.settings.tsv")
    assert s == {"selection": "classifier", "assembly_sequences": "20", "assembly_length": "20000",
                 "assembly_n50": "1000", "assembly_l50": "10", "fragmented_l50": "5", "fragmented": "1",
                 "t_locus": "0.1", "t_locus_source": "fragmented_assembly", "t_iso": "0.3",
                 "rescue_protein_cov": "", "rescue_protein_introns": "", "rescued_loci": "0",
                 "split_loci": "cds", "split_genes": "0"}


def test_fragmented_threshold_also_for_unknown_clade_and_gz_genome(synth):
    write_fasta(synth / "frag.fa.gz", [(f"s{i}", 1000) for i in range(20)])
    _, kept, _ = combine(synth, "unk", "--classifier", str(MODEL), "--unknown-clade",
                         "--genome", str(synth / "frag.fa.gz"), *FRAG, "--settings", str(synth / "u.tsv"))
    assert len(kept) == 3
    assert read_settings(synth / "u.tsv")["t_locus_source"] == "fragmented_assembly"


def test_contiguous_assembly_keeps_the_model_threshold(synth):
    # L50 1 <= 5; also exactly L50 == N is not fragmented (L50 must exceed N)
    write_fasta(synth / "contig.fa", [("chr1", 20000)] + [(f"s{i}", 100) for i in range(10)])
    _, kept, err = combine(synth, "contig", "--classifier", str(MODEL), "--genome", str(synth / "contig.fa"), *FRAG,
                           "--settings", str(synth / "c.tsv"))
    assert len(kept) == 1 and "fragmented (" not in err and "L50 1," in err
    s = read_settings(synth / "c.tsv")
    assert (s["fragmented"], s["t_locus"], s["t_locus_source"], s["assembly_l50"]) == ("0", "0.25", "model", "1")
    write_fasta(synth / "edge.fa", [(f"s{i}", 1000) for i in range(10)])       # L50 5
    combine(synth, "edge", "--classifier", str(MODEL), "--genome", str(synth / "edge.fa"), *FRAG,
            "--settings", str(synth / "e.tsv"))
    s = read_settings(synth / "e.tsv")
    assert (s["assembly_l50"], s["fragmented"], s["t_locus_source"]) == ("5", "0", "model")


def test_t_locus_overrides_the_fragmented_threshold(synth):
    write_fasta(synth / "frag.fa", [(f"s{i}", 1000) for i in range(20)])
    _, kept, _ = combine(synth, "ov", "--classifier", str(MODEL), "--genome", str(synth / "frag.fa"), *FRAG,
                         "--t-locus", "0.5", "--settings", str(synth / "o.tsv"))
    assert len(kept) == 1
    s = read_settings(synth / "o.tsv")
    assert (s["t_locus"], s["t_locus_source"], s["fragmented"]) == ("0.5", "t_locus", "1")


def test_fragmented_options_are_inert_without_genome_or_both_options(synth):
    write_fasta(synth / "frag.fa", [(f"s{i}", 1000) for i in range(20)])
    combine(synth, "nog", "--classifier", str(MODEL), *FRAG, "--settings", str(synth / "n.tsv"))
    s = read_settings(synth / "n.tsv")
    assert "assembly_l50" not in s and s["fragmented"] == "0" and s["t_locus_source"] == "model"
    combine(synth, "not", "--classifier", str(MODEL), "--genome", str(synth / "frag.fa"), "--fragmented-l50", "5",
            "--settings", str(synth / "t.tsv"))        # no --t-locus-fragmented
    s = read_settings(synth / "t.tsv")
    assert s["assembly_l50"] == "10"
    assert s["fragmented"] == "0" and s["t_locus"] == "0.25"


def test_settings_in_tier_rule_mode_write_min_locus(synth):
    write_fasta(synth / "frag.fa", [(f"s{i}", 1000) for i in range(20)])
    combine(synth, "rules", "--genome", str(synth / "frag.fa"), *FRAG, "--min-locus", "7",
            "--settings", str(synth / "r.tsv"))
    s = read_settings(synth / "r.tsv")
    assert s["selection"] == "tier_rules" and float(s["min_locus"]) == 7.0
    assert "t_locus" not in s and "t_iso" not in s
    # the tier rules do not use the fragmented threshold
    assert s["assembly_l50"] == "10" and s["fragmented"] == "1"


# ---------------------------------------------------------------- protein rescue of loci below the locus threshold

def _r(ix, n_exons=2, qcov=0.95, tcov=0.95, pif=1.0):
    s = _s(ix, ["tiberius"], n_exons)
    s.update({"dmd_qcov": qcov, "dmd_tcov": tcov, "prot_intron_frac": pif})
    return s


def test_protein_rescuable_thresholds_are_inclusive():
    assert cg.protein_rescuable(_r(0), 0.9, 0.5)
    assert cg.protein_rescuable(_r(0, qcov=0.9, tcov=0.9, pif=0.5), 0.9, 0.5)
    assert not cg.protein_rescuable(_r(0, qcov=0.89), 0.9, 0.5)
    assert not cg.protein_rescuable(_r(0, tcov=0.89), 0.9, 0.5)
    assert not cg.protein_rescuable(_r(0, pif=0.49), 0.9, 0.5)
    assert cg.protein_rescuable(_r(0, qcov=0.0, tcov=0.0, pif=0.0), 0.0, 0.0)      # 0 is a value, not missing


def test_protein_rescuable_needs_two_cds_segments():
    assert not cg.protein_rescuable(_r(0, n_exons=1), 0.9, 0.5)
    assert cg.protein_rescuable(_r(0, n_exons=2), 0.9, 0.5) and cg.protein_rescuable(_r(0, n_exons=5), 0.9, 0.5)


@pytest.mark.parametrize("key", ["dmd_qcov", "dmd_tcov", "prot_intron_frac"])
@pytest.mark.parametrize("missing", ["", None, "absent"])
def test_protein_rescuable_missing_evidence_never_qualifies(key, missing):
    s = _r(0)
    if missing == "absent":
        del s[key]
    else:
        s[key] = missing
    assert not cg.protein_rescuable(s, 0.9, 0.5)
    assert not cg.protein_rescuable(s, 0.0, 0.0)


def test_selection_rescues_a_locus_below_the_threshold():
    S = [_s(0, ["tiberius"]), _s(1, ["vipsania"]), _s(2, ["td2"])]
    by_cl = {0: S[:]}
    cg.classifier_select(S, by_cl, [0.20, 0.10, 0.05], 0.25, 0.3, set(), rescuable={1})
    assert [s["_ix"] for s in by_cl[0]] == [1, 0, 2]
    assert [(S[i]["rank"], S[i]["keep"], S[i]["why"]) for i in (1, 0, 2)] == [
        (1, True, "primary_protein_rescue"), (2, False, "isoform_rejected"), (3, False, "isoform_rejected")]
    # the other members follow the usual isoform rule: P >= t_iso
    cg.classifier_select(S, {0: S[:]}, [0.24, 0.10, 0.05], 0.25, 0.2, set(), rescuable={1})
    assert (S[0]["keep"], S[0]["why"]) == (True, "isoform_classifier")
    assert (S[2]["keep"], S[2]["why"]) == (False, "isoform_rejected")


def test_selection_rescue_picks_the_best_ranked_rescuable_member():
    S = [_s(0, ["tiberius"]), _s(1, ["vipsania"]), _s(2, ["td2"]), _s(3, ["drusilla"])]
    by_cl = {0: S[:]}
    cg.classifier_select(S, by_cl, [0.20, 0.10, 0.15, 0.05], 0.25, 0.3, set(), rescuable={1, 2, 3})
    assert [s["_ix"] for s in by_cl[0]] == [2, 0, 1, 3]
    assert (S[2]["rank"], S[2]["why"]) == (1, "primary_protein_rescue")
    assert all(S[i]["why"] == "isoform_rejected" for i in (0, 1, 3)) and [S[i]["rank"] for i in (0, 1, 3)] == [2, 3, 4]
    # equal P: input order decides
    cg.classifier_select(S, {0: S[:]}, [0.2, 0.1, 0.1, 0.1], 0.25, 0.3, set(), rescuable={1, 2})
    assert S[1]["why"] == "primary_protein_rescue" and S[2]["why"] == "isoform_rejected"


def test_selection_unselectable_member_is_not_rescued():
    S = [_s(0, ["tiberius"], 1), _s(1, ["vipsania"]), _s(2, ["td2"])]
    cg.classifier_select(S, {0: S[:]}, [0.20, 0.15, 0.10], 0.25, 0.3, {1}, rescuable={1})
    assert (S[1]["keep"], S[1]["why"]) == (False, "dropped_lone_single_exon")
    assert all(not s["keep"] for s in S) and S[0]["why"] == S[2]["why"] == "locus_below_min"
    # another rescuable member of the locus is rescued instead
    cg.classifier_select(S, {0: S[:]}, [0.20, 0.15, 0.10], 0.25, 0.3, {1}, rescuable={1, 2})
    assert (S[2]["keep"], S[2]["why"], S[2]["rank"]) == (True, "primary_protein_rescue", 1)
    assert S[1]["why"] == "dropped_lone_single_exon" and S[0]["why"] == "isoform_rejected"


def test_selection_open_locus_keeps_reason_primary_and_default_is_unchanged():
    S = [_s(0, ["tiberius"]), _s(1, ["vipsania"])]
    cg.classifier_select(S, {0: S[:]}, [0.6, 0.5], 0.25, 0.3, set(), rescuable={1})
    assert (S[0]["why"], S[0]["rank"]) == ("primary", 1)
    assert (S[1]["why"], S[1]["rank"]) == ("isoform_classifier", 2)
    # below the threshold: without a rescuable member (default argument, empty set, other index) nothing changes
    for kw in ({}, {"rescuable": set()}, {"rescuable": {7}}):
        cg.classifier_select(S, {0: S[:]}, [0.2, 0.1], 0.25, 0.3, set(), **kw)
        assert [(s["keep"], s["why"]) for s in S] == [(False, "locus_below_min")] * 2
    # a rescuable member of another locus does not rescue this one
    T = [_s(0, ["tiberius"]), _s(1, ["vipsania"])]
    cg.classifier_select(T, {0: [T[0]], 1: [T[1]]}, [0.2, 0.1], 0.25, 0.3, set(), rescuable={1})
    assert (T[0]["why"], T[1]["why"]) == ("locus_below_min", "primary_protein_rescue")


# end to end: vD (Vipsania alone, two CDS exons, one intron 5301-5400) and tC (single exon) are below the locus
# threshold (--t-locus 0.99 puts every locus there, whatever protein evidence does to P)

RESCUE = ("--rescue-protein-cov", "0.9", "--rescue-protein-introns", "0.5")


def write_diamond(path: Path, rows: list) -> Path:
    """rows: (query, subject, qstart, qend, sstart, send, qlen, slen) in the 14-column layout ('full')."""
    with open(path, "w") as fh:
        for q, sbj, qs, qe, ss, se, ql, sl in rows:
            fh.write("\t".join(map(str, (q, sbj, 90.0, qe - qs + 1, 0, 0, qs, qe, ss, se, 1e-30, 300.0, ql, sl))) + "\n")
    return path


def write_hints(path: Path, introns: list) -> Path:
    with open(path, "w") as fh:
        for a, b in introns:
            fh.write(f"chr1\tb2h\tintron\t{a}\t{b}\t5\t+\t.\tmult=5;pri=4;src=P\n")
    return path


@pytest.fixture
def rescue_inputs(synth: Path):
    """Full-length DIAMOND hit and protein hints for the intron of vD; tC has a hit but cannot be rescued."""
    dmd = write_diamond(synth / "vip.dmd.tsv", [("vD", "sp1_P1", 1, 270, 1, 280, 280, 290)])
    dmd_t = write_diamond(synth / "tib.dmd.tsv", [("tC", "sp1_P2", 1, 200, 1, 200, 200, 200)])
    return ("--diamond-vipsania", str(dmd), "--diamond-tiberius", str(dmd_t), "--diamond-format", "full",
            "--protein-hints", str(write_hints(synth / "hc.gff", [(5301, 5400)])))


def test_rescue_keeps_a_locus_below_the_threshold_end_to_end(synth, rescue_inputs):
    base = ("--classifier", str(MODEL), "--t-locus", "0.99", *rescue_inputs)
    rows, kept, err = combine(synth, "off", *base, "--settings", str(synth / "off.tsv"))
    assert not kept and all(r["reason"] == "locus_below_min" for r in rows.values())
    s = read_settings(synth / "off.tsv")
    assert (s["rescue_protein_cov"], s["rescue_protein_introns"], s["rescued_loci"]) == ("", "", "0")
    assert "rescue_protein_cov" not in err and "primary_protein_rescue: 0" in err

    rows, kept, err = combine(synth, "on", *base, *RESCUE, "--settings", str(synth / "on.tsv"))
    assert rows["vD"]["reason"] == "primary_protein_rescue" and rows["vD"]["kept"] == "1"
    assert {r["reason"] for k, r in rows.items() if k != "vD"} == {"locus_below_min"}      # tC: single exon
    assert len(kept) == 1
    gtf = (synth / "on.gtf").read_text()
    assert "5001" in gtf and "5900" in gtf and "1001" not in gtf and "3001" not in gtf
    s = read_settings(synth / "on.tsv")
    assert (s["rescue_protein_cov"], s["rescue_protein_introns"], s["rescued_loci"]) == ("0.9", "0.5", "1")
    assert s["t_locus"] == "0.99" and s["selection"] == "classifier"
    assert "rescue_protein_cov=0.9 rescue_protein_introns=0.5" in err and "primary_protein_rescue: 1" in err
    assert "locus_below_min: 3" in err


def test_rescue_needs_both_coverages_and_the_introns(synth):
    hints = write_hints(synth / "hc.gff", [(5301, 5400)])
    full = write_diamond(synth / "a.tsv", [("vD", "sp1_P1", 1, 270, 1, 280, 280, 290)])
    short_q = write_diamond(synth / "b.tsv", [("vD", "sp1_P1", 1, 200, 1, 280, 280, 290)])     # qcov 0.71
    short_t = write_diamond(synth / "c.tsv", [("vD", "sp1_P1", 1, 270, 1, 280, 280, 400)])     # tcov 0.70
    common = ("--classifier", str(MODEL), "--t-locus", "0.99", "--diamond-format", "full")

    def run(name, dmd, hint_args, *opts):
        _, kept, err = combine(synth, name, *common, "--diamond-vipsania", str(dmd), *hint_args, *(opts or RESCUE))
        return kept, err
    assert len(run("ok", full, ("--protein-hints", str(hints)))[0]) == 1
    assert not run("sq", short_q, ("--protein-hints", str(hints)))[0]
    assert not run("st", short_t, ("--protein-hints", str(hints)))[0]
    assert not run("nh", full, ())[0]                                                   # no protein hints
    other = write_hints(synth / "other.gff", [(5302, 5399)])                            # hint of another intron
    assert not run("oh", full, ("--protein-hints", str(other)))[0]
    # only the coverage given: --rescue-protein-introns defaults to 0.5
    kept, err = run("def", full, ("--protein-hints", str(hints)), "--rescue-protein-cov", "0.9")
    assert len(kept) == 1 and "rescue_protein_introns=0.5" in err
    # the hit covers 270/280 = 0.964 of the model and 280/290 of the protein: 0.98 rejects it
    assert not run("strict", full, ("--protein-hints", str(hints)), "--rescue-protein-cov", "0.98")[0]


def test_rescue_does_nothing_for_an_open_locus_or_without_the_classifier(synth, rescue_inputs):
    rows, _, err = combine(synth, "open", "--classifier", str(MODEL), *rescue_inputs, *RESCUE,
                           "--settings", str(synth / "open.tsv"))
    assert rows["tA"]["reason"] == "primary"
    n = sum(r["reason"] == "primary_protein_rescue" for r in rows.values())
    assert read_settings(synth / "open.tsv")["rescued_loci"] == str(n) and f"primary_protein_rescue: {n}" in err
    # tier rules: the options are accepted and have no effect, no rescue rows
    rows, _, err = combine(synth, "rules", *rescue_inputs, *RESCUE, "--settings", str(synth / "rules.tsv"))
    assert all(r["reason"] != "primary_protein_rescue" for r in rows.values())
    assert "rescue_protein_cov" not in read_settings(synth / "rules.tsv") and "classifier:" not in err


# ---------------------------------------------------------------- --split-loci

def _seg(ix, a, b, tools=("tiberius",)):
    """structure with one CDS span a..b (two segments)"""
    return {"_ix": ix, "tools": set(tools), "segs": ((a, a + 10), (b - 10, b)), "rank": 0, "keep": False, "why": ""}


def test_split_selects_disjoint_members_as_genes_of_their_own():
    # 0: long model over 100..900 (fusion), 1: 100..300, 2: 600..900; 0 overlaps both, 1 and 2 are disjoint
    S = [_seg(0, 100, 900, ["vipsania"]), _seg(1, 100, 300), _seg(2, 600, 900)]
    p = [0.2, 0.8, 0.6]
    cg.classifier_select(S, {0: S[:]}, p, 0.25, 0.3, set())
    assert [(s["keep"], s["why"]) for s in S] == [(False, "isoform_rejected"), (True, "primary"),
                                                  (True, "isoform_classifier")]
    cg.classifier_select(S, {0: S[:]}, p, 0.25, 0.3, set(), split=True)
    assert [(s["keep"], s["why"]) for s in S] == [(False, "isoform_rejected"), (True, "primary"),
                                                  (True, "primary_split")]
    assert (S[1]["sub"], S[2]["sub"]) == (0, 1) and (S[1]["rank"], S[2]["rank"]) == (1, 3)


def test_split_keeps_a_long_primary_as_one_gene():
    # the long model is the best one: everything overlaps it, nothing is split off
    S = [_seg(0, 100, 900, ["vipsania"]), _seg(1, 100, 300), _seg(2, 600, 900)]
    cg.classifier_select(S, {0: S[:]}, [0.9, 0.8, 0.2], 0.25, 0.3, set(), split=True)
    assert [(s["keep"], s["why"]) for s in S] == [(True, "primary"), (True, "isoform_classifier"),
                                                  (False, "isoform_rejected")]
    assert {s["sub"] for s in S if s["keep"]} == {0}


def test_split_disjoint_member_below_the_threshold_is_dropped():
    S = [_seg(0, 100, 300), _seg(1, 600, 900)]
    cg.classifier_select(S, {0: S[:]}, [0.8, 0.2], 0.25, 0.1, set(), split=True)
    # without the split it would be an isoform (P >= t_iso 0.1); as a locus of its own it is below t_locus
    assert [(s["keep"], s["why"]) for s in S] == [(True, "primary"), (False, "locus_below_min")]
    cg.classifier_select(S, {0: S[:]}, [0.8, 0.2], 0.25, 0.1, set())
    assert (S[1]["keep"], S[1]["why"]) == (True, "isoform_classifier")


def test_split_member_overlapping_only_an_isoform_stays_with_the_gene():
    # 0 primary 100..300; 1 isoform 200..700 (overlaps 0); 2 at 600..900 overlaps the isoform, not the primary:
    # an isoform of the same gene as without the split (P >= t_iso), never a gene of its own
    S = [_seg(0, 100, 300), _seg(1, 200, 700), _seg(2, 600, 900)]
    cg.classifier_select(S, {0: S[:]}, [0.9, 0.5, 0.8], 0.25, 0.3, set(), split=True)
    assert [(s["keep"], s["why"]) for s in S] == [(True, "primary"), (True, "isoform_classifier"),
                                                  (True, "isoform_classifier")]
    assert {s["sub"] for s in S} == {0}
    cg.classifier_select(S, {0: S[:]}, [0.9, 0.5, 0.28], 0.25, 0.3, set(), split=True)
    assert (S[2]["keep"], S[2]["why"]) == (False, "isoform_rejected")
    cg.classifier_select(S, {0: S[:]}, [0.9, 0.5, 0.8], 0.25, 0.3, {2}, split=True)
    assert (S[2]["keep"], S[2]["why"]) == (False, "dropped_lone_single_exon")


def test_split_rescue_and_unselectable_members():
    # first gene below the threshold, rescued; the disjoint second gene is open; an unselectable one is dropped
    S = [_seg(0, 100, 300), _seg(1, 100, 250), _seg(2, 600, 900), _seg(3, 1000, 1200)]
    cg.classifier_select(S, {0: S[:]}, [0.2, 0.1, 0.7, 0.9], 0.25, 0.3, {3}, rescuable={1}, split=True)
    assert (S[2]["why"], S[2]["sub"]) == ("primary", 0)                     # the best open candidate comes first
    assert (S[1]["keep"], S[1]["why"], S[1]["sub"]) == (True, "primary_protein_rescue", 1)
    assert (S[0]["keep"], S[0]["why"]) == (False, "isoform_rejected")
    assert (S[3]["keep"], S[3]["why"]) == (False, "dropped_lone_single_exon")
    # nothing open and nothing rescuable: the whole cluster is below the threshold
    cg.classifier_select(S, {0: S[:]}, [0.2, 0.1, 0.1, 0.9], 0.25, 0.3, {3}, split=True)
    assert [s["why"] for s in S] == ["locus_below_min"] * 3 + ["dropped_lone_single_exon"]


@pytest.fixture
def fusion(tmp_path: Path) -> Path:
    """A Vipsania model fused over two genes (A and D) puts A, the single-exon C between them and D into one
    cluster; D is also predicted on its own by Tiberius and Vipsania."""
    A = [(1001, 1300), (1401, 1600)]
    C = [(3001, 3600)]
    D = [(5001, 5300), (5401, 5900)]
    write_gtf(tmp_path / "tiberius.gtf", "Tiberius", [("tA", A, "+"), ("tC", C, "+"), ("tD", D, "+")])
    write_gtf(tmp_path / "vipsania.gtf", "Vipsania", [("vA", A, "+"), ("vD", D, "+"), ("vF", A + D, "+")])
    write_gtf(tmp_path / "drusilla.gtf", "Drusilla", [("dA", A, "+")])
    write_gtf(tmp_path / "td2.gtf", "TD2", [("xA", A, "+")])
    return tmp_path


def test_split_loci_end_to_end(fusion):
    base = ("--classifier", str(MODEL), "--t-locus", "0", "--t-iso", "1.01")
    rows, kept, err = combine(fusion, "nosplit", *base, "--split-loci", "none", "--settings", str(fusion / "nosplit.tsv"))
    assert len({r["cluster"] for r in rows.values()}) == 1 and kept == {"g1.t1"}
    assert "split_loci" not in read_settings(fusion / "nosplit.tsv") and "primary_split" not in err
    rows, kept, err = combine(fusion, "split", *base, "--split-loci", "span", "--settings", str(fusion / "split.tsv"))
    assert kept == {"g1.t1", "g2.t1", "g3.t1"}
    # one of the three genes is the cluster's first gene (the most probable), the others are split off
    assert sorted(rows[k]["reason"] for k in ("tA", "tC", "tD")) == ["primary", "primary_split", "primary_split"]
    assert rows["vF"]["reason"] == "isoform_rejected"
    s = read_settings(fusion / "split.tsv")
    assert (s["split_loci"], s["split_genes"]) == ("span", "2")
    assert "split_loci=span" in err and "primary_split: 2" in err
    # genes in genome order: A, C, D
    genes = [line.split("\t")[3] for line in (fusion / "split.gtf").read_text().splitlines()
             if line.split("\t")[2] == "gene"]
    assert genes == ["1001", "3001", "5001"]


def test_cds_overlap_needs_shared_coding_bases():
    a = {"segs": ((100, 200), (500, 600))}
    nested = {"segs": ((250, 400),)}                 # inside the intron of a: spans overlap, no shared base
    assert cg.spans_overlap(a, nested) and not cg.cds_overlap(a, nested)
    assert cg.cds_overlap(a, {"segs": ((150, 160),)}) and cg.cds_overlap(a, {"segs": ((300, 500),)})
    assert cg.cds_overlap(a, {"segs": ((600, 700),)}) and not cg.cds_overlap(a, {"segs": ((601, 700),)})
    assert cg.cds_overlap(a, {"segs": ((50, 99), (210, 260), (590, 610))})


def test_split_cds_separates_a_gene_nested_in_an_intron():
    # primary 100..900 with an intron 301..599; 1 lies in the intron; 2 shares the first exon (an isoform)
    prim = {"_ix": 0, "tools": {"tiberius"}, "segs": ((100, 300), (600, 900)), "rank": 0, "keep": False, "why": ""}
    nested = {"_ix": 1, "tools": {"tiberius"}, "segs": ((350, 450), (500, 550)), "rank": 0, "keep": False, "why": ""}
    iso = {"_ix": 2, "tools": {"vipsania"}, "segs": ((150, 300), (650, 900)), "rank": 0, "keep": False, "why": ""}
    S = [prim, nested, iso]
    cg.classifier_select(S, {0: S[:]}, [0.9, 0.6, 0.5], 0.25, 0.3, set(), split="span")
    assert [(s["keep"], s["why"]) for s in S] == [(True, "primary"), (True, "isoform_classifier"),
                                                  (True, "isoform_classifier")]
    cg.classifier_select(S, {0: S[:]}, [0.9, 0.6, 0.5], 0.25, 0.3, set(), split="cds")
    assert [(s["keep"], s["why"]) for s in S] == [(True, "primary"), (True, "primary_split"),
                                                  (True, "isoform_classifier")]
    assert (prim["sub"], nested["sub"], iso["sub"]) == (0, 1, 0)
    # split=True is the span criterion
    cg.classifier_select(S, {0: S[:]}, [0.9, 0.6, 0.5], 0.25, 0.3, set(), split=True)
    assert nested["why"] == "isoform_classifier"


def test_split_loci_cds_end_to_end(tmp_path):
    P = [(1001, 1300), (2001, 2300)]                 # gene with a long intron
    N = [(1501, 1650), (1751, 1900)]                 # gene nested in that intron
    write_gtf(tmp_path / "tiberius.gtf", "Tiberius", [("tP", P, "+"), ("tN", N, "+")])
    write_gtf(tmp_path / "vipsania.gtf", "Vipsania", [("vP", P, "+"), ("vN", N, "+")])
    write_gtf(tmp_path / "drusilla.gtf", "Drusilla", [("dP", P, "+")])
    write_gtf(tmp_path / "td2.gtf", "TD2", [("xP", P, "+")])
    base = ("--classifier", str(MODEL), "--t-locus", "0", "--t-iso", "0")
    _, kept, _ = combine(tmp_path, "span", *base, "--split-loci", "span")
    assert kept == {"g1.t1", "g1.t2"}                # the nested gene is an isoform: spans overlap
    rows, kept, err = combine(tmp_path, "cds", *base, "--split-loci", "cds", "--settings", str(tmp_path / "s.tsv"))
    assert kept == {"g1.t1", "g2.t1"} and "split_loci=cds" in err
    assert read_settings(tmp_path / "s.tsv")["split_loci"] == "cds"
    r = combine(tmp_path, "bad", *base, "--split-loci", "exon", check=False)
    assert r.returncode != 0 and "invalid choice" in r.stderr


def test_split_loci_default_is_cds_and_none_switches_it_off(fusion):
    base = ("--classifier", str(MODEL), "--t-locus", "0", "--t-iso", "1.01")
    outs = {}
    for name, opts in (("default", ()), ("bare", ("--split-loci",)), ("cds", ("--split-loci", "cds")),
                       ("none", ("--split-loci", "none"))):
        _, kept, err = combine(fusion, name, *base, *opts, "--settings", str(fusion / f"{name}.tsv"))
        outs[name] = (fusion / f"{name}.gtf").read_text()
        s = read_settings(fusion / f"{name}.tsv")
        if name == "none":
            assert kept == {"g1.t1"} and "split_loci" not in s and "split_loci=" not in err
        else:
            assert len(kept) == 3 and s["split_loci"] == "cds" and "split_loci=cds" in err
    assert outs["default"] == outs["bare"] == outs["cds"] != outs["none"]
    # the tier rules: no split, and the default does not need the classifier
    _, _, err = combine(fusion, "rules")
    assert "split_loci" not in err and "primary_split" not in err


def test_split_loci_needs_the_classifier(fusion):
    r = combine(fusion, "rules", "--split-loci", check=False)
    assert r.returncode != 0 and "--split-loci needs --classifier" in r.stderr


def test_split_loci_none_reproduces_the_default_on_input_without_split_clusters(synth):
    combine(synth, "a", "--classifier", str(MODEL), "--split-loci", "none")
    combine(synth, "b", "--classifier", str(MODEL))
    # no cluster of the synthetic input holds disjoint candidates: same output
    assert (synth / "a.gtf").read_text() == (synth / "b.gtf").read_text()


def test_base_config_default_of_the_locus_split():
    text = (ROOT / "conf" / "base.config").read_text()
    assert re.search(r"^\s*split_loci\s*=\s*'cds'\s*$", text, re.M)


def test_base_config_defaults_of_the_protein_rescue():
    text = (ROOT / "conf" / "base.config").read_text()
    assert re.search(r"^\s*rescue_protein_cov\s*=\s*0\.9\s*$", text, re.M)
    assert re.search(r"^\s*rescue_protein_introns\s*=\s*0\.5\s*$", text, re.M)
