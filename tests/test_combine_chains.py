"""bin/combine_genesets.py --chains: full-length transcripts (StringTie -L of long reads) as evidence.
Unit tests of read_chains / chain_support / adjust_tier / rank_score and a run of the script
on a tiny synthetic set of four tool GTFs and a chains GTF."""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "bin"
SCRIPT = BIN / "combine_genesets.py"

sys.path.insert(0, str(BIN))
import combine_genesets as cg  # noqa: E402


# ---------------------------------------------------------------- helpers

def write_gtf(path: Path, src: str, recs: list) -> Path:
    """recs: (transcript id, contig, strand, feature, [(start, end), ...])."""
    with open(path, "w") as fh:
        for tid, contig, strand, feat, segs in recs:
            for a, b in segs:
                fh.write(f'{contig}\t{src}\t{feat}\t{a}\t{b}\t.\t{strand}\t0\t'
                         f'gene_id "{tid}"; transcript_id "{tid}";\n')
    return path


def chains_of(tmp_path: Path, recs: list):
    return cg.read_chains([str(write_gtf(tmp_path / "chains.gtf", "StringTie", recs))])


def model(segs, strand="+", contig="chr1") -> dict:
    return {"contig": contig, "strand": strand, "segs": tuple(segs)}


def tx(tid, exons, strand="+", contig="chr1"):
    return (tid, contig, strand, "exon", exons)


# ---------------------------------------------------------------- chain_support

@pytest.fixture
def ch(tmp_path):
    return chains_of(tmp_path, [tx("T1", [(50, 200), (300, 400), (500, 700)])])


def test_multi_exon_full(ch):
    assert cg.chain_support(model([(100, 200), (300, 400), (500, 600)]), ch) == ("full", 1)


def test_multi_exon_full_subrun_of_introns(ch):
    assert cg.chain_support(model([(100, 200), (300, 400)]), ch) == ("full", 1)
    assert cg.chain_support(model([(300, 400), (500, 600)]), ch) == ("full", 1)


def test_multi_exon_partial_when_only_one_intron_matches(ch):
    # first intron (201..299) matches, the second one differs
    assert cg.chain_support(model([(100, 200), (300, 400), (550, 600)]), ch) == ("partial", 0)


def test_intron_spanning_two_transcript_introns_is_none(ch):
    # intron 201..499 is not an intron of the transcript
    assert cg.chain_support(model([(100, 200), (500, 600)]), ch) == ("none", 0)


def test_multi_exon_none(ch):
    assert cg.chain_support(model([(100, 150), (250, 400)]), ch) == ("none", 0)


@pytest.mark.parametrize("kw", [dict(strand="-"), dict(contig="chr2")])
def test_other_strand_or_contig_is_none(ch, kw):
    assert cg.chain_support(model([(100, 200), (300, 400)], **kw), ch) == ("none", 0)


def test_first_segment_starting_before_the_flanking_exon_is_not_full(ch):
    assert cg.chain_support(model([(40, 200), (300, 400)]), ch) == ("partial", 0)


def test_last_segment_ending_after_the_flanking_exon_is_not_full(ch):
    assert cg.chain_support(model([(300, 400), (500, 710)]), ch) == ("partial", 0)


def test_terminal_segment_exactly_at_exon_boundary_is_full(ch):
    assert cg.chain_support(model([(50, 200), (300, 400), (500, 700)]), ch) == ("full", 1)


def test_single_exon_full_partial_none(tmp_path):
    c = chains_of(tmp_path, [tx("S", [(1000, 1500)])])
    assert cg.chain_support(model([(1100, 1400)]), c) == ("full", 1)
    assert cg.chain_support(model([(1000, 1500)]), c) == ("full", 1)
    assert cg.chain_support(model([(900, 1200)]), c) == ("partial", 0)    # overlaps, not inside
    assert cg.chain_support(model([(1100, 1600)]), c) == ("partial", 0)
    assert cg.chain_support(model([(2000, 2100)]), c) == ("none", 0)
    assert cg.chain_support(model([(1100, 1400)], strand="-"), c) == ("none", 0)
    assert cg.chain_support(model([(1100, 1400)], contig="chrX"), c) == ("none", 0)


def test_single_exon_inside_exon_of_multi_exon_transcript(ch):
    assert cg.chain_support(model([(310, 390)]), ch) == ("full", 1)


def test_unstranded_transcript_matches_both_strands(tmp_path):
    c = chains_of(tmp_path, [tx("U", [(50, 200), (300, 400)], strand="."),
                             tx("V", [(1000, 1500)], strand=".")])
    for st in "+-":
        assert cg.chain_support(model([(100, 200), (300, 380)], strand=st), c) == ("full", 1)
        assert cg.chain_support(model([(1100, 1400)], strand=st), c) == ("full", 1)


def test_chain_n_counts_transcripts_with_the_full_chain(tmp_path):
    c = chains_of(tmp_path, [tx("A", [(50, 200), (300, 400), (500, 700)]),
                             tx("B", [(80, 200), (300, 450)]),
                             tx("C", [(50, 200), (310, 400)])])
    assert cg.chain_support(model([(100, 200), (300, 400)]), c) == ("full", 2)
    assert cg.chain_support(model([(100, 200), (310, 400)]), c) == ("full", 1)


def test_empty_chains(tmp_path):
    c = chains_of(tmp_path, [])
    assert cg.chain_support(model([(100, 200), (300, 400)]), c) == ("none", 0)
    assert cg.chain_support(model([(100, 200)]), c) == ("none", 0)


# ---------------------------------------------------------------- tier / rank

def lone_tib(tier, support):
    return {"tools": {"tiberius"}, "tier": tier, "n_introns": 1, "tool_intron_frac": 0.0,
            "rna_frac": 0.0, "rna_cov": None, "chain_support": support}


def test_adjust_tier_up_chain_support():
    cfg = dict(chain_rules=True, agree_rules=False, rna_rules=False, have_rna=False, cov_full=0.9)
    s = lone_tib(3, "full")
    cg.adjust_tier(s, cfg)
    assert s["adj"] == "up_chain_support" and s["tier"] == cg.RESCUE_TIER and s["tier_base"] == 3


@pytest.mark.parametrize("support,rules", [("partial", True), ("none", True), ("", True), ("full", False)])
def test_adjust_tier_no_change_without_full_support_or_rules(support, rules):
    cfg = dict(chain_rules=rules, agree_rules=False, rna_rules=False, have_rna=False, cov_full=0.9)
    s = lone_tib(3, support)
    cg.adjust_tier(s, cfg)
    assert s["tier"] == 3 and not s.get("adj")


def test_rank_score_chain_bonus_only_for_full():
    base = {"tier": 5, "tools": {"tiberius"}, "cds_len": 600, "rna_frac": "", "n_introns": 1}
    full, part = dict(base, chain_support="full"), dict(base, chain_support="partial")
    assert cg.rank_score(full, 0.0, 1.0) == pytest.approx(cg.rank_score(part, 0.0, 1.0) + 1.0)
    assert cg.rank_score(full, 0.0, 0.0) == pytest.approx(cg.rank_score(part, 0.0, 0.0))


# ---------------------------------------------------------------- end to end

P = [(100, 200), (300, 400), (500, 600)]
I = [(100, 200), (500, 603)]               # 2 exons, chain lies in L.1.2
X = [(100, 200), (300, 400)]
Y = [(100, 200), (300, 400), (550, 600)]   # second intron differs: partial only
Z = [(1000, 1500)]                         # single exon, lone Tiberius, 501 nt, minus strand
W = [(100, 400)]


@pytest.fixture(scope="module")
def synth(tmp_path_factory):
    d = tmp_path_factory.mktemp("chains")
    write_gtf(d / "tiberius.gtf", "tib", [("P", "chr1", "+", "CDS", P), ("I", "chr1", "+", "CDS", I),
                                          ("Z", "chr1", "-", "CDS", Z), ("W", "chr2", "+", "CDS", W)])
    write_gtf(d / "vipsania.gtf", "vip", [("P", "chr1", "+", "CDS", P)])
    write_gtf(d / "td2.gtf", "td2", [("P", "chr1", "+", "CDS", P), ("X", "chr1", "+", "CDS", X),
                                     ("Y", "chr1", "+", "CDS", Y)])
    write_gtf(d / "drusilla.gtf", "dru", [("Y", "chr1", "+", "CDS", Y)])
    write_gtf(d / "longreads.gtf", "StringTie", [
        tx("L.1.1", [(50, 200), (300, 400), (500, 700)]),
        tx("L.1.2", [(50, 200), (500, 700)]),
        tx("L.2.1", [(900, 1600)], strand=".")])
    return d


def combine(synth: Path, name: str, *extra: str):
    out, rep = synth / f"{name}.gtf", synth / f"{name}.report.tsv"
    cmd = [sys.executable, str(SCRIPT), "--tiberius", str(synth / "tiberius.gtf"),
           "--vipsania", str(synth / "vipsania.gtf"), "--drusilla", str(synth / "drusilla.gtf"),
           "--td2", str(synth / "td2.gtf"), "-o", str(out), "--report", str(rep), *extra]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    with open(rep) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    by_id = {row["member_ids"].split(",")[0].split(":")[1]: row for row in rows}
    return out, by_id


def test_without_chains_the_columns_are_empty(synth):
    _, rows = combine(synth, "base")
    assert set(rows) == {"P", "I", "X", "Y", "Z", "W"}
    assert all(r["chain_support"] == "" and r["chain_n"] == "" for r in rows.values())
    assert rows["I"]["reason"] == "isoform_rejected"
    assert rows["Z"]["reason"] == "locus_below_min"


def test_with_chains_lone_tiberius_two_exon_model_becomes_isoform_chain(synth):
    _, rows = combine(synth, "chains", "--chains", str(synth / "longreads.gtf"))
    r = rows["I"]
    assert (r["chain_support"], r["chain_n"]) == ("full", "1")
    assert r["adjustment"] == "up_chain_support"
    assert float(r["locus_tier"]) == 6
    assert r["reason"] == "isoform_chain"
    assert rows["P"]["reason"] == "primary"
    assert rows["P"]["chain_support"] == "full"


def test_partial_drusilla_td2_model_and_iso_dt_chains(synth):
    _, rows = combine(synth, "chains2", "--chains", str(synth / "longreads.gtf"))
    assert rows["Y"]["chain_support"] == "partial"
    assert rows["Y"]["reason"] == "isoform_drusilla_td2"
    _, rows = combine(synth, "chains_dt", "--chains", str(synth / "longreads.gtf"), "--iso-dt-chains")
    assert rows["Y"]["reason"] == "isoform_rejected"
    assert rows["I"]["reason"] == "isoform_chain"


def test_single_exon_lone_tiberius_in_unstranded_transcript_opens_the_locus(synth):
    _, rows = combine(synth, "chains3", "--chains", str(synth / "longreads.gtf"))
    z = rows["Z"]
    assert z["chain_support"] == "full" and z["adjustment"] == "up_chain_support"
    assert z["reason"] == "primary"
    assert rows["W"]["chain_support"] == "none" and rows["W"]["reason"] == "locus_below_min"


def test_no_chain_rules_gives_the_same_gtf_as_without_chains(synth):
    base, _ = combine(synth, "base_a")
    off, rows = combine(synth, "off", "--chains", str(synth / "longreads.gtf"), "--no-chain-rules")
    assert base.read_text() == off.read_text()
    assert rows["I"]["reason"] == "isoform_rejected" and rows["Z"]["reason"] == "locus_below_min"
    assert rows["I"]["chain_support"] == "full"   # the support is still reported
