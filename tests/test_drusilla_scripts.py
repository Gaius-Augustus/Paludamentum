"""Unit tests of the Drusilla flow scripts in bin/ on toy inputs.

fix_stop_by_miniprot.py and compute_orf_features.py need pyfaidx,
apply_lgb_model_gtf.py pandas and lightgbm; the tests that need them are
skipped without.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import shutil
import subprocess
import tarfile
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "bin"


def load(name: str):
    spec = importlib.util.spec_from_file_location(name, BIN / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module   # dataclasses look up their module
    spec.loader.exec_module(module)
    return module


def run(script: str, *args: str, cwd: Path) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(BIN / script), *args], capture_output=True, text=True, cwd=cwd)


def rows(text: str) -> list[list[str]]:
    return [line.split("\t") for line in text.splitlines() if line and not line.startswith("#")]


# ---------------------------------------------------------------- gff_to_cds_gtf.py

def test_gff3_to_cds_gtf(tmp_path: Path):
    gff = tmp_path / "a.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "c1\tTiberius\tgene\t1\t300\t.\t+\t.\tID=g1\n"
        "c1\tTiberius\tmRNA\t1\t300\t.\t+\t.\tID=g1.t1;Parent=g1\n"
        "c1\tTiberius\texon\t1\t300\t.\t+\t.\tID=g1.t1.e1;Parent=g1.t1\n"
        "c1\tTiberius\tCDS\t1\t100\t.\t+\t0\tID=g1.t1.c1;Parent=g1.t1\n"
        "c1\tTiberius\tCDS\t201\t300\t.\t+\t2\tID=g1.t1.c2;Parent=g1.t1\n"
        "c1\tStringTie\ttranscript\t500\t600\t.\t-\t.\tID=nc1;Parent=ncg\n"      # no CDS: not written
    )
    proc = run("gff_to_cds_gtf.py", str(gff), cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    out = rows(proc.stdout)
    assert [(r[2], r[3], r[4], r[7], r[8]) for r in out] == [
        ("transcript", "1", "300", ".", 'gene_id "g1"; transcript_id "g1.t1";'),
        ("CDS", "1", "100", "0", 'gene_id "g1"; transcript_id "g1.t1";'),
        ("CDS", "201", "300", "2", 'gene_id "g1"; transcript_id "g1.t1";'),
    ]


def test_gtf_input_keeps_ids(tmp_path: Path):
    gtf = tmp_path / "a.gtf"
    gtf.write_text('c1\tX\tCDS\t1\t99\t.\t+\t0\tgene_id "g"; transcript_id "g.t1";\n')
    out = rows(run("gff_to_cds_gtf.py", str(gtf), cwd=tmp_path).stdout)
    assert out == [["c1", "X", "CDS", "1", "99", ".", "+", "0", 'gene_id "g"; transcript_id "g.t1";']]


# ---------------------------------------------------------------- fix_stop_by_miniprot.py

def genome_with(internal_stop: bool) -> str:
    """0-based: ORF [0,300) without stop; miniprot CDS [0,450); TAA at [450,453)."""
    codons = ["ATG"] + ["GCT"] * 99            # [0,300)
    ext = ["GCT"] * 50                         # [300,450)
    if internal_stop:
        ext[10] = "TAG"                        # in frame at [330,333)
    return "".join(codons + ext) + "TAA" + "GCT" * 20


def fix_stop_inputs(tmp_path: Path, internal_stop: bool) -> Path:
    (tmp_path / "g.fa").write_text(">c1\n" + genome_with(internal_stop) + "\n")
    (tmp_path / "orfs.gtf").write_text("")
    (tmp_path / "partial.gtf").write_text('c1\tDrusilla\tCDS\t1\t300\t.\t+\t0\ttranscript_id "t1"; gene_id "t1";\n')
    (tmp_path / "mp.gff").write_text(
        "c1\tminiprot\tmRNA\t1\t450\t90\t+\t.\tID=MP1;Identity=0.9;StopCodon=0\n"
        "c1\tminiprot\tCDS\t1\t450\t90\t+\t0\tParent=MP1\n"
    )
    return tmp_path / "out.gtf"


@pytest.mark.parametrize("internal_stop", [False, True])
def test_extension_through_an_in_frame_stop_is_dropped(tmp_path: Path, internal_stop: bool):
    pytest.importorskip("pyfaidx")
    out = fix_stop_inputs(tmp_path, internal_stop)
    proc = run("fix_stop_by_miniprot.py", "--orfs", "orfs.gtf", "--partial", "partial.gtf",
               "--miniprot", "mp.gff", "--genome", "g.fa", "--out", str(out), cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    lines = rows(out.read_text())
    if internal_stop:
        assert lines == [], "a CDS with an in-frame TAG must not be written"
        assert "in-frame stop codon: 1" in proc.stderr
    else:
        assert [(r[2], r[3], r[4]) for r in lines] == [("CDS", "1", "453")]
        assert "in-frame stop codon: 0" in proc.stderr


def test_has_internal_stop_reads_the_minus_strand(tmp_path: Path):
    pyfaidx = pytest.importorskip("pyfaidx")
    fs = load("fix_stop_by_miniprot")
    # minus strand CDS [0,9): the reverse complement of c1 is ATGTAGTAA, of c2 ATGGCTTAA
    (tmp_path / "m.fa").write_text(">c1\nTTACTACAT\n>c2\nTTAAGCCAT\n")
    genome = pyfaidx.Fasta(str(tmp_path / "m.fa"), as_raw=True, sequence_always_upper=True)
    assert fs.cds_sequence([(0, 9)], "c1", "-", genome) == "ATGTAGTAA"
    assert fs.has_internal_stop([(0, 9)], "c1", "-", genome)
    assert not fs.has_internal_stop([(0, 9)], "c2", "-", genome)


@pytest.mark.parametrize("intron, written", [
    ("GT" + "A" * 96 + "AG", True),
    ("C" * 100, False),                        # neither GT nor AG: no CDS across it
])
def test_extension_introns_need_canonical_splice_sites(tmp_path: Path, intron: str, written: bool):
    """Partial ORF [0,300) joined to the miniprot exon [400,550): new intron [300,400)."""
    pytest.importorskip("pyfaidx")
    genome = "ATG" + "GCT" * 99 + intron + "GCT" * 50 + "TAA" + "GCT" * 20
    (tmp_path / "g.fa").write_text(">c1\n" + genome + "\n")
    (tmp_path / "orfs.gtf").write_text("")
    (tmp_path / "partial.gtf").write_text('c1\tDrusilla\tCDS\t1\t300\t.\t+\t0\ttranscript_id "t1"; gene_id "t1";\n')
    (tmp_path / "mp.gff").write_text(
        "c1\tminiprot\tmRNA\t1\t550\t90\t+\t.\tID=MP1;Identity=0.9;StopCodon=0\n"
        "c1\tminiprot\tCDS\t1\t300\t90\t+\t0\tParent=MP1\n"
        "c1\tminiprot\tCDS\t401\t550\t90\t+\t0\tParent=MP1\n"
    )
    out = tmp_path / "out.gtf"
    proc = run("fix_stop_by_miniprot.py", "--orfs", "orfs.gtf", "--partial", "partial.gtf",
               "--miniprot", "mp.gff", "--genome", "g.fa", "--out", str(out), cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    got = [(r[3], r[4]) for r in rows(out.read_text())]
    if written:
        assert got == [("1", "300"), ("401", "553")]
        assert "non-canonical new intron: 0" in proc.stderr
    else:
        assert got == []
        assert "non-canonical new intron: 1" in proc.stderr


def _rc_segs(segs, length):
    return [(length - e, length - s) for s, e in segs]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_stop_fix_splices_a_retained_intron_at_the_miniprot_donor(tmp_path: Path, strand: str):
    """+ strand: exon [0,150), intron [150,250) with TAA in frame at 177, exon [250,400), TAA.

    The complete ORF [0,180) reads into the intron to the TAA. The fix must use
    the miniprot intron, not one that starts at the spurious stop (TA...AG).
    """
    pyfaidx = pytest.importorskip("pyfaidx")
    fs = load("fix_stop_by_miniprot")
    intron = "GT" + "C" * 25 + "TAA" + "C" * 68 + "AG"
    plus = "ATG" + "GCT" * 49 + intron + "GCT" * 50 + "TAA" + "GCT" * 20
    seq, orf_segs, mp_segs, want = plus, [(0, 180)], [(0, 150), (250, 400)], [(0, 150), (250, 403)]
    if strand == "-":
        n = len(plus)
        seq = fs._rev_comp(plus)
        orf_segs, mp_segs, want = (sorted(_rc_segs(x, n)) for x in (orf_segs, mp_segs, want))
    (tmp_path / "g.fa").write_text(">c1\n" + seq + "\n")
    genome = pyfaidx.Fasta(str(tmp_path / "g.fa"), as_raw=True, sequence_always_upper=True)
    orf = fs.ORF("t1", "c1", strand, False, "Drusilla")
    orf.segments = orf_segs
    mp = fs.MpAlignment("MP1", "c1", strand, 0.9, False)
    mp.cds_segments = mp_segs
    got, via_hints = fs.try_fix(orf, [mp], genome, {"c1": len(seq)}, 5000, 30)
    assert got == want and not via_hints
    assert not fs.has_internal_stop(got, "c1", strand, genome)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_start_fix_begins_with_the_atg_upstream_of_the_alignment(tmp_path: Path, strand: str):
    """+ strand: ATG at 90, miniprot CDS starts at 96, ORF [150,300) ends with TAA."""
    pyfaidx = pytest.importorskip("pyfaidx")
    fs = load("fix_stop_by_miniprot")
    plus = "C" * 90 + "ATG" + "GCT" * 68 + "TAA" + "C" * 30
    seq, orf_segs, mp_segs, want = plus, [(150, 300)], [(96, 150)], [(90, 300)]
    if strand == "-":
        n = len(plus)
        seq = fs._rev_comp(plus)
        orf_segs, mp_segs, want = (_rc_segs(x, n) for x in (orf_segs, mp_segs, want))
    (tmp_path / "g.fa").write_text(">c1\n" + seq + "\n")
    genome = pyfaidx.Fasta(str(tmp_path / "g.fa"), as_raw=True, sequence_always_upper=True)
    orf = fs.ORF("t1", "c1", strand, False, "Drusilla")
    orf.segments = orf_segs
    mp = fs.MpAlignment("MP1", "c1", strand, 0.9, False)
    mp.cds_segments = mp_segs
    got = fs.try_fix_5prime(orf, [mp], genome, {"c1": len(seq)}, 5000, 30)
    assert got == want
    assert fs.cds_sequence(got, "c1", strand, genome).startswith("ATG")


def test_start_codon_in_a_miniprot_intron_is_not_used():
    pytest.importorskip("pyfaidx")
    fs = load("fix_stop_by_miniprot")
    mp = [(100, 200), (300, 400)]
    assert fs._upstream_ext(mp, "+", 250, 450, 600) is None
    assert fs._upstream_ext(mp, "+", 150, 450, 600) == [(150, 200), (300, 400)]
    assert fs._upstream_ext(mp, "-", 250, 0, 50) is None
    assert fs._upstream_ext(mp, "-", 150, 0, 50) == [(100, 153)]


# ---------------------------------------------------------------- apply_lgb_model_gtf.py

# feature_names of the released model drusilla_lgb_3class_v1, in model order
RELEASED_FEATURES = [
    "n_exons", "cds_length_nt", "dist_upstream_stop_nt", "n_upstream_atgs",
    "n_overlapping_alignments", "best_identity", "best_norm_bitscore",
    "best_target_coverage", "frac_introns_supported", "cds_length_pct",
    "protein_extends_5prime_codons", "protein_extends_3prime_codons",
    "has_protein_support", "has_start_hint", "has_stop_hint", "has_conflict",
    "has_upstream_partner", "has_downstream_partner",
    "lorf_class__LORF_UPSTOP", "lorf_class__LORF_NOUPSTOP", "lorf_class__sORF_UPSTOP",
    "lorf_class__sORF_NOUPSTOP", "lorf_class__upLORF",
    "support_level__fullSupport", "support_level__anySupport", "support_level__noSupport",
]


def load_lgb():
    pytest.importorskip("pandas")
    pytest.importorskip("lightgbm")
    return load("apply_lgb_model_gtf")


def test_feature_matrix_follows_the_model_columns():
    pd = pytest.importorskip("pandas")
    lgb = load_lgb()
    df = pd.DataFrame({
        "transcript_id": ["a", "b"],
        "n_exons": [1, 3], "cds_length_nt": [300, 900], "has_protein_support": [1, 0],
        "lorf_class": ["upLORF", "LORF_UPSTOP"], "support_level": ["fullSupport", "noSupport"],
    })
    names = ["cds_length_nt", "n_exons", "lorf_class__upLORF", "has_protein_support"]
    X = lgb.build_feature_matrix(df, names)
    assert X.shape == (2, 4)
    assert list(X[0]) == [300, 1, 1, 1]      # model order
    assert list(X[1]) == [900, 3, 0, 0]


def test_a_missing_feature_is_an_error():
    pd = pytest.importorskip("pandas")
    lgb = load_lgb()
    df = pd.DataFrame({"transcript_id": ["a"], "n_exons": [1]})
    with pytest.raises(SystemExit, match="best_identity"):
        lgb.build_feature_matrix(df, ["n_exons", "best_identity"])


def test_feature_table_has_every_feature_of_the_released_model():
    """The columns of compute_orf_features.py give all features of the model."""
    pd = pytest.importorskip("pandas")
    pytest.importorskip("pyfaidx")
    lgb = load_lgb()
    features = load("compute_orf_features")
    row = {c: 0 for c in features.COLUMNS}
    row.update(transcript_id="t1", lorf_class="NA", support_level="noSupport")
    X = lgb.build_feature_matrix(pd.DataFrame([row]), RELEASED_FEATURES)
    assert X.shape == (1, len(RELEASED_FEATURES))


def test_empty_feature_table_is_not_scored(tmp_path: Path):
    lgb = load_lgb()
    tsv = tmp_path / "features.tsv"
    tsv.write_text("transcript_id\tn_exons\n")

    class NeverCalled:
        def predict(self, X):
            raise AssertionError("predict on an empty table")

    passing, scores, attrs = lgb.score_features(tsv, NeverCalled(), ["n_exons"], 0.5, "prob_not_wrong")
    assert passing == set() and attrs == {} and len(scores) == 0


def write_model(tmp_path: Path) -> tuple[Path, Path, str]:
    """A tiny 3-class text model on two features, its .json and the release
    tar.gz of both. Returns (directory, archive, sha256 of the archive)."""
    lightgbm = pytest.importorskip("lightgbm")
    import numpy as np
    rng = np.random.default_rng(1)
    X = rng.normal(size=(300, 2))
    y = (X[:, 0] > -0.5).astype(int) + (X[:, 0] > 0.5).astype(int)   # class by the first feature
    booster = lightgbm.train({"objective": "multiclass", "num_class": 3, "verbose": -1,
                              "min_data_in_leaf": 5}, lightgbm.Dataset(X, y), num_boost_round=10)
    d = tmp_path / "m"
    d.mkdir()
    txt = d / "m.txt"
    booster.save_model(str(txt))
    meta = {"name": "m", "model_file": "m.txt",
            "model_sha256": hashlib.sha256(txt.read_bytes()).hexdigest(),
            "feature_names": ["n_exons", "cds_length_nt"],
            "classes": {"0": "wrong", "1": "partial", "2": "correct"}}
    (d / "m.json").write_text(json.dumps(meta))
    tgz = tmp_path / "m.tar.gz"
    with tarfile.open(tgz, "w:gz") as tar:
        tar.add(txt, arcname="m/m.txt")
        tar.add(d / "m.json", arcname="m/m.json")
    return d, tgz, hashlib.sha256(tgz.read_bytes()).hexdigest()


def test_model_loads_from_archive_directory_and_txt(tmp_path: Path):
    lgb = load_lgb()
    d, tgz, sha = write_model(tmp_path)
    for path, checksum in [(tgz, sha), (tgz, None), (d, None), (d / "m.txt", None)]:
        booster, names, meta = lgb.load_model(path, checksum)
        assert names == ["n_exons", "cds_length_nt"] and meta["name"] == "m"
        assert booster.predict([[2.0, 0.0]]).argmax() == 2


def test_model_checksums_are_checked(tmp_path: Path):
    lgb = load_lgb()
    d, tgz, sha = write_model(tmp_path)
    with pytest.raises(SystemExit, match="sha256"):
        lgb.load_model(tgz, "0" * 64)
    # a text model that is not the one of the .json
    (d / "m.txt").write_text((d / "m.txt").read_text().replace("num_class=3", "num_class=3 "))
    with pytest.raises(SystemExit, match="model_sha256"):
        lgb.load_model(d, None)


def test_pickle_is_refused(tmp_path: Path):
    lgb = load_lgb()
    pkl = tmp_path / "lgb_3class_model.pkl"
    pkl.write_bytes(b"not loaded")
    with pytest.raises(SystemExit, match="pickle"):
        lgb.load_model(pkl, None)


def test_filter_scores_and_marks_the_gtf(tmp_path: Path):
    pytest.importorskip("pandas")
    _, tgz, sha = write_model(tmp_path)
    (tmp_path / "features.tsv").write_text(
        "transcript_id\tn_exons\tcds_length_nt\n"
        "good\t2.0\t0\n"
        "bad\t-2.0\t0\n"
    )
    (tmp_path / "in.gtf").write_text(
        'c1\tTiberius\tCDS\t1\t300\t.\t+\t0\ttranscript_id "good"; gene_id "g1";\n'
        'c1\tTiberius\tCDS\t501\t800\t.\t+\t0\ttranscript_id "bad"; gene_id "g2";\n'
    )
    proc = run("apply_lgb_model_gtf.py", "--model", str(tgz), "--sha256", sha,
               "--features", "features.tsv", "--in-gtf", "in.gtf", "--out-gtf", "out.gtf",
               cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    out = rows((tmp_path / "out.gtf").read_text())
    assert len(out) == 1 and 'transcript_id "good"' in out[0][8]
    assert 'lgb_class "correct"' in out[0][8]
    scores = (tmp_path / "out.scores.tsv").read_text().splitlines()
    assert len(scores) == 3 and scores[0].startswith("transcript_id\tprob_wrong")


# ---------------------------------------------------------------- prepare_hint_rescue_loci.py

def gtf_line(feature: str, start: int, end: int, tx: str) -> str:
    return f'chr1\tt\t{feature}\t{start}\t{end}\t.\t+\t0\tgene_id "{tx}.g"; transcript_id "{tx}";\n'


def rescue_loci(tmp_path: Path, orf_cds: list[tuple[int, int]], *extra: str) -> list[str]:
    """Locus ids of one partial transcript with a two-intron chain and one ORF."""
    (tmp_path / "genome.fa").write_text(">chr1\n" + "ACGT" * 500 + "\n")
    subprocess.run(["samtools", "faidx", "genome.fa"], cwd=tmp_path, check=True)
    (tmp_path / "partial.gtf").write_text(gtf_line("transcript", 500, 1500, "p1"))
    (tmp_path / "correct.gtf").write_text("")
    (tmp_path / "chained.gff").write_text("".join(
        f"chr1\tc\tintron\t{s}\t{e}\t.\t+\t.\tchain_id=c1;al_score=10\n"
        for s, e in [(700, 800), (1000, 1100)]))
    (tmp_path / "orfs.gtf").write_text("".join(gtf_line("CDS", s, e, "o1") for s, e in orf_cds))
    res = run("prepare_hint_rescue_loci.py", "--partial_gtf", "partial.gtf", "--correct_gtf", "correct.gtf",
              "--chained_hints", "chained.gff", "--orfs_gtf", "orfs.gtf", "--genome", "genome.fa",
              "--outdir", "rescue", "--flank", "100", *extra, cwd=tmp_path)
    assert res.returncode == 0, res.stderr
    manifest = (tmp_path / "rescue" / "loci_manifest.tsv").read_text().splitlines()[1:]
    return [line.split("\t")[0] for line in manifest]


AGREEING_ORF = [(500, 699), (801, 999), (1101, 1400)]   # has both chain introns


@pytest.mark.skipif(shutil.which("samtools") is None, reason="needs samtools")
def test_orf_filter_is_off_by_default(tmp_path: Path):
    assert rescue_loci(tmp_path, AGREEING_ORF) == ["locus_0000000"]


@pytest.mark.skipif(shutil.which("samtools") is None, reason="needs samtools")
def test_orf_filter_skips_loci_whose_chain_introns_an_orf_has(tmp_path: Path):
    assert rescue_loci(tmp_path, AGREEING_ORF, "--orf_filter") == []


@pytest.mark.skipif(shutil.which("samtools") is None, reason="needs samtools")
def test_orf_filter_keeps_loci_with_a_chain_intron_the_orfs_lack(tmp_path: Path):
    assert rescue_loci(tmp_path, [(500, 699), (801, 1400)], "--orf_filter") == ["locus_0000000"]
