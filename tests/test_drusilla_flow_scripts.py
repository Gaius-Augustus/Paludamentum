"""Unit tests of the Drusilla flow scripts in bin/ on toy inputs: the StringTie
pre-filter, the ORF feature table of the LightGBM filter, and both steps of
the hint rescue.

compute_orf_features.py needs pyfaidx, prepare_hint_rescue_loci.py samtools;
the tests that need them are skipped without.
"""
from __future__ import annotations

import csv
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "bin"


def run(script: str, *args: str, cwd: Path) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(BIN / script), *args], capture_output=True, text=True, cwd=cwd)


def rows(text: str) -> list[list[str]]:
    return [line.split("\t") for line in text.splitlines() if line and not line.startswith("#")]


def table(path: Path) -> dict[str, dict[str, str]]:
    with open(path) as fh:
        return {r[next(iter(r))]: r for r in csv.DictReader(fh, delimiter="\t")}


# ---------------------------------------------------------------- filter_stringtie_gtf.py

def stringtie_tx(tx: str, exons: list[tuple[int, int]], cov: str | None, tpm: str | None) -> str:
    """A StringTie transcript with its exons; cov or TPM None leaves the attribute out."""
    attrs = f'gene_id "{tx}.g"; transcript_id "{tx}";'
    tattrs = attrs + (f' cov "{cov}";' if cov else "") + (f' TPM "{tpm}";' if tpm else "")
    lines = [f"c1\tStringTie\ttranscript\t{exons[0][0]}\t{exons[-1][1]}\t1000\t+\t.\t{tattrs}\n"]
    lines += [f'c1\tStringTie\texon\t{s}\t{e}\t1000\t+\t.\t{attrs} exon_number "{i}";\n'
              for i, (s, e) in enumerate(exons, 1)]
    return "".join(lines)


STRINGTIE = {
    # transcript: (exons, cov, TPM, reason with the default thresholds)
    "short":    ([(1, 100), (201, 300)], "10", "5", "length"),       # 200 nt
    "low_cov":  ([(1, 200), (301, 500)], "2.9", "5", "cov"),
    "no_cov":   ([(1, 400)], None, "5", "cov"),
    "low_tpm":  ([(1, 400)], "5", "0.99", "tpm"),
    "edge":     ([(1, 100), (201, 400)], "3", "1", "pass"),         # 300 nt, cov 3, TPM 1
    "long":     ([(1, 1000), (2001, 4000)], "4", "0.5", "pass_long_relaxed"),   # 3000 nt
    "long_low": ([(1, 3000)], "4", "0.49", "tpm_long"),
}


def write_stringtie(tmp_path: Path) -> Path:
    gtf = tmp_path / "stringtie.gtf"
    gtf.write_text("# StringTie version 2.2.1\n" + "".join(
        stringtie_tx(tx, exons, cov, tpm) for tx, (exons, cov, tpm, _) in STRINGTIE.items()))
    return gtf


def test_stringtie_filter_rules(tmp_path: Path):
    gtf = write_stringtie(tmp_path)
    proc = run("filter_stringtie_gtf.py", "--in-gtf", str(gtf), "--out-gtf", "out.gtf",
               "--out-tsv", "decisions.tsv", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    decisions = table(tmp_path / "decisions.tsv")
    assert {tx: d["reason"] for tx, d in decisions.items()} == {tx: v[3] for tx, v in STRINGTIE.items()}
    assert {tx for tx, d in decisions.items() if d["kept"] == "1"} == {"edge", "long"}
    assert decisions["long"]["length"] == "3000" and decisions["no_cov"]["cov"] == "None"
    # the kept transcripts with all their lines, in input order, and the header
    out = (tmp_path / "out.gtf").read_text()
    assert out == ("# StringTie version 2.2.1\n" + stringtie_tx("edge", *STRINGTIE["edge"][:3])
                   + stringtie_tx("long", *STRINGTIE["long"][:3]))
    assert "transcripts total=7 kept=2 (28.6%) dropped=5" in proc.stderr


def test_stringtie_filter_thresholds_are_options(tmp_path: Path):
    gtf = write_stringtie(tmp_path)
    proc = run("filter_stringtie_gtf.py", "--in-gtf", str(gtf), "--out-gtf", "out.gtf",
               "--min-length", "200", "--min-cov", "2", "--min-tpm", "0.9",
               "--long-length", "5000", "--min-tpm-long", "0", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    # the decision log defaults to <out-gtf>.decisions.tsv
    decisions = table(tmp_path / "out.gtf.decisions.tsv")
    kept = {tx for tx, d in decisions.items() if d["kept"] == "1"}
    assert kept == {"short", "low_cov", "low_tpm", "edge"}
    # no transcript is long now: both 3000 nt ones need TPM 0.9
    assert decisions["long"]["reason"] == decisions["long_low"]["reason"] == "tpm"


# ---------------------------------------------------------------- compute_orf_features.py

def orf_genome() -> str:
    """3000 nt of GCC repeats (no stop, no ATG in any frame on either strand) with:
    in frame 0 upstream of ORF a (ATG at 0-based 300) a TAA at 270 and an ATG at 285;
    upstream of ORF b (minus strand, ATG at 1497-1500) a TTA (TAA) at 1557."""
    seq = list("GCC" * 1000)
    for pos, codon in [(270, "TAA"), (285, "ATG"), (1557, "TTA")]:
        seq[pos:pos + 3] = codon
    return "".join(seq)


def orf_gtf_line(tx: str, strand: str, start: int, end: int, extra: str = "") -> str:
    return f'chr1\tDrusilla\tCDS\t{start}\t{end}\t.\t{strand}\t0\ttranscript_id "{tx}"; gene_id "{tx}.g";{extra}\n'


def write_orf_inputs(tmp_path: Path) -> list[str]:
    (tmp_path / "genome.fa").write_text(">chr1\n" + orf_genome() + "\n")
    (tmp_path / "proteins.faa").write_text(">prot1\n" + "M" * 120 + "\n>prot2\n" + "M" * 80 + "\n")
    # a: + strand, two exons [300,399) [500,650); b: - strand, one exon; c: + strand,
    # 1350 nt (in frame) downstream of a
    (tmp_path / "orfs.gtf").write_text(
        orf_gtf_line("a", "+", 301, 399) + orf_gtf_line("a", "+", 501, 650)
        + orf_gtf_line("b", "-", 1201, 1500, ' lorf_class "upLORF";')
        + orf_gtf_line("c", "+", 2001, 2300)
    )
    # MP1 covers a and extends it by 10 codons at the 5' end; MP2 starts 20 codons
    # downstream of MP1, with a lower identity
    (tmp_path / "miniprot.gff").write_text(
        "##gff-version 3\n"
        "chr1\tminiprot\tmRNA\t271\t650\t200\t+\t.\tID=MP1;Rank=1;Identity=0.9000;Target=prot1 1 100\n"
        "chr1\tminiprot\tCDS\t271\t399\t200\t+\t0\tParent=MP1;Rank=1;Identity=0.9000;Target=prot1 1 43\n"
        "chr1\tminiprot\tCDS\t501\t650\t200\t+\t0\tParent=MP1;Rank=1;Identity=0.9000;Target=prot1 44 100\n"
        "chr1\tminiprot\tmRNA\t331\t650\t90\t+\t.\tID=MP2;Rank=1;Identity=0.5000;Target=prot2 1 73\n"
        "chr1\tminiprot\tCDS\t331\t399\t90\t+\t0\tParent=MP2\n"
        "chr1\tminiprot\tCDS\t501\t650\t90\t+\t0\tParent=MP2\n"
    )
    # the intron, start and stop codon of a
    (tmp_path / "hc.gff").write_text(
        "chr1\tminiprothint\tintron\t400\t500\t.\t+\t.\tsrc=P\n"
        "chr1\tminiprothint\tstart_codon\t301\t303\t.\t+\t0\tsrc=P\n"
        "chr1\tminiprothint\tstop_codon\t648\t650\t.\t+\t0\tsrc=P\n"
    )
    (tmp_path / "ref.tmap").write_text(
        "ref_gene_id\tref_id\tclass_code\tqry_gene_id\tqry_id\n"
        "r1\tr1.t1\t=\ta.g\ta\n"
    )
    return ["--orfs-gtf", "orfs.gtf", "--genome", "genome.fa", "--out", "features.tsv"]


def test_orf_features(tmp_path: Path):
    pytest.importorskip("pyfaidx")
    args = write_orf_inputs(tmp_path)
    proc = run("compute_orf_features.py", *args, "--miniprot-gff", "miniprot.gff",
               "--proteins-fasta", "proteins.faa", "--hints-gff", "hc.gff",
               "--ref-tmap", "ref.tmap", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    f = table(tmp_path / "features.tsv")
    assert list(f) == ["a", "b", "c"]

    def cols(tx: str, *names: str) -> dict[str, str]:
        return {n: f[tx][n] for n in names}

    assert cols("a", "strand", "n_exons", "cds_length_nt", "lorf_class",
                "dist_upstream_stop_nt", "n_upstream_atgs", "gffcompare_class") == {
        "strand": "+", "n_exons": "2", "cds_length_nt": "249", "lorf_class": "NA",
        "dist_upstream_stop_nt": "30", "n_upstream_atgs": "1", "gffcompare_class": "=",
    }
    assert cols("a", "has_protein_support", "n_overlapping_alignments", "best_identity",
                "best_score", "best_norm_bitscore", "best_target_coverage", "best_protein_coverage",
                "protein_extends_5prime_codons", "protein_extends_3prime_codons",
                "has_conflict", "conflict_identity_delta") == {
        "has_protein_support": "1", "n_overlapping_alignments": "2", "best_identity": "0.9000",
        "best_score": "200.0000",
        "best_norm_bitscore": f"{200 / ((129 + 150) // 3):.4f}",   # score per aligned aa
        "best_target_coverage": "1.0000",
        "best_protein_coverage": f"{99 / 120:.4f}",               # Target 1-100 of 120 aa
        "protein_extends_5prime_codons": "10", "protein_extends_3prime_codons": "0",
        "has_conflict": "1", "conflict_identity_delta": "-0.4000",
    }
    assert cols("a", "n_introns_supported", "frac_introns_supported", "has_start_hint",
                "has_stop_hint", "support_level") == {
        "n_introns_supported": "1", "frac_introns_supported": "1.0000",
        "has_start_hint": "1", "has_stop_hint": "1", "support_level": "fullSupport",
    }
    # minus strand: upstream is to the right, read reverse complemented
    assert cols("b", "strand", "lorf_class", "dist_upstream_stop_nt", "n_upstream_atgs",
                "has_protein_support", "best_identity", "frac_introns_supported", "support_level",
                "gffcompare_class") == {
        "strand": "-", "lorf_class": "upLORF", "dist_upstream_stop_nt": "60", "n_upstream_atgs": "0",
        "has_protein_support": "0", "best_identity": "NA", "frac_introns_supported": "NA",
        "support_level": "noSupport", "gffcompare_class": "NA",
    }
    # ranks among all ORFs: lengths 249, 300, 300; alignments 2, 0, 0
    assert [f[t]["cds_length_pct"] for t in "abc"] == ["0.0000", "0.3333", "0.3333"]
    assert [f[t]["n_overlapping_alignments_pct"] for t in "abc"] == ["0.6667", "0.0000", "0.0000"]
    # c starts 1350 nt (in frame) after the end of a
    assert [(f[t]["has_upstream_partner"], f[t]["has_downstream_partner"]) for t in "abc"] == [
        ("0", "1"), ("0", "0"), ("1", "0")]


def test_orf_features_without_evidence(tmp_path: Path):
    """Without miniprot alignments and hints, the protein and hint features are empty."""
    pytest.importorskip("pyfaidx")
    args = write_orf_inputs(tmp_path)
    proc = run("compute_orf_features.py", *args, cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    features = load_features()
    f = table(tmp_path / "features.tsv")
    assert next(iter(f.values())).keys() == set(features.COLUMNS)   # no gffcompare_class
    for tx in "abc":
        assert f[tx]["has_protein_support"] == "0" and f[tx]["best_identity"] == "NA"
        assert f[tx]["n_introns_supported"] == "0" and f[tx]["support_level"] == "noSupport"
    assert f["a"]["dist_upstream_stop_nt"] == "30"


def load_features():
    import importlib.util
    spec = importlib.util.spec_from_file_location("compute_orf_features", BIN / "compute_orf_features.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["compute_orf_features"] = module   # dataclasses look up their module
    spec.loader.exec_module(module)
    return module


def test_split_gene_partners_need_the_frame_and_a_short_gap():
    pytest.importorskip("pyfaidx")
    features = load_features()
    orf = lambda tid, s, e, strand="+": features.ORFRec(tid, "chr1", strand, [(s, e)])
    orfs = [
        orf("a", 0, 300), orf("b", 600, 900),        # gap 300: partners
        orf("c", 1001, 1300),                         # gap 101: other frame
        orf("d", 9000, 9300),                         # gap 7700 to c: too far
        orf("e", 9200, 9500),                         # overlaps d
        orf("f", 9800, 9900, "-"),                    # other strand
    ]
    up, down = features.split_gene_flags(orfs, max_gap=5000)
    assert up == {"b"} and down == {"a"}


# ---------------------------------------------------------------- filter_and_merge_rescue_gtf.py

def tiberius_tx(locus: str, tx: str, strand: str, cds: list[tuple[int, int]]) -> str:
    """Tiberius GTF of one transcript on a rescue locus (local coordinates)."""
    gene = tx.split(".")[0]
    attrs = f'transcript_id "{tx}"; gene_id "{gene}";'
    lines = [f"{locus}\tTiberius\tgene\t{cds[0][0]}\t{cds[-1][1]}\t.\t{strand}\t.\t{gene}\n",
             f"{locus}\tTiberius\ttranscript\t{cds[0][0]}\t{cds[-1][1]}\t.\t{strand}\t.\t{attrs}\n"]
    lines += [f"{locus}\tTiberius\tCDS\t{s}\t{e}\t.\t{strand}\t0\t{attrs}\n" for s, e in cds]
    return "".join(lines)


def test_rescue_predictions_are_filtered_mapped_and_deduplicated(tmp_path: Path):
    (tmp_path / "manifest.tsv").write_text(
        "locus_id\tchr\tbed_start\tbed_end\tstrand\tchain_id\n"
        "locus_0000000\tchr1\t1000\t2000\t+\tc1\n"
        "locus_0000001\tchr1\t1050\t2050\t+\tc2\n"
        "locus_0000002\tchr2\t0\t1000\t-\tnone\n"
        "locus_0000004\tchr3\t0\t1000\t+\tc4\n"
    )
    hint = lambda locus, f, s, e, strand: f"{locus}\tminiprothint\t{f}\t{s}\t{e}\t.\t{strand}\t.\tchain_id=x;\n"
    (tmp_path / "hints.gff").write_text(
        # locus 0: mostly + introns; the hints span 101-800
        hint("locus_0000000", "intron", 201, 300, "+") + hint("locus_0000000", "intron", 501, 600, "+")
        + hint("locus_0000000", "intron", 701, 800, "-") + hint("locus_0000000", "start_codon", 101, 103, "+")
        + hint("locus_0000001", "intron", 151, 250, "+")
        # locus 4: the hints are on the other strand than the locus
        + hint("locus_0000004", "intron", 101, 200, "-") + hint("locus_0000004", "intron", 301, 400, "-")
    )
    (tmp_path / "raw.gtf").write_text(
        tiberius_tx("locus_0000000", "g1.t1", "+", [(101, 200), (301, 500), (601, 700)])
        + tiberius_tx("locus_0000000", "g2.t1", "-", [(301, 400)])         # strand of the locus
        + tiberius_tx("locus_0000000", "g3.t1", "+", [(901, 990)])         # outside the hints
        # the same CDS as g1.t1 of locus 0 in genome coordinates
        + tiberius_tx("locus_0000001", "g1.t1", "+", [(51, 150), (251, 450), (551, 650)])
        + tiberius_tx("locus_0000002", "g1.t1", "-", [(11, 310)])          # locus without hints
        + tiberius_tx("locus_0000003", "g1.t1", "+", [(11, 310)])          # not in the manifest
        + tiberius_tx("locus_0000004", "g1.t1", "+", [(101, 400)])         # strand of the hints
    )
    proc = run("filter_and_merge_rescue_gtf.py", "raw.gtf", "hints.gff", "manifest.tsv", "out.gtf",
               cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    # gene lines have no transcript_id and are not written
    attrs = lambda n: f'transcript_id "rescued_t{n}"; gene_id "rescued_g{n}";'
    assert rows((tmp_path / "out.gtf").read_text()) == [
        ["chr1", "Tiberius", "transcript", "1101", "1700", ".", "+", ".", attrs(1)],
        ["chr1", "Tiberius", "CDS", "1101", "1200", ".", "+", "0", attrs(1)],
        ["chr1", "Tiberius", "CDS", "1301", "1500", ".", "+", "0", attrs(1)],
        ["chr1", "Tiberius", "CDS", "1601", "1700", ".", "+", "0", attrs(1)],
        ["chr2", "Tiberius", "transcript", "11", "310", ".", "-", ".", attrs(2)],
        ["chr2", "Tiberius", "CDS", "11", "310", ".", "-", "0", attrs(2)],
    ]
    assert ("2 wrong-strand | 1 outside hint region | 1 no manifest entry | "
            "1 accepted from chain-less loci") in proc.stderr
    assert "1 duplicate structures removed → 2 unique transcripts" in proc.stderr


# ---------------------------------------------------------------- prepare_hint_rescue_loci.py

needs_samtools = pytest.mark.skipif(shutil.which("samtools") is None, reason="needs samtools")


def chain_hint(feature: str, start: int, end: int, strand: str, chain: str, score: float) -> str:
    """A chainedHints.py line (genome coordinates)."""
    return (f"chr1\tminiprothint\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t"
            f"src=P;al_score={score};chain_id={chain};\n")


CHAINS = [
    chain_hint("start_codon", 501, 503, "+", "A", 0) + chain_hint("intron", 551, 600, "+", "A", 5),
    # B scores higher than A at the merged locus of p1 and p2
    chain_hint("intron", 561, 620, "+", "B", 8) + chain_hint("intron", 661, 700, "+", "B", 8)
    + chain_hint("stop_codon", 898, 900, "+", "B", 0),
    # the best chain, but on the minus strand
    chain_hint("intron", 561, 620, "-", "E", 100),
    # not at any locus
    chain_hint("intron", 1701, 1750, "+", "D", 50),
]


@needs_samtools
def test_rescue_loci(tmp_path: Path):
    rng = random.Random(0)
    genome = "".join(rng.choice("ACGT") for _ in range(2000))
    (tmp_path / "genome.fa").write_text(">chr1\n" + genome + "\n")
    subprocess.run(["samtools", "faidx", "genome.fa"], cwd=tmp_path, check=True)
    tx = lambda s, e, strand: f'chr1\tTiberius\ttranscript\t{s}\t{e}\t.\t{strand}\t.\ttranscript_id "t"; gene_id "g";\n'
    (tmp_path / "partial.gtf").write_text(
        tx(21, 80, "+")            # locus [0, 180): the flank ends at the contig start
        + tx(501, 700, "+") + tx(751, 900, "+")   # flanks overlap: one locus [400, 1000)
        + tx(601, 700, "-")        # the same region on the other strand: its own locus
        + tx(1501, 1600, "+")      # a kept transcript on the same strand overlaps it
        + tx(1951, 2000, "+")      # locus [1850, 2000): the flank ends at the contig end
    )
    (tmp_path / "correct.gtf").write_text(tx(1551, 1700, "+") + tx(501, 560, "-"))
    (tmp_path / "chained.gff").write_text("".join(CHAINS))
    (tmp_path / "orfs.gtf").write_text("")
    proc = run("prepare_hint_rescue_loci.py", "--partial_gtf", "partial.gtf", "--correct_gtf", "correct.gtf",
               "--chained_hints", "chained.gff", "--orfs_gtf", "orfs.gtf", "--genome", "genome.fa",
               "--outdir", "rescue", "--flank", "100", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    out = tmp_path / "rescue"
    manifest = rows((out / "loci_manifest.tsv").read_text())
    assert manifest == [
        ["locus_id", "chr", "bed_start", "bed_end", "strand", "chain_id"],
        ["locus_0000000", "chr1", "0", "180", "+", "none"],
        ["locus_0000001", "chr1", "400", "1000", "+", "B"],
        ["locus_0000002", "chr1", "1850", "2000", "+", "none"],
        ["locus_0000003", "chr1", "500", "800", "-", "E"],
    ]
    assert "1 skipped: tib_correct overlaps same strand" in proc.stderr
    # the sequence of each locus
    fasta = {}
    for block in (out / "combined_loci.fa").read_text().split(">")[1:]:
        name, *lines = block.splitlines()
        assert all(len(line) <= 60 for line in lines)
        fasta[name] = "".join(lines)
    starts = {m[0]: int(m[2]) for m in manifest[1:]}
    assert {name: (starts[name], starts[name] + len(seq)) for name, seq in fasta.items()} == {
        m[0]: (int(m[2]), int(m[3])) for m in manifest[1:]}
    assert all(seq == genome[starts[name]:starts[name] + len(seq)] for name, seq in fasta.items())
    # all hints of the chosen chain in locus coordinates; loci without a chain have none
    hints = rows((out / "combined_hints.gff").read_text())
    assert [(h[0], h[2], h[3], h[4], h[6], h[8]) for h in hints] == [
        ("locus_0000001", "intron", "161", "220", "+", "src=P;al_score=8;chain_id=B;"),
        ("locus_0000001", "intron", "261", "300", "+", "src=P;al_score=8;chain_id=B;"),
        ("locus_0000001", "stop_codon", "498", "500", "+", "src=P;al_score=0;chain_id=B;"),
        ("locus_0000003", "intron", "61", "120", "-", "src=P;al_score=100;chain_id=E;"),
    ]
    for h in hints:
        local, genomic = int(h[3]), int(h[3]) + starts[h[0]]
        assert fasta[h[0]][local - 1:int(h[4])] == genome[genomic - 1:int(h[4]) + starts[h[0]]]
    assert "4 rescue loci written (2 with chain hints, 2 ab-initio)" in proc.stderr


@needs_samtools
def test_no_rescue_loci(tmp_path: Path):
    """All partial transcripts are covered by kept ones: empty outputs, no error."""
    (tmp_path / "genome.fa").write_text(">chr1\n" + "ACGT" * 100 + "\n")
    subprocess.run(["samtools", "faidx", "genome.fa"], cwd=tmp_path, check=True)
    line = 'chr1\tTiberius\ttranscript\t101\t200\t.\t+\t.\ttranscript_id "t"; gene_id "g";\n'
    (tmp_path / "partial.gtf").write_text(line)
    (tmp_path / "correct.gtf").write_text(line)
    (tmp_path / "chained.gff").write_text("")
    proc = run("prepare_hint_rescue_loci.py", "--partial_gtf", "partial.gtf", "--correct_gtf", "correct.gtf",
               "--chained_hints", "chained.gff", "--orfs_gtf", "partial.gtf", "--genome", "genome.fa",
               "--outdir", "rescue", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert (tmp_path / "rescue" / "combined_loci.fa").read_text() == ""
    assert (tmp_path / "rescue" / "loci_manifest.tsv").read_text().count("\n") == 1
