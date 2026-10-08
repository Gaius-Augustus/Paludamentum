"""Unit tests of the pipeline scripts in bin/ on toy inputs. No tools needed."""
from __future__ import annotations

import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "bin"


def load(name: str):
    spec = importlib.util.spec_from_file_location(name, BIN / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run(script: str, *args: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(BIN / script), *args], capture_output=True, text=True, cwd=cwd)


def rows(text: str) -> list[list[str]]:
    return [line.split("\t") for line in text.splitlines() if line and not line.startswith("#")]


# ---------------------------------------------------------------- split_genome_fasta.py

def test_chunks_are_balanced_and_capped():
    sg = load("split_genome_fasta")
    info = [("", f"s{i}", 50_000_000) for i in range(60)]   # 60 x 50 Mb, max 20 chunks
    groups = sg.make_groups(info, 20_000_000, 20)
    assert len(groups) == 20
    assert {len(g) for g in groups} == {3}                  # 150 Mb each, not one 2 Gb chunk
    assert sg.make_groups([("", "a", 5), ("", "b", 5)], 100, 20) == [["a", "b"]]
    seven = [("", f"s{i}", 1000) for i in range(7)]
    assert len(sg.make_groups(seven, 1000, 3)) == 3


# ---------------------------------------------------------------- check_stop_codons.py

def test_passing_gff_keeps_gene_and_transcript_lines(tmp_path: Path):
    (tmp_path / "g.fa").write_text(">c1\nATGAAATAAATGCCCCCC\n")
    (tmp_path / "in.gff").write_text(
        "##gff-version 3\n"
        "c1\tx\tgene\t1\t18\t.\t+\t.\tID=g1\n"
        "c1\tx\ttranscript\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n"
        "c1\tx\texon\t1\t9\t.\t+\t.\tID=t1.e1;Parent=t1\n"
        "c1\tx\tCDS\t1\t9\t.\t+\t0\tID=t1.c1;Parent=t1\n"
        "c1\tx\ttranscript\t10\t18\t.\t+\t.\tID=t2;Parent=g1\n"
        "c1\tx\tCDS\t10\t18\t.\t+\t0\tID=t2.c1;Parent=t2\n"       # no stop codon
        "c1\tx\tgene\t10\t18\t.\t-\t.\tID=g2\n"
        "c1\tx\ttranscript\t10\t18\t.\t-\t.\tID=t3;Parent=g2\n"
        "c1\tx\tCDS\t10\t18\t.\t-\t0\tID=t3.c1;Parent=t3\n"
    )
    proc = run("check_stop_codons.py", "in.gff", "g.fa", "--write-passing-gff", "out.gff", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    out = rows((tmp_path / "out.gff").read_text())
    assert [(r[2], r[8]) for r in out] == [
        ("gene", "ID=g1"), ("transcript", "ID=t1;Parent=g1"),
        ("exon", "ID=t1.e1;Parent=t1"), ("CDS", "ID=t1.c1;Parent=t1"),
    ]


# ---------------------------------------------------------------- extend_cds_with_stop_codon.py

# ATGAAA, stop TAA split as T | intron GTAAGTTTAG | AA, then CCC
SPLIT_STOP_GENOME = "ATGAAAT" + "GTAAGTTTAG" + "AACCC"
SPLIT_STOP_GENOME_RC = SPLIT_STOP_GENOME[::-1].translate(str.maketrans("ACGT", "TGCA"))


def split_stop_gff(strand: str, cds: list[tuple[int, int]], stops: list[tuple[int, int]]) -> str:
    lo = min(s for s, _ in cds + stops)
    hi = max(e for _, e in cds + stops)
    text = (f"c1\tx\tgene\t{lo}\t{hi}\t.\t{strand}\t.\tID=g1\n"
            f"c1\tx\tmRNA\t{lo}\t{hi}\t.\t{strand}\t.\tID=t1;Parent=g1\n")
    text += "".join(f"c1\tx\tCDS\t{s}\t{e}\t.\t{strand}\t0\tID=cds.t1;Parent=t1\n" for s, e in cds)
    text += "".join(f"c1\tx\tstop_codon\t{s}\t{e}\t.\t{strand}\t0\tParent=t1\n" for s, e in stops)
    return text


@pytest.mark.parametrize("strand, genome, cds, stops, want_cds, want_exons", [
    # stop codon split by an intron: one CDS segment on each side, no CDS in the intron
    ("+", SPLIT_STOP_GENOME, [(1, 6)], [(7, 7), (18, 19)],
     [(1, 7, "0"), (18, 19, "2")], [(1, 7), (18, 19)]),
    ("-", SPLIT_STOP_GENOME_RC, [(17, 22)], [(16, 16), (4, 5)], [(4, 5, "2"), (16, 22, "0")], [(4, 5), (16, 22)]),
    # stop codon next to the CDS, and one already inside it
    ("+", "ATGAAATAACCC", [(1, 6)], [(7, 9)], [(1, 9, "0")], [(1, 9)]),
    ("+", "ATGAAATAACCC", [(1, 9)], [(7, 9)], [(1, 9, "0")], [(1, 9)]),
])
def test_stop_codon_is_added_to_the_cds(tmp_path: Path, strand, genome, cds, stops, want_cds, want_exons):
    (tmp_path / "g.fa").write_text(f">c1\n{genome}\n")
    (tmp_path / "in.gff").write_text(split_stop_gff(strand, cds, stops))
    proc = run("extend_cds_with_stop_codon.py", "in.gff", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    out = rows(proc.stdout)
    assert [(int(r[3]), int(r[4]), r[7]) for r in out if r[2] == "CDS"] == want_cds
    assert [(int(r[3]), int(r[4])) for r in out if r[2] == "exon"] == want_exons
    (tmp_path / "ext.gff").write_text(proc.stdout)
    proc = run("check_stop_codons.py", "ext.gff", "g.fa", "--write-passing-gff", "out.gff", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert any(r[2] == "CDS" for r in rows((tmp_path / "out.gff").read_text()))   # valid stop, kept


# ---------------------------------------------------------------- merge_annotations.py

def merge(tmp_path: Path, *contents: str) -> subprocess.CompletedProcess:
    inputs = []
    for i, text in enumerate(contents):
        path = tmp_path / f"in{i}.gtf"
        path.write_text(text)
        inputs.append(str(path))
    return run("merge_annotations.py", "--mode", "full", *inputs)


HOST = (
    'chr1\tTiberius\texon\t100\t200\t.\t+\t.\tgene_id "h"; transcript_id "h.t1";\n'
    'chr1\tTiberius\tCDS\t100\t200\t.\t+\t0\tgene_id "h"; transcript_id "h.t1";\n'
    'chr1\tTiberius\texon\t900\t1000\t.\t+\t.\tgene_id "h"; transcript_id "h.t1";\n'
    'chr1\tTiberius\tCDS\t900\t1000\t.\t+\t1\tgene_id "h"; transcript_id "h.t1";\n'
)
NESTED = (
    'chr1\tTiberius\texon\t400\t500\t.\t+\t.\tgene_id "n"; transcript_id "n.t1";\n'
    'chr1\tTiberius\tCDS\t400\t500\t.\t+\t0\tgene_id "n"; transcript_id "n.t1";\n'
)
HOST_WITH_UTR = (
    'chr1\tTransDecoder\texon\t50\t200\t.\t+\t.\tgene_id "x"; transcript_id "x.t1";\n'
    'chr1\tTransDecoder\tfive_prime_UTR\t50\t99\t.\t+\t.\tgene_id "x"; transcript_id "x.t1";\n'
    'chr1\tTransDecoder\tCDS\t100\t200\t.\t+\t0\tgene_id "x"; transcript_id "x.t1";\n'
    'chr1\tTransDecoder\texon\t900\t1000\t.\t+\t.\tgene_id "x"; transcript_id "x.t1";\n'
    'chr1\tTransDecoder\tCDS\t900\t1000\t.\t+\t1\tgene_id "x"; transcript_id "x.t1";\n'
)


def test_gene_nested_in_an_intron_is_its_own_gene(tmp_path: Path):
    out = rows(merge(tmp_path, HOST + NESTED).stdout)
    genes = [r for r in out if r[2] == "gene"]
    assert [(r[3], r[4]) for r in genes] == [("100", "1000"), ("400", "500")]
    mrnas = [r for r in out if r[2] == "mRNA"]
    assert len(mrnas) == 2 and len({r[8].split(";")[1] for r in mrnas}) == 2   # different Parents


def test_same_cds_with_and_without_utr_is_one_transcript(tmp_path: Path):
    out = rows(merge(tmp_path, HOST, HOST_WITH_UTR).stdout)
    assert len([r for r in out if r[2] == "mRNA"]) == 1
    assert [r for r in out if r[2] == "five_prime_UTR"]          # the UTR version is kept
    assert [r[3] for r in out if r[2] == "gene"] == ["50"]


def test_overlapping_exons_still_form_one_gene(tmp_path: Path):
    other = 'chr1\tX\texon\t150\t300\t.\t+\t.\tgene_id "o"; transcript_id "o.t1";\n'
    out = rows(merge(tmp_path, HOST, other).stdout)
    assert len([r for r in out if r[2] == "gene"]) == 1
    assert len([r for r in out if r[2] in ("mRNA", "transcript")]) == 2


def test_no_transcripts_is_an_empty_gff3_not_an_error(tmp_path: Path):
    proc = merge(tmp_path, "")
    assert proc.returncode == 0
    assert proc.stdout == "##gff-version 3\n"
    assert "no transcripts" in proc.stderr


# ---------------------------------------------------------------- rank_species_from_diamond.py

def diamond_line(q: str, s: str) -> str:
    # qseqid sseqid pident length evalue bitscore qlen slen
    return f"{q}\t{s}\t90\t100\t1e-50\t200\t100\t100\n"


def test_species_ranking_only_for_orthodb_ids(tmp_path: Path):
    hits = tmp_path / "hits.tsv"
    hits.write_text(diamond_line("q1", "9606_0:00001") + diamond_line("q2", "9606_0:00002") + diamond_line("q3", "10090_0:00003"))
    proc = run("rank_species_from_diamond.py", str(hits), "1", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert (tmp_path / "top_species.txt").read_text() == "9606\n"

    hits.write_text(diamond_line("q1", "sp|P12345|NAME_HUMAN") + diamond_line("q2", "sp|P12346|OTHER_HUMAN"))
    proc = run("rank_species_from_diamond.py", str(hits), "5", cwd=tmp_path)
    assert proc.returncode == 0
    # no grouping by the text before the first '_': every id is its own entry
    assert set((tmp_path / "top_species.txt").read_text().split()) == {"sp|P12345|NAME_HUMAN", "sp|P12346|OTHER_HUMAN"}

    hits.write_text("")
    proc = run("rank_species_from_diamond.py", str(hits), "5", cwd=tmp_path)
    assert proc.returncode == 0 and (tmp_path / "top_species.txt").read_text() == ""
    assert "empty" in proc.stderr


# ---------------------------------------------------------------- hc_module.py

def test_transcript_gff3_has_distinct_gene_and_mrna_ids(tmp_path: Path):
    pytest.importorskip("Bio")
    sys.path.insert(0, str(BIN))
    hc_module = load("hc_module")
    pep = tmp_path / "hc.pep"
    pep.write_text(">STRG.1.1.p1 STRG.1.1.p1 GENE.STRG.1.1~~STRG.1.1.p1 ORF type:complete len:4 (+),score=1 STRG.1.1:10-21(+)\nMKR*\n")
    gtf = tmp_path / "st.gtf"
    gtf.write_text('c1\tStringTie\texon\t1\t50\t.\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.1";\n')
    out = tmp_path / "out.gff3"
    hc_module.from_pep_file_to_gff3(str(pep), str(gtf), str(out))
    lines = rows(out.read_text())
    gene = next(r for r in lines if r[2] == "gene")
    mrna = next(r for r in lines if r[2] == "mRNA")
    gene_id = gene[8].split(";")[0].removeprefix("ID=")
    assert mrna[8].startswith("ID=STRG.1.1.p1;Parent=" + gene_id)
    assert gene_id != "STRG.1.1.p1"


# ---------------------------------------------------------------- filterIntronsFindStrand.pl

@pytest.mark.skipif(shutil.which("perl") is None, reason="perl not found")
def test_intron_strand_with_described_fasta_header(tmp_path: Path):
    (tmp_path / "g.fa").write_text(">chr1 Homo sapiens chromosome 1, GRCh38\nAAAGTAAGTCCCCCCCCCAGGGG\n")
    (tmp_path / "i.gff").write_text("chr1\tb2h\tintron\t4\t20\t1\t.\t.\tmult=3;src=E\n")
    proc = subprocess.run(["perl", str(BIN / "filterIntronsFindStrand.pl"), "g.fa", "i.gff", "--score"],
                          cwd=tmp_path, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert rows(proc.stdout) == [["chr1", "b2h", "intron", "4", "20", "3", "+", ".", "mult=3;src=E"]]
    assert "does not match" not in proc.stderr


# ---------------------------------------------------------------- merge_intron_hints.py

def hint(seq: str, start: int, end: int, mult: int, strand: str = "+") -> str:
    """A hint line as BAM2HINTS (and pyVARUS hints.gff) writes it."""
    attr = (f"mult={mult};" if mult > 1 else "") + "pri=4;src=E"
    return f"{seq}\tb2h\tintron\t{start}\t{end}\t{mult}\t{strand}\t.\t{attr}\n"


def test_intron_hints_are_summed_as_on_the_merged_bam(tmp_path: Path):
    (tmp_path / "a.gff").write_text(
        hint("chr2", 50, 90, 2, "-") + hint("chr1", 121, 220, 1) + hint("chr1", 300, 400, 3)
    )
    (tmp_path / "b.gff").write_text(
        "# comment\n\n" + hint("chr1", 121, 220, 2) + hint("chr1", 100, 220, 1) + hint("chr3", 5, 60, 1)
        + hint("chr2", 50, 90, 1, "+")
    )
    proc = run("merge_intron_hints.py", "a.gff", "b.gff", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout == (
        hint("chr2", 50, 90, 1, "+") + hint("chr2", 50, 90, 2, "-")
        + hint("chr1", 100, 220, 1) + hint("chr1", 121, 220, 3) + hint("chr1", 300, 400, 3)
        + hint("chr3", 5, 60, 1)
    )


def test_one_hint_file_is_unchanged(tmp_path: Path):
    text = hint("chr1", 121, 220, 3) + hint("chr1", 500, 800, 1, "-") + hint("chr2", 10, 70, 1)
    (tmp_path / "a.gff").write_text(text)
    proc = run("merge_intron_hints.py", "a.gff", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout == text


def test_malformed_hint_line_is_an_error(tmp_path: Path):
    (tmp_path / "a.gff").write_text("chr1\tb2h\tintron\t1\t2\n")
    proc = run("merge_intron_hints.py", "a.gff", cwd=tmp_path)
    assert proc.returncode != 0
    assert "a.gff:1" in proc.stderr


# ---------------------------------------------------------------- compleasm_wrapper.py

# Stand-in for compleasm.py: its own URLError class and the hash download of
# Downloader.download_file_version_document, which writes into the library path
FAKE_COMPLEASM = '''
import sys, urllib.request
class URLError(OSError):
    pass
lib = sys.argv[sys.argv.index("-L") + 1]
try:
    urllib.request.urlretrieve(sys.argv[-1], lib + "/file_versions.tsv.hash")
    print("downloaded")
except URLError:
    print("offline: using cached file_versions.tsv")
'''


@pytest.mark.parametrize("mode", ["run", "protein"])
def test_compleasm_wrapper_never_downloads_in_run_and_protein_mode(tmp_path, mode):
    # tasks sharing a library path overwrote file_versions.tsv.hash while
    # another task read it (IndexError); in these modes nothing is written
    fake = tmp_path / "compleasm.py"
    fake.write_text(FAKE_COMPLEASM)
    lib = tmp_path / "lib"
    lib.mkdir()
    remote = tmp_path / "remote.hash"
    remote.write_text("0123456789abcdef  file_versions.tsv\n")
    r = run("compleasm_wrapper.py", str(fake), mode, "-L", str(lib), remote.as_uri())
    assert r.returncode == 0, r.stderr
    assert "offline" in r.stdout
    assert not (lib / "file_versions.tsv.hash").exists()


def test_compleasm_wrapper_downloads_in_download_mode(tmp_path):
    fake = tmp_path / "compleasm.py"
    fake.write_text(FAKE_COMPLEASM)
    lib = tmp_path / "lib"
    lib.mkdir()
    remote = tmp_path / "remote.hash"
    remote.write_text("0123456789abcdef  file_versions.tsv\n")
    r = run("compleasm_wrapper.py", str(fake), "download", "vertebrata", "-L", str(lib), remote.as_uri())
    assert r.returncode == 0, r.stderr
    assert "downloaded" in r.stdout
    assert (lib / "file_versions.tsv.hash").read_text() == remote.read_text()


def test_compleasm_wrapper_extracts_members_in_archive_order(tmp_path, monkeypatch):
    # compleasm passes the HMMs in set order; in a .tar.gz each backward seek
    # decompresses the archive again, so the wrapper sorts them by offset
    import tarfile
    names = [f"lin/hmms/g{i}.hmm" for i in range(5)]
    archive = tmp_path / "lin.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        for n in names:
            f = tmp_path / n.replace("/", "_")
            f.write_text(n)
            tar.add(f, arcname=n)
    seen = []
    monkeypatch.setattr(tarfile.TarFile, "extractall",
                        lambda self, path=".", members=None, **kw: seen.extend(m.name for m in members))
    load("compleasm_wrapper").patch_tarfile_member_order()
    with tarfile.open(archive) as tar:
        tar.extractall(tmp_path / "out", members=[tar.getmember(n) for n in reversed(names)])
    assert seen == names
