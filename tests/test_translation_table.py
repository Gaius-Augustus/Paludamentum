"""params.translation_table: the scripts that read the genetic code
(gff3_lib.py, sanity_filter_gff3.py, translate_gff.py), the launcher rules and
the wiring of a run with another table than 1 (docs/translation_table.md)."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "bin"
sys.path.insert(0, str(BIN))
sys.path.insert(0, str(ROOT))
import gff3_lib  # noqa: E402
from paludamentum import launcher  # noqa: E402


def run(script: str, *args: str, cwd: Path) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(BIN / script), *args], capture_output=True, text=True, cwd=cwd)


def fasta(path: Path) -> dict[str, str]:
    out, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = ""
        else:
            out[name] += line
    return out


# ---------------------------------------------------------------- the tables

def test_table_1_is_the_standard_code():
    assert gff3_lib.codon_table(1) == gff3_lib.CODON_TABLE
    assert gff3_lib.stop_codons(1) == {"TAA", "TAG", "TGA"}


@pytest.mark.parametrize("table, stops, sense", [
    (6, {"TGA"}, {"TAA": "Q", "TAG": "Q"}),
    (10, {"TAA", "TAG"}, {"TGA": "C"}),
    (15, {"TAA", "TGA"}, {"TAG": "Q"}),
    (29, {"TGA"}, {"TAA": "Y", "TAG": "Y"}),
    (30, {"TGA"}, {"TAA": "E", "TAG": "E"}),
])
def test_reassigned_stop_codons(table, stops, sense):
    assert gff3_lib.stop_codons(table) == stops
    codons = gff3_lib.codon_table(table)
    for codon, amino in sense.items():
        assert codons[codon] == amino
    # nothing else differs from the standard code
    assert {c for c in codons if codons[c] != gff3_lib.CODON_TABLE[c]} == set(sense)


def test_translate_reads_the_table():
    assert gff3_lib.translate("ATGTAAGCTTAGTGA") == "M*A**"
    assert gff3_lib.translate("ATGTAAGCTTAGTGA", 6) == "MQAQ*"
    assert gff3_lib.translate("ATGTAAGCTTAGTGA", 30) == "MEAE*"


def test_unsupported_table():
    with pytest.raises(gff3_lib.GFF3Error, match="translation table 11 is not supported"):
        gff3_lib.codon_table(11)


def test_the_three_lists_of_tables_agree():
    """bin/gff3_lib.py, the launcher and lib_nf/functions.nf name the same tables."""
    tables = sorted(gff3_lib.REASSIGNED_STOPS)
    assert tables == list(launcher.TRANSLATION_TABLES)
    text = (ROOT / "lib_nf" / "functions.nf").read_text()
    assert "return [" + ", ".join(str(t) for t in tables) + "]" in text


# ---------------------------------------------------------------- toy gene

#          ATG GCT TAA GCT TAG GCT GCT ... TGA    in-frame TAA and TAG, stop TGA
CDS = "ATGGCTTAAGCTTAGGCT" + "GCT" * 14 + "TGA"   # 63 nt
GENOME = "C" * 100 + CDS + "C" * 100 + gff3_lib.revcomp(CDS) + "C" * 100
GFF3 = """##gff-version 3
c1\tx\tgene\t101\t163\t.\t+\t.\tID=g1
c1\tx\tmRNA\t101\t163\t.\t+\t.\tID=g1.t1;Parent=g1
c1\tx\texon\t101\t163\t.\t+\t.\tID=e1;Parent=g1.t1
c1\tx\tCDS\t101\t163\t.\t+\t0\tID=c1;Parent=g1.t1
c1\tx\tgene\t264\t326\t.\t-\t.\tID=g2
c1\tx\tmRNA\t264\t326\t.\t-\t.\tID=g2.t1;Parent=g2
c1\tx\texon\t264\t326\t.\t-\t.\tID=e2;Parent=g2.t1
c1\tx\tCDS\t264\t326\t.\t-\t0\tID=c2;Parent=g2.t1
"""
PROTEIN_6 = "MAQAQA" + "A" * 14


@pytest.fixture()
def toy(tmp_path: Path) -> Path:
    (tmp_path / "genome.fa").write_text(f">c1\n{GENOME}\n")
    (tmp_path / "in.gff3").write_text(GFF3)
    return tmp_path


def sanity(toy: Path, *extra: str) -> list[list[str]]:
    proc = run("sanity_filter_gff3.py", "--gff3", "in.gff3", "--genome", "genome.fa", "--out", "out.gff3",
               "--report", "r.tsv", *extra, cwd=toy)
    assert proc.returncode == 0, proc.stderr
    return [line.split("\t") for line in (toy / "r.tsv").read_text().splitlines() if not line.startswith("#")][1:]


def test_sanity_filter_table_1_removes_the_ciliate_genes(toy: Path):
    assert sorted(r[1:] for r in sanity(toy)) == [["g1.t1", "removed", "internal_stop"],
                                                  ["g2.t1", "removed", "internal_stop"]]
    assert "mRNA" not in (toy / "out.gff3").read_text()


def test_sanity_filter_table_6_keeps_them(toy: Path):
    assert sanity(toy, "--translation-table", "6") == []
    assert (toy / "out.gff3").read_text().count("\tmRNA\t") == 2


def test_sanity_filter_extends_the_cds_by_a_stop_codon_of_the_table_only(tmp_path: Path):
    # CDS without its stop codon, followed by TAA (table 6: Gln, no stop) or by TGA
    body = "ATG" + "GCT" * 19
    for name, after in (("taa", "TAA"), ("tga", "TGA")):
        d = tmp_path / name
        d.mkdir()
        (d / "genome.fa").write_text(f">c1\n{'C' * 100}{body}{after}{'C' * 100}\n")
        (d / "in.gff3").write_text(
            "##gff-version 3\nc1\tx\tgene\t101\t170\t.\t+\t.\tID=g\nc1\tx\tmRNA\t101\t170\t.\t+\t.\tID=g.t1;Parent=g\n"
            "c1\tx\texon\t101\t170\t.\t+\t.\tID=e;Parent=g.t1\nc1\tx\tCDS\t101\t160\t.\t+\t0\tID=c;Parent=g.t1\n")
    assert [r[3] for r in sanity(tmp_path / "taa")] == ["extended_stop"]
    assert [r[3] for r in sanity(tmp_path / "taa", "--translation-table", "6")] == ["no_stop"]
    assert [r[3] for r in sanity(tmp_path / "tga", "--translation-table", "6")] == ["extended_stop"]


def test_sanity_filter_refuses_an_unsupported_table(toy: Path):
    proc = run("sanity_filter_gff3.py", "--gff3", "in.gff3", "--genome", "genome.fa", "--out", "out.gff3",
               "--report", "r.tsv", "--translation-table", "11", cwd=toy)
    assert proc.returncode == 1 and "translation table 11 is not supported" in proc.stderr


def test_translate_gff_both_strands(toy: Path):
    proc = run("translate_gff.py", "--gff", "in.gff3", "--genome", "genome.fa", "--translation-table", "6",
               "--out", "p.fa", cwd=toy)
    assert proc.returncode == 0, proc.stderr
    assert fasta(toy / "p.fa") == {"g1.t1": PROTEIN_6, "g2.t1": PROTEIN_6}


def test_translate_gff_reads_gtf_spliced_cds_and_phase(tmp_path: Path):
    # two CDS segments on the minus strand; the first one (in transcription order) has phase 1
    cds = "G" + "CAATAAGCT" + "GCT" * 10 + "TGA"      # phase 1: G | CAA TAA GCT ...
    first, second = cds[:20], cds[20:]
    genome = "C" * 50 + gff3_lib.revcomp(second) + "T" * 30 + gff3_lib.revcomp(first) + "C" * 50
    s2, e2 = 51, 50 + len(second)
    s1, e1 = e2 + 31, e2 + 30 + len(first)
    (tmp_path / "genome.fa").write_text(f">c1\n{genome}\n")
    (tmp_path / "in.gtf").write_text(
        f'c1\tdrusilla\tCDS\t{s2}\t{e2}\t.\t-\t.\ttranscript_id "tx1"; gene_id "g1";\n'
        f'c1\tdrusilla\tCDS\t{s1}\t{e1}\t.\t-\t1\ttranscript_id "tx1"; gene_id "g1";\n')
    proc = run("translate_gff.py", "--gff", "in.gtf", "--genome", "genome.fa", "--translation-table", "6",
               "--out", "p.fa", cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert fasta(tmp_path / "p.fa") == {"tx1": "QQA" + "A" * 10}


GFFREAD = shutil.which("gffread")


@pytest.mark.skipif(GFFREAD is None, reason="gffread not found")
def test_translate_gff_writes_what_gffread_writes_in_table_1(toy: Path):
    """PROTEIN_FROM_GFF and FINALIZE use gffread -y -S with table 1 and translate_gff.py otherwise."""
    subprocess.run([GFFREAD, "in.gff3", "-g", "genome.fa", "-y", "gffread.fa", "-S"], cwd=toy, check=True)
    proc = run("translate_gff.py", "--gff", "in.gff3", "--genome", "genome.fa", "--translation-table", "1",
               "--out", "p.fa", cwd=toy)
    assert proc.returncode == 0, proc.stderr
    assert fasta(toy / "p.fa") == fasta(toy / "gffread.fa")
    assert "*" in fasta(toy / "p.fa")["g1.t1"]


# ---------------------------------------------------------------- launcher

def params(**kw) -> dict:
    base = {"genome": "g.fa", "tiberius": {}, "vipsania": {}}
    base.update(kw)
    return base


def test_launcher_default_table_is_1():
    assert launcher.translation_table(params()) == 1
    assert launcher.translation_table(params(translation_table="6")) == 6


@pytest.mark.parametrize("value", [11, 4, "six", 0])
def test_launcher_refuses_an_unsupported_table(value):
    with pytest.raises(SystemExit, match="translation_table"):
        launcher.translation_table(params(translation_table=value))


def test_vipsania_needs_no_model_with_another_table():
    assert launcher.gene_finder_models(params()) == (False, False)
    assert launcher.gene_finder_models(params(translation_table=6)) == (False, True)
    evidence = {"proteins": ["p.fa"], "rnaseq_bam": ["a.bam"]}
    assert launcher.check_mode(params(translation_table=6, **evidence)) == "rnaseq"
    with pytest.raises(SystemExit, match="vipsania.model"):
        launcher.check_mode(params(**evidence))


def test_tiberius_model_is_refused_with_another_table():
    with pytest.raises(SystemExit, match="Tiberius has no model for this genetic code"):
        launcher.check_mode(params(translation_table=6, proteins=["p.fa"], rnaseq_bam=["a.bam"],
                                   tiberius={"model_cfg": "diatoms"}))


def test_no_tiberius_model_is_derived_with_another_table():
    p = params(translation_table=6, proteins=["p.fa"], rnaseq_bam=["a.bam"], vipsania={"model": "Fungi"})
    assert launcher.complete_models(p) == []
    assert not p["tiberius"].get("model_cfg")


# ---------------------------------------------------------------- the modules

@pytest.mark.parametrize("module, text", [
    ("modules/transdecoder.nf", '-@ ${task.cpus}${gcode}'),
    ("modules/transdecoder.nf", 'def gcode = table == 1 ? \'\' : " -G ${table}"'),
    ("modules/drusilla.nf", 'extra += " --translation-table ${table}"'),
    ("modules/proteins.nf", '${params.tools.miniprot} -t ${task.cpus}${gcode} --aln'),
    ("modules/proteins.nf", 'def gcode = table == 1 ? \'\' : " -T ${table}"'),
    ("modules/finalize.nf", '--translation-table ${translationTable(params)}'),
    ("modules/vipsania.nf", 'vipsania train "\\$config" --translation_table ${table}'),
    ("modules/vipsania.nf", '--translation_table ${table} \\\\\n        -o "vipsania.${genome.name}.gtf"'),
])
def test_the_tools_get_the_table(module: str, text: str):
    assert text in (ROOT / module).read_text()


def test_no_protein_file_is_written_by_gffread_alone():
    """Every protein FASTA of the pipeline comes from proteinsCommand (gffread in table 1 only)."""
    for module in sorted((ROOT / "modules").glob("*.nf")) + sorted((ROOT / "subworkflows").glob("*.nf")):
        for line in module.read_text().splitlines():
            if "gffread" in line and " -y " in line:
                pytest.fail(f"{module.name}: {line.strip()}")


# ---- compleasm and BUSCO in genome mode: miniprot -T through a wrapper --------


def fake_miniprot(directory: Path, with_table: bool) -> Path:
    directory.mkdir()
    tool = directory / "miniprot"
    usage = "    -T INT       NCBI translation table" if with_table else "    -t INT       threads"
    tool.write_text(f'#!/bin/sh\nif [ $# -eq 0 ]; then echo "{usage}" >&2; exit 1; fi\necho "$@"\n')
    tool.chmod(0o755)
    return tool


def run_wrapper_script(tmp_path: Path, path: str, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", str(ROOT / "bin" / "miniprot_table_wrapper.sh"), "6", str(tmp_path / "w"), *args],
                          env={"PATH": path}, capture_output=True, text=True)


def test_miniprot_table_wrapper_adds_the_table(tmp_path: Path) -> None:
    fake_miniprot(tmp_path / "tools", with_table=True)
    proc = run_wrapper_script(tmp_path, f"{tmp_path / 'tools'}:/usr/bin:/bin")
    assert proc.returncode == 0, proc.stderr
    out = subprocess.run([str(tmp_path / "w" / "miniprot"), "--trans", "g.fa", "p.fa"], capture_output=True, text=True)
    assert out.stdout.strip() == "-T 6 --trans g.fa p.fa"


def test_miniprot_table_wrapper_finds_miniprot_beside_compleasm(tmp_path: Path) -> None:
    fake_miniprot(tmp_path / "kit", with_table=True)
    proc = run_wrapper_script(tmp_path, "/usr/bin:/bin", str(tmp_path / "kit"))
    assert proc.returncode == 0, proc.stderr
    assert str(tmp_path / "kit" / "miniprot") in (tmp_path / "w" / "miniprot").read_text()


def test_miniprot_table_wrapper_refuses_a_miniprot_without_the_option(tmp_path: Path) -> None:
    fake_miniprot(tmp_path / "tools", with_table=False)
    proc = run_wrapper_script(tmp_path, f"{tmp_path / 'tools'}:/usr/bin:/bin")
    assert proc.returncode != 0
    assert "has no option -T" in proc.stderr
    assert not (tmp_path / "w" / "miniprot").exists()


def test_completeness_of_the_genome_passes_the_table() -> None:
    text = (ROOT / "modules" / "completeness.nf").read_text()
    assert text.count("miniprot_table_wrapper.sh ${table}") == 2
    assert "--miniprot_execute_path" in text
