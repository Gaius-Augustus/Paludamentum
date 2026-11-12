"""Every read type as FASTQ, SRA accessions, BAM and/or StringTie assembly, alone and combined (mixed mode).

Stub runs: which steps a combination of inputs starts, and what it publishes. The rules under test
(main.nf): a read type whose assembly is given (rnaseq_stringtie, isoseq_stringtie) is not assembled and
its reads only give the BAM; a mixed run assembles both BAMs together (stringtie --mix) when both read
types are assembled here, and otherwise puts the two assemblies side by side (CONCAT_ASSEMBLIES); a read
type without alignments gets its intron hints from its assembly (ASSEMBLY_HINTS).

Needs a ``nextflow`` executable (or NEXTFLOW_BIN), like test_stub_run.py.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from test_bam_stub_run import tasks, write_bam
from test_stub_run import (BOTH_FINDERS, ISOSEQ, NEXTFLOW, PROTEINS, SHORT_READS, assert_ok, final_files,
                           run_pipeline)

pytestmark = [
    pytest.mark.nextflow,
    pytest.mark.skipif(NEXTFLOW is None, reason="nextflow executable not found"),
]

BASE = {**BOTH_FINDERS, **PROTEINS}


def gtf(tmp_path: Path, name: str) -> str:
    path = tmp_path / name
    path.write_text('c1\tStringTie\ttranscript\t1\t100\t1000\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.1"; '
                    'cov "5.0"; FPKM "1.0"; TPM "2.0";\n'
                    'c1\tStringTie\texon\t1\t100\t1000\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.1"; '
                    'exon_number "1"; cov "5.0";\n')
    return str(path)


def started(tmp_path: Path, *processes: str) -> dict[str, int]:
    return {p: tasks(tmp_path, p) for p in processes}


ASSEMBLERS = ("STRINGTIE_ASSEMBLE_RNA", "STRINGTIE_ASSEMBLE_ISO", "STRINGTIE_ASSEMBLE_MIX", "CONCAT_ASSEMBLIES")


# ---- Iso-Seq BAM -----------------------------------------------------------------

def test_isoseq_bam_alone_is_mode_isoseq_without_minimap2(tmp_path: Path) -> None:
    bam = write_bam(tmp_path / "iso.bam", "coordinate")
    proc, published = run_pipeline(tmp_path, {**BASE, "isoseq_bam": [bam]})
    assert_ok(proc)
    assert "Running mode: isoseq" in proc.stdout
    assert started(tmp_path, "MINIMAP2_MAP", "SAMTOOLS_SORT_ISO", "STRINGTIE_ASSEMBLE_ISO", "BAM2HINTS_ISO", "BAM2COV") == \
        {"MINIMAP2_MAP": 0, "SAMTOOLS_SORT_ISO": 0, "STRINGTIE_ASSEMBLE_ISO": 1, "BAM2HINTS_ISO": 1, "BAM2COV": 1}, proc.stdout
    assert final_files("paludamentum") <= published
    assert "**minimap2**" not in (tmp_path / "out" / "citations.md").read_text()
    assert "Provided Iso-Seq alignments (BAM)" in (tmp_path / "out" / "methods.md").read_text()


def test_unsorted_isoseq_bam_is_sorted_and_merged_with_the_reads(tmp_path: Path) -> None:
    bam = write_bam(tmp_path / "iso.bam", "unsorted")
    proc, _ = run_pipeline(tmp_path, {**BASE, **ISOSEQ, "isoseq_bam": [bam]})
    assert_ok(proc)
    assert started(tmp_path, "MINIMAP2_MAP", "SAMTOOLS_SORT_ISO", "SAMTOOLS_MERGE", "STRINGTIE_ASSEMBLE_ISO") == \
        {"MINIMAP2_MAP": 1, "SAMTOOLS_SORT_ISO": 1, "SAMTOOLS_MERGE": 1, "STRINGTIE_ASSEMBLE_ISO": 1}, proc.stdout


def test_short_and_long_bam_are_a_mixed_run_with_one_joint_assembly(tmp_path: Path) -> None:
    """The request: a short-read BAM and a long-read BAM. Both are alignments, so stringtie --mix assembles
    them together, and the Iso-Seq BAM alone gives the full-length transcripts (chains)."""
    params = {**BASE, "rnaseq_bam": [write_bam(tmp_path / "short.bam", "coordinate")],
              "isoseq_bam": [write_bam(tmp_path / "long.bam", "coordinate")]}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "Running mode: mixed" in proc.stdout
    assert started(tmp_path, "HISAT2_BUILD", "MINIMAP2_MAP", *ASSEMBLERS) == \
        {"HISAT2_BUILD": 0, "MINIMAP2_MAP": 0, "STRINGTIE_ASSEMBLE_RNA": 0, "STRINGTIE_ASSEMBLE_ISO": 1,
         "STRINGTIE_ASSEMBLE_MIX": 1, "CONCAT_ASSEMBLIES": 0}, proc.stdout
    assert tasks(tmp_path, "STRAND_FIX_BAM") == 1      # the short-read BAM is strand-fixed for stringtie --mix
    assert tasks(tmp_path, "BAM2HINTS_RNA") == 1 and tasks(tmp_path, "BAM2HINTS_ISO") == 1
    assert tasks(tmp_path, "ASSEMBLY_HINTS") == 0
    assert final_files("paludamentum") <= published


# ---- assemblies per read type ----------------------------------------------------

def test_short_and_long_assembly_are_a_mixed_run_side_by_side(tmp_path: Path) -> None:
    """The request: a StringTie assembly of the short reads and one of the long reads, no reads at all.
    Nothing is aligned or assembled; the two assemblies are concatenated, and the intron hints and the
    coverage track come from that one assembly."""
    params = {**BASE, "rnaseq_stringtie": gtf(tmp_path, "short.gtf"), "isoseq_stringtie": gtf(tmp_path, "long.gtf")}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "Running mode: mixed" in proc.stdout
    assert started(tmp_path, "HISAT2_BUILD", "MINIMAP2_MAP", "BAM2COV", *ASSEMBLERS) == \
        {"HISAT2_BUILD": 0, "MINIMAP2_MAP": 0, "BAM2COV": 0, "STRINGTIE_ASSEMBLE_RNA": 0, "STRINGTIE_ASSEMBLE_ISO": 0,
         "STRINGTIE_ASSEMBLE_MIX": 0, "CONCAT_ASSEMBLIES": 1}, proc.stdout
    assert tasks(tmp_path, "ASSEMBLY_HINTS") == 1
    assert final_files("paludamentum") <= published
    assert {"intermediate/stringtie_short_and_long.gtf", "intermediate/assembly_introns.gff",
            "intermediate/assembly_cov.bg.gz", "hintsfile.gff"} <= published
    methods = (tmp_path / "out" / "methods.md").read_text()
    assert "the provided StringTie assembly of the short reads together with the provided StringTie assembly of the Iso-Seq reads" in methods
    assert "transcripts of the provided Iso-Seq assembly" in methods        # chains
    assert "Intron hints were derived from the assembly" in methods


@pytest.mark.parametrize("short,long_,expected", [
    # short reads are assembled here, the Iso-Seq assembly is given: its intron hints come from the assembly
    ("reads", "assembly", {"STRINGTIE_ASSEMBLE_RNA": 1, "STRINGTIE_ASSEMBLE_ISO": 0, "STRINGTIE_ASSEMBLE_MIX": 0,
                           "CONCAT_ASSEMBLIES": 1, "BAM2HINTS_RNA": 1, "BAM2HINTS_ISO": 0, "ASSEMBLY_HINTS": 1,
                           "STRAND_FIX_BAM": 1}),
    ("bam", "assembly", {"STRINGTIE_ASSEMBLE_RNA": 1, "STRINGTIE_ASSEMBLE_ISO": 0, "STRINGTIE_ASSEMBLE_MIX": 0,
                         "CONCAT_ASSEMBLIES": 1, "BAM2HINTS_RNA": 1, "BAM2HINTS_ISO": 0, "ASSEMBLY_HINTS": 1,
                         "STRAND_FIX_BAM": 1}),
    # the short-read assembly is given, the Iso-Seq reads are assembled here
    ("assembly", "reads", {"STRINGTIE_ASSEMBLE_RNA": 0, "STRINGTIE_ASSEMBLE_ISO": 1, "STRINGTIE_ASSEMBLE_MIX": 0,
                           "CONCAT_ASSEMBLIES": 1, "BAM2HINTS_RNA": 0, "BAM2HINTS_ISO": 1, "ASSEMBLY_HINTS": 1,
                           "STRAND_FIX_BAM": 0}),
    ("assembly", "bam", {"STRINGTIE_ASSEMBLE_RNA": 0, "STRINGTIE_ASSEMBLE_ISO": 1, "STRINGTIE_ASSEMBLE_MIX": 0,
                         "CONCAT_ASSEMBLIES": 1, "BAM2HINTS_RNA": 0, "BAM2HINTS_ISO": 1, "ASSEMBLY_HINTS": 1,
                         "STRAND_FIX_BAM": 0}),
    # an assembly and the reads it was made from: the reads only give the BAM (hints and coverage)
    ("assembly+reads", "assembly+reads", {"STRINGTIE_ASSEMBLE_RNA": 0, "STRINGTIE_ASSEMBLE_ISO": 0,
                                          "STRINGTIE_ASSEMBLE_MIX": 0, "CONCAT_ASSEMBLIES": 1, "BAM2HINTS_RNA": 1,
                                          "BAM2HINTS_ISO": 1, "ASSEMBLY_HINTS": 0, "STRAND_FIX_BAM": 0}),
    ("assembly+bam", "reads", {"STRINGTIE_ASSEMBLE_RNA": 0, "STRINGTIE_ASSEMBLE_ISO": 1, "STRINGTIE_ASSEMBLE_MIX": 0,
                               "CONCAT_ASSEMBLIES": 1, "BAM2HINTS_RNA": 1, "BAM2HINTS_ISO": 1, "ASSEMBLY_HINTS": 0,
                               "STRAND_FIX_BAM": 0}),
    # no assembly given: both read types are alignments, one joint assembly
    # (no short-read assembly of its own: nothing would read it)
    ("reads", "bam", {"STRINGTIE_ASSEMBLE_RNA": 0, "STRINGTIE_ASSEMBLE_ISO": 1, "STRINGTIE_ASSEMBLE_MIX": 1,
                      "CONCAT_ASSEMBLIES": 0, "BAM2HINTS_RNA": 1, "BAM2HINTS_ISO": 1, "ASSEMBLY_HINTS": 0,
                      "STRAND_FIX_BAM": 1}),
    ("bam+reads", "bam+reads", {"STRINGTIE_ASSEMBLE_RNA": 0, "STRINGTIE_ASSEMBLE_ISO": 1, "STRINGTIE_ASSEMBLE_MIX": 1,
                                "CONCAT_ASSEMBLIES": 0, "BAM2HINTS_RNA": 1, "BAM2HINTS_ISO": 1, "ASSEMBLY_HINTS": 0,
                                "STRAND_FIX_BAM": 2}),
])
def test_mixed_combinations(short: str, long_: str, expected: dict, tmp_path: Path) -> None:
    params = dict(BASE)
    for kind in short.split("+"):
        params.update({"reads": SHORT_READS, "bam": {"rnaseq_bam": [write_bam(tmp_path / "short.bam", "coordinate")]},
                       "assembly": {"rnaseq_stringtie": gtf(tmp_path, "short.gtf")}}[kind])
    for kind in long_.split("+"):
        params.update({"reads": ISOSEQ, "bam": {"isoseq_bam": [write_bam(tmp_path / "long.bam", "coordinate")]},
                       "assembly": {"isoseq_stringtie": gtf(tmp_path, "long.gtf")}}[kind])
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "Running mode: mixed" in proc.stdout
    assert started(tmp_path, *expected) == expected, proc.stdout
    assert tasks(tmp_path, "BAM2COV") == 1          # every combination here has at least one BAM
    assert final_files("paludamentum") <= published


def test_sra_accessions_combine_with_bam_and_assembly(tmp_path: Path) -> None:
    params = {**BASE, "rnaseq_sra_paired": ["SRR000001"], "rnaseq_bam": [write_bam(tmp_path / "short.bam", "coordinate")],
              "isoseq_sra": ["SRR000002"], "isoseq_stringtie": gtf(tmp_path, "long.gtf")}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert "Running mode: mixed" in proc.stdout
    assert started(tmp_path, "DOWNLOAD_SRA_PAIRED", "DOWNLOAD_SRA_ISOSEQ", "MINIMAP2_MAP", "STRINGTIE_ASSEMBLE_ISO",
                   "CONCAT_ASSEMBLIES") == \
        {"DOWNLOAD_SRA_PAIRED": 1, "DOWNLOAD_SRA_ISOSEQ": 1, "MINIMAP2_MAP": 1, "STRINGTIE_ASSEMBLE_ISO": 0,
         "CONCAT_ASSEMBLIES": 1}, proc.stdout
    assert final_files("paludamentum") <= published


def test_joint_assembly_with_given_chains_assembles_nothing_else(tmp_path: Path) -> None:
    """stringtie --mix run with chains given: neither read type needs an assembly of its own."""
    params = {**BASE, **SHORT_READS, **ISOSEQ, "chains": [gtf(tmp_path, "chains.gtf")]}
    proc, published = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert started(tmp_path, *ASSEMBLERS, "STRAND_FIX_BAM") == \
        {"STRINGTIE_ASSEMBLE_RNA": 0, "STRINGTIE_ASSEMBLE_ISO": 0, "STRINGTIE_ASSEMBLE_MIX": 1, "CONCAT_ASSEMBLIES": 0,
         "STRAND_FIX_BAM": 1}, proc.stdout
    assert final_files("paludamentum") <= published


def test_pyvarus_directory_in_a_mixed_run_names_the_other_routes(tmp_path: Path) -> None:
    varus = tmp_path / "varus_short"
    varus.mkdir()
    (varus / "stringtie.gtf").write_text("")
    (varus / "hints.gff").write_text("")
    (varus / "VARUS.manifest.tsv").write_text("#genome_md5=0\n")
    proc, _ = run_pipeline(tmp_path, {**BASE, "rnaseq_varus": [str(varus)], **ISOSEQ})
    assert proc.returncode != 0
    out = " ".join((proc.stdout + proc.stderr).split())
    for hint in ("mixed_varus", "A/VARUS.bam as rnaseq_bam", "B/VARUS.bam as isoseq_bam",
                 "A/stringtie.gtf as rnaseq_stringtie", "B/stringtie.gtf as isoseq_stringtie"):
        assert hint in out, (hint, out[-1500:])


# ---- one read type with its assembly ----------------------------------------------

@pytest.mark.parametrize("key,mode,assembler", [
    ("rnaseq_stringtie", "rnaseq", "STRINGTIE_ASSEMBLE_RNA"), ("isoseq_stringtie", "isoseq", "STRINGTIE_ASSEMBLE_ISO"),
])
def test_one_assembly_of_one_read_type(key: str, mode: str, assembler: str, tmp_path: Path) -> None:
    proc, published = run_pipeline(tmp_path, {**BASE, key: gtf(tmp_path, "a.gtf")})
    assert_ok(proc)
    assert f"Running mode: {mode}" in proc.stdout
    assert started(tmp_path, assembler, "CONCAT_ASSEMBLIES", "ASSEMBLY_HINTS") == \
        {assembler: 0, "CONCAT_ASSEMBLIES": 0, "ASSEMBLY_HINTS": 1}, proc.stdout
    assert final_files("paludamentum") <= published


def test_isoseq_assembly_with_its_bam(tmp_path: Path) -> None:
    params = {**BASE, "isoseq_stringtie": gtf(tmp_path, "long.gtf"),
              "isoseq_bam": [write_bam(tmp_path / "long.bam", "coordinate")]}
    proc, _ = run_pipeline(tmp_path, params)
    assert_ok(proc)
    assert started(tmp_path, "STRINGTIE_ASSEMBLE_ISO", "BAM2HINTS_ISO", "BAM2COV", "ASSEMBLY_HINTS") == \
        {"STRINGTIE_ASSEMBLE_ISO": 0, "BAM2HINTS_ISO": 1, "BAM2COV": 1, "ASSEMBLY_HINTS": 0}, proc.stdout


# ---- combinations that are refused -----------------------------------------------

@pytest.mark.parametrize("extra,message", [
    (lambda t: {"stringtie": gtf(t, "all.gtf"), "rnaseq_stringtie": gtf(t, "short.gtf")},
     "stringtie is the one assembly of all reads"),
    (lambda t: {"rnaseq_stringtie": [gtf(t, "a.gtf"), gtf(t, "b.gtf")]}, "rnaseq_stringtie: one StringTie assembly"),
    (lambda t: {"isoseq_stringtie": [gtf(t, "a.gtf"), gtf(t, "b.gtf")]}, "isoseq_stringtie: one StringTie assembly"),
    (lambda t: {"isoseq_stringtie": str(t / "missing.gtf")}, "isoseq_stringtie: no such file"),
])
def test_refused_assembly_inputs(extra, message: str, tmp_path: Path) -> None:
    proc, _ = run_pipeline(tmp_path, {**BASE, **extra(tmp_path)})
    assert proc.returncode != 0
    assert message in proc.stdout + proc.stderr
