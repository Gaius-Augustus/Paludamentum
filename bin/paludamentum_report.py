#!/usr/bin/env python3
"""
Write the HTML report of a Paludamentum run: one self-contained file (CSS
inline, PNGs embedded as base64 data URIs, no external resources, light and
dark colour scheme by prefers-color-scheme).

STAGED_DIR holds the published files of the run (any of them may be missing;
the section of a missing or empty file is skipped):

    <stem>.gff3 .gtf _proteins.fa _cds.fa _with_ncRNA.gff3 _go.gff3 _with_ncRNA_go.gff3
    citations.md
    qc/sanity_filter.tsv qc/utr_report.tsv qc/completeness.tsv qc/gene_set_statistics.txt
    qc/isoform_and_exon_structure.png qc/transcript_lengths.png qc/introns_per_gene.png
    qc/evidence_support.png qc/gene_support.tsv qc/software_versions.tsv
    qc/omark_summary.txt qc/gffcompare.stats
    qc/fantasia/fantasia_summary.txt qc/fantasia/fantasia_go_categories.png
    ncrna/rRNA.gff3 ncrna/tRNAs.gff3 ncrna/ncRNAs_infernal.gff3 ncrna/lncRNAs.gff3

Sections, in this order: run summary (--run-info), output files, gene set
statistics, completeness, evidence support, sanity filter, UTRs, ncRNA,
OMArk, gffcompare, FANTASIA, software versions, references (citations.md).

If qc/completeness.tsv has rows, a stacked horizontal bar chart of it
(complete single-copy, complete duplicated, fragmented, missing; one bar per
row) is written to --completeness-png (default: completeness.png next to the
report) and embedded. Needs matplotlib; without it the chart is skipped.

The stem (file prefix, e.g. tiberius_evidence) is the 'stem' of run_info.json,
else the prefix of the one '<stem>_proteins.fa' or '<stem>.gff3' in STAGED_DIR.

run_info.json (optional, every key optional): {"version", "mode",
"genefinder", "model", "hc", "stem", "busco_lineage", "outdir"}; further
keys are listed as they are.

Usage:
    paludamentum_report.py --dir STAGED_DIR --out report.html [--run-info run_info.json]
                           [--completeness-png completeness.png]
"""
from __future__ import annotations

import argparse
import base64
import csv
import glob
import html
import json
import os
import re
import sys
from collections import Counter, OrderedDict

# --------------------------------------------------------------------------- files

STEM_FILES = [
    ("{stem}.gff3", "final annotation (GFF3)"),
    ("{stem}.gtf", "final annotation (GTF)"),
    ("{stem}_proteins.fa", "protein sequences of all transcripts"),
    ("{stem}_cds.fa", "coding sequences of all transcripts"),
    ("{stem}_with_ncRNA.gff3", "final annotation plus ncRNA genes"),
    ("{stem}_go.gff3", "final annotation with GO terms (Ontology_term) from FANTASIA"),
    ("{stem}_with_ncRNA_go.gff3", "final annotation plus ncRNA genes, with GO terms"),
]
OTHER_FILES = [
    ("citations.md", "references of the software and data this run used"),
    ("hintsfile.gff", "extrinsic hints (introns, protein alignments)"),
    ("params.yaml", "parameters of the run"),
    ("qc/sanity_filter.tsv", "transcripts removed or changed by the sanity filter, and why"),
    ("qc/utr_report.tsv", "UTRs added from transcript assemblies, per transcript"),
    ("qc/completeness.tsv", "BUSCO and compleasm completeness of genome and proteome"),
    ("qc/completeness.png", "completeness chart"),
    ("qc/gene_set_statistics.txt", "gene set statistics"),
    ("qc/isoform_and_exon_structure.png", "isoforms per gene and exons per transcript"),
    ("qc/transcript_lengths.png", "transcript and CDS lengths"),
    ("qc/introns_per_gene.png", "introns per gene"),
    ("qc/evidence_support.png", "evidence support of the transcripts"),
    ("qc/gene_support.tsv", "evidence support per transcript"),
    ("qc/software_versions.tsv", "software versions and container images"),
    ("qc/omark_summary.txt", "OMArk summary (completeness and consistency of the proteome)"),
    ("qc/gffcompare.stats", "gffcompare accuracy against the reference annotation"),
    ("qc/fantasia/fantasia_summary.txt", "FANTASIA functional annotation summary"),
    ("qc/fantasia/fantasia_go_categories.png", "GO categories of the FANTASIA annotation"),
    ("ncrna/rRNA.gff3", "rRNA genes (barrnap)"),
    ("ncrna/tRNAs.gff3", "tRNA genes (tRNAscan-SE)"),
    ("ncrna/ncRNAs_infernal.gff3", "Rfam ncRNA families (Infernal)"),
    ("ncrna/lncRNAs.gff3", "lncRNA genes"),
]
NCRNA_FILES = ["ncrna/rRNA.gff3", "ncrna/tRNAs.gff3", "ncrna/ncRNAs_infernal.gff3", "ncrna/lncRNAs.gff3"]
SUBFEATURE_TYPES = {"exon", "CDS", "five_prime_UTR", "three_prime_UTR", "UTR", "intron", "start_codon",
                    "stop_codon", "noncoding_exon", "pseudogenic_exon"}
GENE_TYPES = {"gene", "ncRNA_gene", "pseudogene", "rRNA_gene", "tRNA_gene"}
GFFCOMPARE_LEVELS = ["Base level", "Exon level", "Intron level", "Intron chain level", "Transcript level",
                     "Locus level"]
RUN_INFO_LABELS = OrderedDict([
    ("version", "Paludamentum version"), ("mode", "Mode"), ("genefinder", "Gene finder"),
    ("model", "Model"), ("hc", "High-confidence genes"), ("stem", "Output prefix"),
    ("busco_lineage", "BUSCO lineage"), ("outdir", "Output directory"),
])
UNMATCHED = {"", "none", "-", ".", "na", "no"}


def warn(msg: str) -> None:
    print(f"paludamentum_report: {msg}", file=sys.stderr)


def nonempty(path: str) -> bool:
    try:
        return os.path.isfile(path) and os.path.getsize(path) > 0
    except OSError:
        return False


def read_text(path: str) -> str | None:
    """Text of a non-empty file, else None."""
    if not nonempty(path):
        return None
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            text = fh.read()
    except OSError as err:
        warn(f"cannot read {path}: {err}")
        return None
    return text if text.strip() else None


def read_tsv(path: str, expected: list[str] | None = None) -> tuple[list[str], list[str], list[list[str]]]:
    """(comment lines without '#', header, rows) of a TSV with '#' comments.

    With `expected`, the first non-comment line is a header only if it names
    one of the expected columns; otherwise `expected` is the header.
    """
    text = read_text(path)
    if text is None:
        return [], [], []
    comments: list[str] = []
    header: list[str] = []
    rows: list[list[str]] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        if line.startswith("#"):
            comments.append(re.sub(r"^# ?", "", line.rstrip()))   # indentation after '# ' kept
            continue
        fields = line.rstrip("\r\n").split("\t")
        if not header and not rows:
            if expected is None or any(f.strip() in expected for f in fields):
                header = [f.strip() for f in fields]
                continue
            header = list(expected)
        rows.append(fields)
    return comments, header, rows


def column(header: list[str], name: str, default: int) -> int:
    return header.index(name) if name in header else default


def field(row: list[str], idx: int) -> str:
    return row[idx].strip() if 0 <= idx < len(row) else ""


def to_float(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def human_size(n: int) -> str:
    size = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024 or unit == "TB":
            return f"{int(size)} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{n} B"


def find_stem(staged: str, run_info: dict) -> str | None:
    stem = run_info.get("stem")
    if isinstance(stem, str) and stem:
        return stem
    prot = sorted(glob.glob(os.path.join(staged, "*_proteins.fa")))
    if len(prot) == 1:
        return os.path.basename(prot[0])[:-len("_proteins.fa")]
    gffs = [os.path.basename(p)[:-len(".gff3")] for p in sorted(glob.glob(os.path.join(staged, "*.gff3")))]
    gffs = [g for g in gffs if not g.endswith(("_with_ncRNA", "_go", "_with_ncRNA_go"))]
    return gffs[0] if len(gffs) == 1 else None


def load_run_info(path: str | None) -> dict:
    if not path:
        return {}
    text = read_text(path)
    if text is None:
        warn(f"run info {path} missing or empty")
        return {}
    try:
        data = json.loads(text)
    except ValueError as err:
        warn(f"run info {path} is not valid JSON ({err}), ignored")
        return {}
    if not isinstance(data, dict):
        warn(f"run info {path} is not a JSON object, ignored")
        return {}
    return data


# --------------------------------------------------------------------------- HTML helpers

def esc(value) -> str:
    return html.escape(str(value), quote=True)


def table(header: list[str] | None, rows: list[list], numeric: set[int] | None = None) -> str:
    """An HTML table; columns in `numeric` are right-aligned. No header row if header is None."""
    numeric = numeric or set()

    def cls(i: int) -> str:
        return ' class="num"' if i in numeric else ""

    out = ['<div class="scroll"><table>']
    if header is not None:
        out.append("<thead><tr>" + "".join(f"<th{cls(i)}>{esc(h)}</th>" for i, h in enumerate(header))
                   + "</tr></thead>")
    out.append("<tbody>")
    for row in rows:
        out.append("<tr>" + "".join(f"<td{cls(i)}>{esc(v)}</td>" for i, v in enumerate(row)) + "</tr>")
    out.append("</tbody></table></div>")
    return "\n".join(out)


def pre(text: str) -> str:
    return f"<pre>{esc(text.rstrip())}</pre>"


def png(path: str, alt: str) -> str:
    """An <img> with the PNG embedded as a data URI, or '' if missing or empty."""
    if not nonempty(path):
        return ""
    try:
        with open(path, "rb") as fh:
            data = base64.b64encode(fh.read()).decode("ascii")
    except OSError as err:
        warn(f"cannot read {path}: {err}")
        return ""
    return f"<figure><img src=\"data:image/png;base64,{data}\" alt=\"{esc(alt)}\"></figure>"


def fmt_num(value: float) -> str:
    return f"{value:,.1f}" if value != int(value) else f"{int(value):,}"


# --------------------------------------------------------------------------- sections
# Each section function returns its HTML body or '' (section skipped).

def section_run_summary(run_info: dict) -> str:
    rows = [[label, run_info[key]] for key, label in RUN_INFO_LABELS.items()
            if run_info.get(key) not in (None, "")]
    rows += [[key, value] for key, value in run_info.items()
             if key not in RUN_INFO_LABELS and value not in (None, "")]
    rows = [[k, json.dumps(v) if isinstance(v, (dict, list)) else v] for k, v in rows]
    return table(None, rows) if rows else ""


def section_output_files(staged: str, stem: str | None) -> str:
    entries = []
    if stem:
        entries += [(p.format(stem=stem), d) for p, d in STEM_FILES]
    entries += OTHER_FILES
    rows = []
    for rel, desc in entries:
        path = os.path.join(staged, rel)
        if os.path.isfile(path):
            rows.append([rel, desc, human_size(os.path.getsize(path))])
    return table(["File", "Content", "Size"], rows, numeric={2}) if rows else ""


def section_gene_set_statistics(staged: str) -> str:
    parts = []
    text = read_text(os.path.join(staged, "qc", "gene_set_statistics.txt"))
    if text:
        parts.append(pre(text))
    for name, alt in (("isoform_and_exon_structure.png", "Isoforms per gene and exons per transcript"),
                      ("transcript_lengths.png", "Transcript and CDS lengths"),
                      ("introns_per_gene.png", "Introns per gene")):
        parts.append(png(os.path.join(staged, "qc", name), alt))
    return "\n".join(p for p in parts if p)


def read_completeness(path: str) -> tuple[list[str], list[str], list[list[str]]]:
    comments, header, rows = read_tsv(path, expected=["assessment"])
    rows = [r for r in rows if len(r) >= 8]
    return comments, header, rows


def completeness_chart(rows: list[list[str]], out_png: str) -> bool:
    """Stacked horizontal bar chart of the completeness rows; True if written."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        warn("matplotlib not available, completeness chart skipped")
        return False
    parts = [("Complete, single-copy", 3, "#2a78d6"), ("Complete, duplicated", 4, "#1baf7a"),
             ("Fragmented", 5, "#eda100"), ("Missing", 6, "#c3c2b7")]
    labels, values = [], []
    for row in rows:
        vals = [to_float(row[idx]) or 0.0 for _, idx, _ in parts]
        n = row[7].strip()
        labels.append(f"{row[0]} {row[1]}" + (f"\n(n = {n})" if n else ""))
        values.append(vals)
    if not values:
        return False
    fig, ax = plt.subplots(figsize=(9, 1.0 + 0.75 * len(values)), dpi=150)
    ys = list(range(len(values)))
    left = [0.0] * len(values)
    for j, (name, _, colour) in enumerate(parts):
        widths = [v[j] for v in values]
        drawn = [(y, x0, w) for y, x0, w in zip(ys, left, widths) if w > 0]   # no white sliver for 0 %
        ax.barh([d[0] for d in drawn], [d[2] for d in drawn], left=[d[1] for d in drawn], height=0.6,
                color=colour, edgecolor="white", linewidth=1.5)
        for y, x0, w in zip(ys, left, widths):
            if w >= 6:
                ax.text(x0 + w / 2, y, f"{w:.1f}%", ha="center", va="center", fontsize=8, color="#0b0b0b")
        left = [a + b for a, b in zip(left, widths)]
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of BUSCO groups", fontsize=9)
    ax.tick_params(axis="x", labelsize=8, colors="#52514e")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color("#bdbcb6")
    ax.tick_params(axis="y", length=0)
    from matplotlib.patches import Patch   # explicit handles: a part with no bar still gets its colour
    ax.legend(handles=[Patch(facecolor=c, label=name) for name, _, c in parts], loc="lower center",
              bbox_to_anchor=(0.5, 1.0), ncol=4, frameon=False, fontsize=8)
    fig.tight_layout()
    try:
        os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
        fig.savefig(out_png, facecolor="white")
    except OSError as err:
        warn(f"cannot write {out_png}: {err}")
        return False
    finally:
        plt.close(fig)
    return True


def section_completeness(staged: str, out_png: str) -> str:
    comments, header, rows = read_completeness(os.path.join(staged, "qc", "completeness.tsv"))
    if not rows:
        return ""
    parts = []
    if comments:
        parts.append("<ul class=\"notes\">" + "".join(f"<li>{esc(c.strip())}</li>" for c in comments) + "</ul>")
    head = ["Assessment", "Genome or proteome", "Complete %", "Single %", "Duplicated %", "Fragmented %",
            "Missing %", "n"]
    parts.append(table(head, [r[:8] for r in rows], numeric={2, 3, 4, 5, 6, 7}))
    parts.append("<p class=\"note\">compleasm: fragmented = F + I (both fragment classes of compleasm).</p>"
                 if any(r[0].strip() == "compleasm" for r in rows) else "")
    if completeness_chart(rows, out_png):
        parts.append(png(out_png, "Completeness: complete single-copy, complete duplicated, fragmented, missing"))
    return "\n".join(p for p in parts if p)


def is_support_column(name: str) -> bool:
    """Count columns of gene_support.tsv (introns_sup_rnaseq, ...); not the pct_ duplicates."""
    name = name.lower()
    return ("_sup_" in name or "support" in name) and not name.startswith("pct")


def section_evidence_support(staged: str) -> str:
    parts = []
    comments, header, rows = read_tsv(os.path.join(staged, "qc", "gene_support.tsv"))
    if comments:
        parts.append(pre("\n".join(comments)))
    if header or rows:
        parts.append(f"<p>{len(rows):,} transcripts in <code>qc/gene_support.tsv</code>.</p>")
        table_rows = []
        for idx, name in enumerate(header):
            if not is_support_column(name):
                continue
            values = [field(r, idx) for r in rows]
            nums = [to_float(v) for v in values if v not in ("", "NA", ".")]
            if not nums or any(v is None for v in nums):
                continue
            positive = sum(1 for v in nums if v > 0)
            share = f"{100.0 * positive / len(rows):.1f}" if rows else "0"
            table_rows.append([name, f"{positive:,}", share])
        if table_rows:
            parts.append(table(["Column", "Transcripts with value > 0", "% of transcripts"], table_rows,
                               numeric={1, 2}))
    parts.append(png(os.path.join(staged, "qc", "evidence_support.png"), "Evidence support"))
    return "\n".join(p for p in parts if p)


def section_sanity_filter(staged: str) -> str:
    expected = ["gene_id", "transcript_id", "action", "reason"]
    comments, header, rows = read_tsv(os.path.join(staged, "qc", "sanity_filter.tsv"), expected=expected)
    if not comments and not rows:
        return ""
    parts = []
    if comments:
        parts.append(pre("\n".join(comments)))
    a_idx, r_idx = column(header, "action", 2), column(header, "reason", 3)
    counts = Counter((field(r, a_idx), field(r, r_idx)) for r in rows)
    if counts:
        parts.append(table(["Action", "Reason", "Transcripts"],
                           [[a, r, f"{n:,}"] for (a, r), n in sorted(counts.items())], numeric={2}))
    else:
        parts.append("<p>No transcript was removed or changed.</p>")
    return "\n".join(parts)


def section_utrs(staged: str) -> str:
    expected = ["transcript_id", "matched_by", "five_prime_utr_bp", "three_prime_utr_bp"]
    comments, header, rows = read_tsv(os.path.join(staged, "qc", "utr_report.tsv"), expected=expected)
    if not rows:
        return pre("\n".join(comments)) if comments else ""
    m_idx = column(header, "matched_by", 1)
    f_idx, t_idx = column(header, "five_prime_utr_bp", 2), column(header, "three_prime_utr_bp", 3)
    parts = []
    if comments:
        parts.append(pre("\n".join(comments)))
    counts = Counter(field(r, m_idx) or "(empty)" for r in rows)
    parts.append(table(["Matched by", "Transcripts"], [[k, f"{n:,}"] for k, n in sorted(counts.items())],
                       numeric={1}))
    matched = [r for r in rows if field(r, m_idx).lower() not in UNMATCHED]
    if matched:
        stats = []
        for label, idx in (("5' UTR", f_idx), ("3' UTR", t_idx)):
            vals = [v for v in (to_float(field(r, idx)) for r in matched) if v is not None]
            mean = sum(vals) / len(vals) if vals else 0.0
            stats.append([label, f"{sum(1 for v in vals if v > 0):,}", f"{mean:,.1f}"])
        parts.append(f"<p>{len(matched):,} of {len(rows):,} transcripts matched an assembled transcript.</p>")
        parts.append(table(["UTR", "Matched transcripts with this UTR", "Mean length of matched (bp)"],
                           stats, numeric={1, 2}))
    return "\n".join(parts)


def gff3_attr(attrs: str, key: str) -> str | None:
    for part in attrs.split(";"):
        if part.startswith(key + "="):
            return part[len(key) + 1:]
    return None


def count_ncrna(path: str) -> list[tuple[str, int]]:
    """(label, count) of genes and transcript-level features of an ncRNA GFF3."""
    text = read_text(path)
    if text is None:
        return []
    genes = 0
    types: Counter = Counter()
    rrna_names: Counter = Counter()
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 9:
            continue
        ftype, attrs = cols[2], cols[8]
        if ftype in GENE_TYPES:
            genes += 1
        elif ftype not in SUBFEATURE_TYPES:
            types[ftype] += 1
            if ftype == "rRNA":
                name = gff3_attr(attrs, "Name")
                if name:
                    rrna_names[name] += 1
    out = []
    if genes:
        out.append(("genes", genes))
    out += sorted(types.items())
    out += [(f"rRNA {name}", n) for name, n in sorted(rrna_names.items())]
    return out


def section_ncrna(staged: str) -> str:
    rows = []
    for rel in NCRNA_FILES:
        for label, n in count_ncrna(os.path.join(staged, rel)):
            rows.append([rel, label, f"{n:,}"])
    if not rows:
        return ""
    return ("<p>Genes and transcript-level features (column 3 type) per file; for rRNA also by name.</p>\n"
            + table(["File", "Type", "Count"], rows, numeric={2}))


def section_text(path: str) -> str:
    text = read_text(path)
    return pre(text) if text else ""


def section_gffcompare(staged: str) -> str:
    text = read_text(os.path.join(staged, "qc", "gffcompare.stats"))
    if text is None:
        return ""
    level_re = re.compile(r"^\s*(" + "|".join(re.escape(lv) for lv in GFFCOMPARE_LEVELS)
                          + r"):\s*([\d.]+|-)\s*\|\s*([\d.]+|-)\s*\|")
    rows = []
    for line in text.splitlines():
        m = level_re.match(line)
        if m:
            rows.append([m.group(1), m.group(2), m.group(3)])
    parts = []
    if rows:
        parts.append(table(["Level", "Sensitivity %", "Precision %"], rows, numeric={1, 2}))
    parts.append(f"<details><summary>Full gffcompare.stats</summary>{pre(text)}</details>")
    return "\n".join(parts)


def section_fantasia(staged: str) -> str:
    parts = [section_text(os.path.join(staged, "qc", "fantasia", "fantasia_summary.txt")),
             png(os.path.join(staged, "qc", "fantasia", "fantasia_go_categories.png"), "FANTASIA GO categories")]
    return "\n".join(p for p in parts if p)


def section_software_versions(staged: str) -> str:
    text = read_text(os.path.join(staged, "qc", "software_versions.tsv"))
    if text is None:
        return ""
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    if not lines:
        return ""
    data = list(csv.reader(lines, delimiter="\t"))
    header, rows = data[0], data[1:]
    if not rows:
        return ""
    return table(header, rows)


URL_RE = re.compile(r"https?://[^\s<>\"']+")


def inline_md(text: str) -> str:
    """Escaped text with **bold**, *italic*, `code` and URLs as links."""
    out = []
    pos = 0
    for m in URL_RE.finditer(text):
        url = m.group(0)
        trail = ""
        while url and url[-1] in ".,;:)]":
            trail = url[-1] + trail
            url = url[:-1]
        out.append(_inline_style(text[pos:m.start()]))
        out.append(f"<a href=\"{esc(url)}\">{esc(url)}</a>{esc(trail)}")
        pos = m.end()
    out.append(_inline_style(text[pos:]))
    return "".join(out)


def _inline_style(text: str) -> str:
    s = esc(text)
    s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
    s = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", s)
    s = re.sub(r"(?<![\*\w])\*(?=\S)(.+?)(?<=\S)\*(?![\*\w])", r"<em>\1</em>", s)
    return s


def render_citations(text: str) -> str:
    out: list[str] = []
    para: list[str] = []
    items: list[str] = []

    def flush() -> None:
        if para:
            out.append("<p>" + inline_md(" ".join(para)) + "</p>")
            para.clear()
        if items:
            out.append("<ul class=\"refs\">" + "".join(f"<li>{inline_md(i)}</li>" for i in items) + "</ul>")
            items.clear()

    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            flush()
        elif line.startswith("- "):
            if para:
                flush()
            items.append(line[2:].strip())
        elif re.match(r"^#{1,6}\s", line):
            flush()
            tag = "h2" if line.startswith("# ") else "h3"
            out.append(f"<{tag}>" + inline_md(line.lstrip("#").strip()) + f"</{tag}>")
        else:
            if items:
                flush()
            para.append(line)
    flush()
    return "\n".join(out)


def section_references(staged: str) -> str:
    text = read_text(os.path.join(staged, "citations.md"))
    return render_citations(text) if text else ""


# --------------------------------------------------------------------------- page

CSS = """
:root {
  color-scheme: light;
  --bg: #fcfcfb; --surface: #ffffff; --text: #0b0b0b; --muted: #52514e;
  --border: #dedcd5; --accent: #2a78d6; --code-bg: #f3f2ee; --stripe: #f7f6f2;
}
@media (prefers-color-scheme: dark) {
  :root {
    color-scheme: dark;
    --bg: #141413; --surface: #1d1d1c; --text: #f2f1ec; --muted: #b5b4ab;
    --border: #3a3936; --accent: #6aa7f0; --code-bg: #262624; --stripe: #222220;
  }
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--bg); color: var(--text);
  font: 15px/1.55 system-ui, -apple-system, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; }
main { max-width: 1040px; margin: 0 auto; padding: 24px 16px 64px; }
h1 { font-size: 1.7rem; margin: 0 0 4px; }
h1 .stem { color: var(--muted); font-weight: 500; }
h2 { font-size: 1.25rem; margin: 40px 0 12px; padding-bottom: 6px; border-bottom: 1px solid var(--border); }
h3 { font-size: 1.05rem; margin: 20px 0 8px; }
a { color: var(--accent); overflow-wrap: anywhere; }
nav { margin: 12px 0 8px; color: var(--muted); font-size: 0.9rem; }
nav a { margin-right: 12px; white-space: nowrap; }
.scroll { overflow-x: auto; }
table { border-collapse: collapse; margin: 8px 0 16px; font-size: 0.9rem; background: var(--surface); }
th, td { border: 1px solid var(--border); padding: 4px 10px; text-align: left; vertical-align: top; }
th { background: var(--code-bg); font-weight: 600; }
tbody tr:nth-child(even) td { background: var(--stripe); }
td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
pre { background: var(--code-bg); border: 1px solid var(--border); border-radius: 6px; padding: 10px 12px;
  overflow-x: auto; font-size: 0.85rem; line-height: 1.4; }
code { background: var(--code-bg); padding: 1px 4px; border-radius: 4px; font-size: 0.9em; }
figure { margin: 12px 0 20px; }
figure img { display: block; max-width: 100%; height: auto; background: #ffffff; border-radius: 6px;
  border: 1px solid var(--border); padding: 6px; }
details { margin: 8px 0 16px; }
summary { cursor: pointer; color: var(--muted); }
.note, .notes { color: var(--muted); font-size: 0.9rem; }
ul.refs li { margin-bottom: 6px; }
footer { margin-top: 48px; color: var(--muted); font-size: 0.85rem; }
"""


def build_report(staged: str, run_info: dict, completeness_png: str) -> str:
    stem = find_stem(staged, run_info)
    qc = os.path.join(staged, "qc")
    sections = [
        ("run-summary", "Run summary", lambda: section_run_summary(run_info)),
        ("output-files", "Output files", lambda: section_output_files(staged, stem)),
        ("gene-set-statistics", "Gene set statistics", lambda: section_gene_set_statistics(staged)),
        ("completeness", "Completeness", lambda: section_completeness(staged, completeness_png)),
        ("evidence-support", "Evidence support", lambda: section_evidence_support(staged)),
        ("sanity-filter", "Sanity filter", lambda: section_sanity_filter(staged)),
        ("utrs", "UTRs", lambda: section_utrs(staged)),
        ("ncrna", "ncRNA", lambda: section_ncrna(staged)),
        ("omark", "OMArk", lambda: section_text(os.path.join(qc, "omark_summary.txt"))),
        ("gffcompare", "gffcompare", lambda: section_gffcompare(staged)),
        ("fantasia", "FANTASIA", lambda: section_fantasia(staged)),
        ("software-versions", "Software versions", lambda: section_software_versions(staged)),
        ("references", "References", lambda: section_references(staged)),
    ]
    bodies = []
    for anchor, title, func in sections:
        try:
            body = func()
        except Exception as err:   # a broken input must not stop the report
            warn(f"section '{title}' skipped: {type(err).__name__}: {err}")
            body = ""
        if body:
            bodies.append((anchor, title, body))

    title = "Paludamentum report" + (f" – {stem}" if stem else "")
    out = ["<!DOCTYPE html>", "<html lang=\"en\">", "<head>", "<meta charset=\"utf-8\">",
           "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
           "<meta name=\"color-scheme\" content=\"light dark\">",
           f"<title>{esc(title)}</title>", f"<style>{CSS}</style>", "</head>", "<body>", "<main>",
           "<h1>Paludamentum report" + (f" <span class=\"stem\">{esc(stem)}</span>" if stem else "") + "</h1>"]
    if bodies:
        out.append("<nav>" + "".join(f"<a href=\"#{a}\">{esc(t)}</a>" for a, t, _ in bodies) + "</nav>")
    else:
        out.append("<p class=\"note\">No result files were found.</p>")
    for anchor, title_, body in bodies:
        out.append(f"<section id=\"{anchor}\">\n<h2>{esc(title_)}</h2>\n{body}\n</section>")
    out += ["<footer>Written by paludamentum_report.py.</footer>", "</main>", "</body>", "</html>"]
    return "\n".join(out) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dir", required=True, help="directory with the staged result files")
    parser.add_argument("--out", required=True, help="output HTML")
    parser.add_argument("--run-info", help="run_info.json")
    parser.add_argument("--completeness-png",
                        help="where to write the completeness chart (default: completeness.png next to --out)")
    args = parser.parse_args(argv)

    if not os.path.isdir(args.dir):
        warn(f"{args.dir} is not a directory; the report will be empty")
    completeness_png = args.completeness_png or os.path.join(os.path.dirname(os.path.abspath(args.out)),
                                                             "completeness.png")
    report = build_report(args.dir, load_run_info(args.run_info), completeness_png)
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(report)
    print(f"paludamentum_report: wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
