#!/usr/bin/env python3
# Port of BRAKER4 scripts/gene_set_statistics.py (commit 3535ed3), see the docstring.
# Copyright (c) 2025 Katharina Hoff. MIT License, see LICENSE-BRAKER4.
"""
Gene set statistics and plots of a GFF3 gene set.

Writes to --out-dir:
  gene_set_statistics.txt          genes, transcripts, mean transcripts per
                                   gene and exons per transcript, single- and
                                   multi-exon transcripts, mean and median CDS
                                   length and genomic span
  isoform_and_exon_structure.png   transcripts per gene, single- vs multi-exon
  transcript_lengths.png           CDS length and genomic span histograms
  introns_per_gene.png             introns per gene histogram
  evidence_support.png             only with --support (gene_support.tsv of
                                   gene_support.py): fraction of supported
                                   introns per multi-exon transcript and
                                   full/partial/no support

Exons of a transcript are its exon features, its CDS segments if it has no
exons; introns are the gaps between them. The CDS length is the sum of the
CDS segments (of the exons if the transcript has no CDS); the genomic span
runs from the first to the last exon. Introns per gene = the introns of all
its transcripts divided by the number of transcripts, rounded down. Medians
are the upper median (the element at n // 2 of the sorted list). The plots
need matplotlib (Agg backend); without it only the text file is written.

Port of BRAKER4 scripts/gene_set_statistics.py (commit 3535ed3) to GFF3;
changes: GFF3 input read with gff3_lib instead of the GTF regexes; introns
per gene from the exons or, without exons, the CDS segments (BRAKER4: exon
features only); title line "Paludamentum Gene Set Statistics"; percentages
NA instead of a ZeroDivisionError on an empty gene set; PNG only (no PDF);
every plot is always written (an empty panel says "no data") so that the
output files do not depend on the input; the evidence support plot leaves
out single-exon transcripts (BRAKER4 counted them as "No support"); no numpy;
messages on stderr; options --gff3/--out-dir/--support.

Usage:
    gene_set_statistics.py --gff3 final.gff3 --out-dir DIR [--support gene_support.tsv]
"""
from __future__ import annotations

import argparse
import math
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gff3_lib import GFF3Error, read_gff3  # noqa: E402

PLOT_FILES = ["isoform_and_exon_structure.png", "transcript_lengths.png", "introns_per_gene.png"]
SUPPORT_PLOT = "evidence_support.png"


def upper_median(values):
    return sorted(values)[len(values) // 2]


def compute_statistics(ann):
    """Summary statistics of the gene set (as BRAKER4's compute_statistics)."""
    tx_per_gene = []
    exon_counts = []
    single_exon_tx = 0
    multi_exon_tx = 0
    cds_lengths = []      # sum of the CDS segment lengths
    genomic_spans = []    # first exon start to last exon end
    introns_per_gene = []

    for gene in ann.genes:
        tx_per_gene.append(len(gene.transcripts))
        gene_introns = 0
        for tx in gene.transcripts:
            exons = [(e.start, e.end) for e in tx.exons]
            cds = [(c.start, c.end) for c in tx.cds]
            parts = exons if exons else cds
            n_exons = len(parts)
            if n_exons <= 1:
                single_exon_tx += 1
            else:
                multi_exon_tx += 1
            exon_counts.append(n_exons)

            cds_len = sum(e - s + 1 for s, e in (cds if cds else exons))
            if cds_len > 0:
                cds_lengths.append(cds_len)
            if parts:
                genomic_spans.append(max(e for _s, e in parts) - min(s for s, _e in parts) + 1)
            if n_exons > 1:
                gene_introns += n_exons - 1
        if gene.transcripts:
            introns_per_gene.append(gene_introns // len(gene.transcripts))

    def mean(values):
        return sum(values) / len(values) if values else 0

    return {
        "n_genes": len(ann.genes),
        "n_transcripts": sum(tx_per_gene),
        "tx_per_gene": tx_per_gene,
        "tx_per_gene_counts": Counter(tx_per_gene),
        "single_exon_tx": single_exon_tx,
        "multi_exon_tx": multi_exon_tx,
        "exon_counts": exon_counts,
        "cds_lengths": cds_lengths,
        "genomic_spans": genomic_spans,
        "introns_per_gene": introns_per_gene,
        "mean_tx_per_gene": mean(tx_per_gene),
        "mean_exons": mean(exon_counts),
        "mean_cds_length": mean(cds_lengths),
        "mean_genomic_span": mean(genomic_spans),
    }


def write_summary(stats, outdir):
    """gene_set_statistics.txt, laid out as BRAKER4's."""
    path = os.path.join(outdir, "gene_set_statistics.txt")
    n_tx = stats["single_exon_tx"] + stats["multi_exon_tx"]

    def share(count):
        return f"{count / n_tx * 100:.1f}%" if n_tx else "NA"

    with open(path, "w", encoding="utf-8") as f:
        f.write("Paludamentum Gene Set Statistics\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Genes:                  {stats['n_genes']:>8,}\n")
        f.write(f"Transcripts:            {stats['n_transcripts']:>8,}\n")
        f.write(f"Mean transcripts/gene:  {stats['mean_tx_per_gene']:>8.2f}\n")
        f.write(f"Mean exons/transcript:  {stats['mean_exons']:>8.2f}\n")
        f.write("\n")
        f.write(f"Single-exon transcripts:{stats['single_exon_tx']:>8,} ({share(stats['single_exon_tx'])})\n")
        f.write(f"Multi-exon transcripts: {stats['multi_exon_tx']:>8,} ({share(stats['multi_exon_tx'])})\n")
        f.write("\n")
        f.write(f"Mean CDS length:        {stats['mean_cds_length']:>8,.0f} bp\n")
        f.write(f"Mean genomic span:      {stats['mean_genomic_span']:>8,.0f} bp\n")
        if stats["cds_lengths"]:
            f.write(f"Median CDS length:      {upper_median(stats['cds_lengths']):>8,} bp\n")
        if stats["genomic_spans"]:
            f.write(f"Median genomic span:    {upper_median(stats['genomic_spans']):>8,} bp\n")
    return path


def parse_support_tsv(tsv_path):
    """Rows of gene_support.tsv as dicts ('#' lines skipped, first other line is the header)."""
    data = []
    header = None
    with open(tsv_path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n\r")
            if not line.strip() or line.startswith("#"):
                continue
            if header is None:
                header = line.split("\t")
                continue
            cols = line.split("\t")
            if len(cols) >= len(header):
                data.append(dict(zip(header, cols)))
    return data


def support_fractions(support_data):
    """Fraction of introns supported by any evidence, per multi-exon transcript."""
    if not support_data:
        return []
    col_names = list(support_data[0].keys())
    intron_cols = [c for c in col_names if "intron" in c.lower() and ("support" in c.lower() or "sup" in c.lower())]
    sup_col = next((c for c in intron_cols if "any" in c.lower()), intron_cols[0] if intron_cols else None)
    n_col = next((c for c in col_names if c in ("n_introns", "num_introns")), None)
    if sup_col is None or n_col is None:
        raise ValueError("no intron support columns (num_introns, introns_sup_any) in the support table")
    fractions = []
    for row in support_data:
        n = int(row[n_col])
        if n > 0:
            fractions.append(int(row[sup_col]) / n)
    return fractions


def generate_plots(stats, outdir, support_data=None):
    """All plots as PNG (plot_style.py: page width, 200 dpi, common font sizes); {name: path}."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FuncFormatter, LogLocator, MaxNLocator
        import plot_style as ps
    except ImportError:
        print("gene_set_statistics.py: WARNING: matplotlib not available, skipping plots", file=sys.stderr)
        return {}

    plots = {}
    ps.apply(plt)
    pair = {"figsize": (ps.FIG_WIDTH, 3.6), "gridspec_kw": {"width_ratios": [1.6, 1]}}
    thousands = FuncFormatter(lambda x, _pos: f"{x:,.0f}")

    def counts_axis(axis):
        """whole-number ticks with thousands separators (counts)"""
        axis.set_major_locator(MaxNLocator(integer=True, steps=[1, 2, 5, 10]))
        axis.set_major_formatter(thousands)

    def short_bp(x, _pos=None):
        """100, 200, 500, 1k, 2k, 10k, 1M: compact labels of a log axis in bp"""
        for div, suffix in ((1e6, "M"), (1e3, "k")):
            if x >= div:
                return f"{x / div:g}{suffix}"
        return f"{x:g}"

    def no_data(ax, title):
        ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes, color=ps.GREY)
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])

    def save(fig, name):
        fig.tight_layout()
        path = os.path.join(outdir, f"{name}.png")
        fig.savefig(path, dpi=ps.DPI)
        plt.close(fig)
        plots[name] = path

    def pct(part, total):
        p = 100 * part / total
        return "<0.1%" if 0 < p < 0.05 else f"{p:.1f}%"

    def pie(ax, sizes, names, colours):
        total = sum(sizes)
        kept = [(s, n, c) for s, n, c in zip(sizes, names, colours) if s > 0]
        ax.pie([k[0] for k in kept], labels=[f"{n}\n{s:,} ({pct(s, total)})" for s, n, _ in kept],
               colors=[k[2] for k in kept], startangle=90, counterclock=False, labeldistance=1.12,
               wedgeprops={"edgecolor": "white", "linewidth": 1.5}, textprops={"fontsize": 10})
        ax.set_aspect("equal")

    def count_bars(ax, values, cap, colour, annotate):
        """Bars of the integer values 0/1..cap; a value above cap is counted in the bar 'cap', labelled '>=cap'."""
        counts = Counter(min(v, cap) for v in values)
        lo, hi = min(counts), max(counts)
        xs = list(range(lo, hi + 1))
        heights = [counts.get(x, 0) for x in xs]
        bars = ax.bar(xs, heights, width=0.8, color=colour, edgecolor="white")
        if annotate:
            total = len(values)
            for bar, h in zip(bars, heights):
                if h > 0:
                    ax.annotate(f"{h:,}\n{pct(h, total)}", (bar.get_x() + bar.get_width() / 2, h),
                                xytext=(0, 3), textcoords="offset points", ha="center", va="bottom",
                                fontsize=ps.ANNOTATION_SIZE, color=ps.MUTED)
            ax.set_ylim(0, max(heights) * 1.25)
        step = 1 if len(xs) <= 16 else 2
        ticks = [x for x in xs if (x - lo) % step == 0 or x == hi]
        if len(ticks) > 1 and ticks[-1] - ticks[-2] < step:
            ticks.pop(-2)
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"\u2265{t}" if t == cap and max(values) > cap else str(t) for t in ticks])
        counts_axis(ax.yaxis)

    def median_line(ax, median, label):
        ax.axvline(median, color=ps.MEDIAN, linestyle="--", linewidth=1, label=label)
        ax.legend(loc="upper right")

    # --- 1+2. Isoform distribution + exon structure (side by side) ---
    fig, (ax1, ax2) = plt.subplots(1, 2, **pair)
    if stats["tx_per_gene"]:
        count_bars(ax1, stats["tx_per_gene"], 15, ps.BLUE, annotate=True)
        ax1.set_xlabel("Transcripts per gene")
        ax1.set_ylabel("Number of genes")
        ax1.set_title("Transcripts per gene")
    else:
        no_data(ax1, "Transcripts per gene")
    sizes = [stats["multi_exon_tx"], stats["single_exon_tx"]]
    if sum(sizes) > 0:
        pie(ax2, sizes, ["Multi-exon", "Single-exon"], [ps.BLUE, ps.ORANGE])
        ax2.set_title(f"Exon structure ({sum(sizes):,} transcripts)")
    else:
        no_data(ax2, "Exon structure")
    save(fig, "isoform_and_exon_structure")

    # --- 3. Transcript length distributions (log-scale histograms) ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(ps.FIG_WIDTH, 3.6))
    for ax, values, xlabel, title in (
            (ax1, stats["cds_lengths"], "CDS length (bp, log scale)", "CDS length"),
            (ax2, stats["genomic_spans"], "Genomic span (bp, log scale)", "Genomic span (incl. introns)")):
        values = [v for v in values if v > 0]
        if values:
            lo, hi = min(values), max(values)
            hi = max(hi, lo * 1.01)
            edges = [lo * (hi / lo) ** (i / 50) for i in range(51)]
            ax.hist(values, bins=edges, color=ps.BLUE, edgecolor="white", linewidth=0.5)
            ax.set_xscale("log")
            # labels at 1, 2 and 5 of each decade (1 and 3 over more than 3 decades, only 1 over
            # more than 5), unlabelled minor ticks between
            decades = math.log10(hi / lo)
            subs = (1.0, 2.0, 5.0) if decades <= 3 else (1.0, 3.0) if decades <= 5 else (1.0,)
            ax.xaxis.set_major_locator(LogLocator(base=10, subs=subs))
            ax.xaxis.set_major_formatter(FuncFormatter(short_bp))
            ax.xaxis.set_minor_locator(LogLocator(base=10, subs=[k for k in range(1, 10) if k not in subs]))
            ax.xaxis.set_minor_formatter(FuncFormatter(lambda _x, _pos: ""))
            counts_axis(ax.yaxis)
            ax.set_xlabel(xlabel)
            ax.set_ylabel("Number of transcripts")
            ax.set_title(title)
            median = upper_median(values)
            median_line(ax, median, f"median {median:,} bp")
        else:
            no_data(ax, title)
    save(fig, "transcript_lengths")

    # --- 4. Introns per gene ---
    fig, ax = plt.subplots(figsize=(ps.FIG_WIDTH, 3.2))
    if stats["introns_per_gene"]:
        count_bars(ax, stats["introns_per_gene"], 30, ps.BLUE, annotate=False)
        ax.set_xlabel("Introns per gene")
        ax.set_ylabel("Number of genes")
        ax.set_title("Introns per gene")
        median = upper_median(stats["introns_per_gene"])
        median_line(ax, median, f"median {median}")
    else:
        no_data(ax, "Introns per gene")
    save(fig, "introns_per_gene")

    # --- 5. Evidence support (intron support of the multi-exon transcripts) ---
    if support_data is not None:
        fractions = support_fractions(support_data)
        fig, axes = plt.subplots(1, 2, **pair)
        if fractions:
            axes[0].hist(fractions, bins=[i / 20 for i in range(21)], color=ps.BLUE, edgecolor="white")
            axes[0].set_xlabel("Fraction of introns supported")
            axes[0].set_ylabel("Number of transcripts")
            axes[0].set_title("Intron evidence support")
            axes[0].set_xlim(-0.02, 1.02)
            counts_axis(axes[0].yaxis)
            full = sum(1 for f in fractions if f >= 1.0)
            partial = sum(1 for f in fractions if 0 < f < 1.0)
            none = sum(1 for f in fractions if f == 0)
            pie(axes[1], [full, partial, none], ["Full support", "Partial", "No support"],
                [ps.GREEN, ps.ORANGE, ps.RED])
            axes[1].set_title(f"Evidence support ({len(fractions):,} multi-exon transcripts)")
        else:
            no_data(axes[0], "Intron evidence support")
            no_data(axes[1], "Evidence support (no multi-exon transcripts)")
        save(fig, "evidence_support")

    return plots


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff3", required=True, help="gene set (GFF3)")
    parser.add_argument("--out-dir", required=True, help="output directory (created if missing)")
    parser.add_argument("--support", default=None, help="gene_support.tsv of gene_support.py")
    args = parser.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    try:
        ann = read_gff3(args.gff3)
    except GFF3Error as exc:
        print(f"gene_set_statistics.py: {args.gff3}: {exc}", file=sys.stderr)
        return 1
    stats = compute_statistics(ann)
    summary_path = write_summary(stats, args.out_dir)

    support_data = None
    if args.support:
        try:
            support_data = parse_support_tsv(args.support)
            support_fractions(support_data)       # fail early on a table without the support columns
        except (OSError, ValueError) as exc:
            print(f"gene_set_statistics.py: {args.support}: {exc}", file=sys.stderr)
            return 1

    plots = generate_plots(stats, args.out_dir, support_data)
    print(f"gene_set_statistics.py: {stats['n_genes']} genes, {stats['n_transcripts']} transcripts; "
          f"{summary_path}" + "".join(f", {p}" for p in plots.values()), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
