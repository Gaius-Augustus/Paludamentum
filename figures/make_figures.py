#!/usr/bin/env python3
"""Writes figures/overview.svg and figures/drusilla_flow.svg."""
import sys
from xml.sax.saxutils import escape

SANS = "Helvetica, Arial, sans-serif"
MONO = "ui-monospace, SFMono-Regular, Menlo, Consolas, 'DejaVu Sans Mono', monospace"
INK, SUB = "#1f2933", "#52606d"
COL = {  # fill, stroke
    "green":  ("#e3f1e0", "#4f9440"),
    "blue":   ("#dfeaf8", "#4577ba"),
    "red":    ("#f9dedc", "#c0504d"),
    "yellow": ("#fdf1cc", "#b8921a"),
    "orange": ("#fde6cf", "#cf7f1f"),
    "purple": ("#ebe1f2", "#8a63a8"),
    "grey":   ("#eceff3", "#5b6675"),
}


class Fig:
    def __init__(self, w, h, title, desc):
        self.w, self.h, self.title, self.desc = w, h, title, desc
        self.body = []

    def add(self, s):
        self.body.append(s)

    def text(self, x, y, s, size=12.5, weight="normal", fill=SUB, anchor="start", mono=False, style=""):
        fam = MONO if mono else SANS
        self.add(f'<text x="{x}" y="{y}" font-family="{fam}" font-size="{size}" font-weight="{weight}" '
                 f'fill="{fill}" text-anchor="{anchor}"{style}>{escape(s)}</text>')

    def box(self, x, y, w, h, color, title, lines, num=None):
        fill, stroke = COL[color]
        self.add(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="10" fill="{fill}" stroke="{stroke}" stroke-width="1.6"/>')
        tx = x + 14
        if num is not None:
            self.add(f'<circle cx="{x + 24}" cy="{y + 25}" r="11" fill="{stroke}"/>')
            self.text(x + 24, y + 29.5, str(num), 13, "bold", "#ffffff", "middle")
            tx = x + 43
        self.text(tx, y + 30, title, 15, "bold", INK)
        for i, ln in enumerate(lines):
            self.text(x + 14, y + 53 + 17 * i, ln)

    def doc(self, x, y, w, h, color, title, sub, mono=False, tag=None):
        fill, stroke = COL[color]
        b = y + h - 7
        d = (f'M{x},{y + 6} q0,-6 6,-6 h{w - 12} q6,0 6,6 V{b} '
             f'q{-w / 4},-15 {-w / 2},0 t{-w / 2},0 Z')
        self.add(f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="1.6" stroke-linejoin="round"/>')
        self.text(x + w / 2, y + 31, title, 15, "bold", INK, "middle")
        self.text(x + w / 2, y + 52, sub, 11.5 if mono else 12.5, anchor="middle", mono=mono)
        if tag:
            self.text(x + w / 2, y + 69, tag, 11, fill=stroke, anchor="middle", weight="bold")

    def arrow(self, pts, color, dashed=False, head=True):
        stroke = COL[color][1]
        d = "M" + " L".join(f"{x},{y}" for x, y in pts)
        dash = ' stroke-dasharray="6 4"' if dashed else ""
        end = f' marker-end="url(#a-{color})"' if head else ""
        self.add(f'<path d="{d}" fill="none" stroke="{stroke}" stroke-width="2" stroke-linejoin="round"{dash}{end}/>')

    def dot(self, x, y, color):
        self.add(f'<circle cx="{x}" cy="{y}" r="3.5" fill="{COL[color][1]}"/>')

    def svg(self):
        markers = "".join(
            f'<marker id="a-{k}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
            f'orient="auto-start-reverse"><path d="M0,0.8 L10,5 L0,9.2 Z" fill="{v[1]}"/></marker>'
            for k, v in COL.items())
        return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {self.w} {self.h}" width="{self.w}" '
                f'height="{self.h}" role="img" aria-labelledby="t d">\n'
                f'<title id="t">{escape(self.title)}</title>\n<desc id="d">{escape(self.desc)}</desc>\n'
                f'<defs>{markers}</defs>\n'
                f'<rect x="1" y="1" width="{self.w - 2}" height="{self.h - 2}" rx="14" fill="#ffffff" stroke="#d9dee4"/>\n'
                + "\n".join(self.body) + "\n</svg>\n")


def overview():
    f = Fig(1160, 606, "Paludamentum overview",
            "Inputs (genome, proteins, short reads, Iso-Seq reads), the five processing steps "
            "(ab initio prediction, protein evidence, transcript evidence, high-confidence genes, "
            "integration) and the output files of Paludamentum.")
    A, B, C, D = 30, 300, 620, 920
    wA, wB, wC, wD = 210, 260, 240, 210
    h, r = 88, [80, 198, 316, 434]
    for x, w, s in ((A, wA, "INPUTS"), (B, wB + C + wC - B - wB, "PROCESSING"), (D, wD, "OUTPUTS")):
        f.text(x, 52, s, 13, "bold", SUB, style=' letter-spacing="1.5"')
        f.add(f'<path d="M{x},62 h{w}" stroke="#d9dee4" stroke-width="1.5"/>')

    f.doc(A, r[0], wA, h, "green", "Genome", "FASTA", tag="required")
    f.doc(A, r[1], wA, h, "blue", "Proteins", "FASTA and/or OrthoDB v12", tag="optional")
    f.doc(A, r[2], wA, h, "red", "Short reads", "FASTQ, BAM, SRA or VARUS", tag="optional")
    f.doc(A, r[3], wA, h, "yellow", "Iso-Seq reads", "FASTQ, SRA or VARUS", tag="optional")

    f.box(B, r[0], wB, h, "green", "Ab initio prediction",
          ["Tiberius or Vipsania,", "genome chunks in parallel on GPUs"], 1)
    f.box(B, r[1], wB, h, "blue", "Protein evidence",
          ["miniprot, miniprot-boundary-scorer,", "miniprothint"], 2)
    f.box(B, r[2], wB, r[3] + h - r[2], "orange", "Transcript evidence",
          ["HISAT2 for short reads,", "minimap2 for Iso-Seq reads", "",
           "Libraries with a low alignment", "rate are dropped.", "",
           "StringTie assembles transcripts."], 3)

    f.box(C, r[0], wC, h, "grey", "Integration",
          ["HC genes are merged with the", "ab initio predictions"], 5)
    f.box(C, r[1], wC, r[3] + h - r[1], "purple", "High-confidence genes",
          ["HC genes: ORFs of the assembled", "transcripts that the protein", "evidence supports", "",
           "TransDecoder or TD2,", "DIAMOND", "",
           "Vertebrate models:", "Drusilla ORFs and a filter for the", "ab initio predictions", "",
           "Without reads: from the protein", "alignments alone"], 4)

    f.doc(D, r[0], wD, h, "grey", "Final annotation", "<tool>_evidence.gff3", mono=True)
    f.doc(D, r[1], wD, h, "grey", "Protein sequences", "<tool>_evidence_proteins.fa", mono=True)
    y0 = r[2]
    f.add(f'<rect x="{D}" y="{y0}" width="{wD}" height="{r[3] + h - y0}" rx="10" fill="#fafbfc" '
          f'stroke="#b6bec8" stroke-width="1.4" stroke-dasharray="5 4"/>')
    f.text(D + 14, y0 + 28, "Also written", 14, "bold", INK)
    for i, (name, what) in enumerate((("<tool>_ab_initio.gff3", "ab initio predictions"),
                                      ("intermediate/hc.gff3", "HC genes"),
                                      ("hintsfile.gff", "protein and transcript hints"))):
        f.text(D + 14, y0 + 58 + 48 * i, name, 11.5, fill=INK, mono=True)
        f.text(D + 14, y0 + 75 + 48 * i, what)

    for i, c in enumerate(("green", "blue", "red", "yellow")):
        y = r[i] + h / 2
        f.arrow([(A + wA, y), (B, y)], c)
    f.arrow([(B + wB, r[0] + h / 2), (C, r[0] + h / 2)], "green")
    f.arrow([(B + wB, r[1] + h / 2), (C, r[1] + h / 2)], "blue")
    ym = (r[2] + r[3] + h) / 2
    f.arrow([(B + wB, ym), (C, ym)], "orange")
    f.arrow([(C + wC / 2, r[1]), (C + wC / 2, r[0] + h)], "purple")
    f.arrow([(C + wC, r[0] + h / 2), (D, r[0] + h / 2)], "grey")
    f.arrow([(D + wD / 2, r[0] + h - 4), (D + wD / 2, r[1])], "grey")

    f.text(A, 556, "Genome only: step 1 runs alone, for example to parallelize a gene finder over several GPUs.")
    f.text(A, 577, "Reads need proteins: HC genes are built from transcript and protein evidence together. "
                   "<tool> is tiberius or vipsania.")
    return f


def drusilla():
    f = Fig(1160, 510, "Drusilla flow of Paludamentum",
            "For vertebrate models: one StringTie assembly, a transcript filter, Drusilla ORFs and a "
            "codon fix give the high-confidence genes; a LightGBM filter and a hint rescue treat the "
            "ab initio predictions; both are merged into the final annotation.")
    w, h = 160, 88
    X = [20, 202, 384, 566, 748]
    M, wM = 968, 168
    ya, yb, yc = 44, 190, 336

    f.doc(X[0], ya, w, h, "red", "Read alignments", "short reads, Iso-Seq")
    f.box(X[1], ya, w, h, "orange", "StringTie", ["one assembly of", "all reads"], 1)
    f.box(X[2], ya, w, h, "orange", "Transcript filter", ["length, coverage,", "TPM"], 2)
    f.box(X[3], ya, w, h, "purple", "Drusilla", ["ORFs, also truncated", "at transcript ends"], 3)
    f.box(X[4], ya, w, h, "purple", "Codon fix", ["stop and start codons;", "isoforms collapsed"], 4)

    f.box(X[3], yb, w, h, "blue", "Protein evidence", ["miniprot alignments,", "miniprothint hints"])

    f.doc(X[1], yc, w, h, "green", "Ab initio genes", "Tiberius or Vipsania")
    f.box(X[2], yc, w, h, "green", "LightGBM filter", ["wrong, partial", "or correct"], 5)
    f.box(X[4], yc, w, h, "green", "Hint rescue", ["Tiberius with", "protein hints"], 6)

    f.box(M, yb, wM, h, "grey", "Merge", ["ORFs, kept and rescued", "ab initio genes"], 7)
    f.doc(M, yc, wM, h, "grey", "Final annotation", "<tool>_evidence.gff3", mono=True)

    cy = lambda y: y + h / 2
    for i, c in enumerate(("red", "orange", "orange", "purple")):
        f.arrow([(X[i] + w, cy(ya)), (X[i + 1], cy(ya))], c)
    # Drusilla ORFs to the merge
    xm = M + wM / 2
    f.arrow([(X[4] + w, cy(ya)), (xm, cy(ya)), (xm, yb)], "purple")
    f.text((X[4] + w + xm) / 2 + 4, cy(ya) - 9, "HC genes", anchor="middle")
    # protein evidence to codon fix, hint rescue and filter
    xt = X[4] + w / 2
    f.arrow([(X[3] + w, cy(yb)), (xt, cy(yb)), (xt, ya + h)], "blue")
    f.arrow([(xt, cy(yb)), (xt, yc)], "blue")
    f.dot(xt, cy(yb), "blue")
    xl = X[2] + w / 2
    f.arrow([(X[3], cy(yb)), (xl, cy(yb)), (xl, yc)], "blue")
    # ab initio lane
    f.arrow([(X[1] + w, cy(yc)), (X[2], cy(yc))], "green")
    f.arrow([(X[2] + w, cy(yc)), (X[4], cy(yc))], "green")
    f.text((X[2] + w + X[4]) / 2, cy(yc) - 26, "partial, no transcript", anchor="middle")
    f.text((X[2] + w + X[4]) / 2, cy(yc) - 9, "kept at the locus", anchor="middle")
    xj, yk = 938, yc + h + 30
    f.arrow([(X[4] + w, cy(yc)), (xj, cy(yc))], "green", head=False)
    f.arrow([(xl, yc + h), (xl, yk), (xj, yk), (xj, cy(yb)), (M, cy(yb))], "green")
    f.dot(xj, cy(yc), "green")
    f.text((xl + xj) / 2, yk + 19, "correct: kept", anchor="middle")
    f.arrow([(xm, yb + h), (xm, yc)], "grey")
    return f


out = sys.argv[1]
for name, fig in (("overview", overview()), ("drusilla_flow", drusilla())):
    with open(f"{out}/{name}.svg", "w") as fh:
        fh.write(fig.svg())
