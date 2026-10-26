#!/usr/bin/env python3
"""Writes figures/overview.svg, figures/combined_flow.svg, figures/insect_flow.svg and
figures/postprocess_flow.svg."""
import re
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


def italic(s):
    """ab initio is always set in italics."""
    return re.sub(r"\b([Aa]b initio)\b", r'<tspan font-style="italic">\1</tspan>', s)


class Fig:
    def __init__(self, w, h, title, desc):
        self.w, self.h, self.title, self.desc = w, h, title, desc
        self.body = []

    def add(self, s):
        self.body.append(s)

    def text(self, x, y, s, size=12.5, weight="normal", fill=SUB, anchor="start", mono=False, style=""):
        fam = MONO if mono else SANS
        self.add(f'<text x="{x}" y="{y}" font-family="{fam}" font-size="{size}" font-weight="{weight}" '
                 f'fill="{fill}" text-anchor="{anchor}"{style}>{italic(escape(s))}</text>')

    def box(self, x, y, w, h, color, title, lines, num=None, dashed=False, fill=None):
        fill, stroke = (fill or COL[color][0]), COL[color][1]
        dash = ' stroke-dasharray="6 4"' if dashed else ""
        self.add(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="10" fill="{fill}" stroke="{stroke}" stroke-width="1.6"{dash}/>')
        tx = x + 14
        if num is not None:
            self.add(f'<circle cx="{x + 24}" cy="{y + 25}" r="11" fill="{stroke}"/>')
            self.text(x + 24, y + 29.5, str(num), 13, "bold", "#ffffff", "middle")
            tx = x + 43
        self.text(tx, y + 30, title, 15, "bold", INK)
        for i, ln in enumerate(lines):
            self.text(x + 14, y + 53 + 17 * i, ln)

    def doc(self, x, y, w, h, color, title, sub, mono=False, tag=None, dashed=False, fill=None):
        fill, stroke = (fill or COL[color][0]), COL[color][1]
        b = y + h - 7
        d = (f'M{x},{y + 6} q0,-6 6,-6 h{w - 12} q6,0 6,6 V{b} '
             f'q{-w / 4},-15 {-w / 2},0 t{-w / 2},0 Z')
        dash = ' stroke-dasharray="6 4"' if dashed else ""
        self.add(f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="1.6" stroke-linejoin="round"{dash}/>')
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
    """Five-step overview of the unified pipeline."""
    f = Fig(1400, 606, "Paludamentum overview",
            "Inputs (genome, proteins, short reads, Iso-Seq reads), the five processing steps "
            "(ab initio predictions with both gene finders, protein evidence, transcript evidence, "
            "combined high-confidence genes with combine_genesets, post-processing and quality control) "
            "and the output files of Paludamentum.")
    A, B, C, E, D = 30, 290, 580, 870, 1170
    wA, wB, wC, wE, wD = 220, 260, 270, 280, 200
    h, r = 88, [80, 198, 316, 434]
    for x, w, s in ((A, wA, "INPUTS"), (B, E + wE - B, "PROCESSING"), (D, wD, "OUTPUTS")):
        f.text(x, 52, s, 13, "bold", SUB, style=' letter-spacing="1.5"')
        f.add(f'<path d="M{x},62 h{w}" stroke="#d9dee4" stroke-width="1.5"/>')

    f.doc(A, r[0], wA, h, "green", "Genome", "FASTA", tag="required")
    f.doc(A, r[1], wA, h, "blue", "Proteins", "FASTA and/or OrthoDB v12",
          tag="required with reads")
    f.doc(A, r[2], wA, h, "red", "Short reads", "FASTQ, BAM, SRA or pyVARUS", tag="optional")
    f.doc(A, r[3], wA, h, "yellow", "Iso-Seq reads", "FASTQ, SRA or pyVARUS", tag="optional")

    f.box(B, r[0], wB, h, "green", "Ab initio predictions",
          ["Tiberius and Vipsania", "(Vipsania --finetune) in parallel on GPUs"], 1)
    f.box(B, r[1], wB, h, "blue", "Protein evidence",
          ["miniprot, miniprot-boundary-scorer,", "miniprothint"], 2)
    f.box(B, r[2], wB, r[3] + h - r[2], "orange", "Transcript evidence",
          ["HISAT2 for short reads,", "minimap2 for Iso-Seq reads", "",
           "Libraries with a low alignment", "rate are dropped.", "",
           "StringTie: one assembly", "of all reads"], 3)

    f.box(C, r[0], wC, r[3] + h - r[0], "purple", "Combined high-confidence genes",
          ["TD2 (shorten + DIAMOND +", "revise) and Drusilla predict",
           "ORFs on the one StringTie", "assembly.", "",
           "DIAMOND runs on all four", "protein sets (Tiberius, Vipsania,",
           "revised TD2, Drusilla),", "bedtools genomecov writes the",
           "depth, bam2hints the introns.", "",
           "combine_genesets.py", "integrates the four gene sets",
           "with the RNA and protein", "evidence into the final",
           "annotation. Insect clades", "add --drop-lone-single-exon",
           "tiberius."], 4)

    f.box(E, r[0], wE, r[3] + h - r[0], "grey", "Post-processing",
          ["and quality control", "",
           "Sanity filter: broken CDS", "and internal stops removed", "",
           "UTRs from the StringTie", "assembly", "",
           "GFF3 (validated), GTF,", "protein and CDS FASTA", "",
           "Optional: BUSCO, compleasm,", "OMArk, gffcompare, ncRNA,", "GO terms (FANTASIA-Lite)", "",
           "report.html sums it up"], 5)

    f.doc(D, r[0], wD, h, "grey", "Final annotation", "combined_evidence.gff3", mono=True)
    f.doc(D, r[1], wD, h, "grey", "Sequences, GTF", "combined_evidence_*.fa", mono=True)
    y0 = r[2]
    f.add(f'<rect x="{D}" y="{y0}" width="{wD}" height="{r[3] + h - y0}" rx="10" fill="#fafbfc" '
          f'stroke="#b6bec8" stroke-width="1.4" stroke-dasharray="5 4"/>')
    f.text(D + 14, y0 + 28, "Also written", 14, "bold", INK)
    for i, (name, what) in enumerate((("intermediate/", "the four gene sets,"),
                                      ("", "combine_genesets report"),
                                      ("qc/, report.html", "quality control"),
                                      ("hintsfile.gff", "all hints"))):
        f.text(D + 14, y0 + 58 + 38 * i, name, 11.5, fill=INK, mono=True)
        f.text(D + 14, y0 + 75 + 38 * i, what)

    # input -> processing arrows
    for i, c in enumerate(("green", "blue", "red", "yellow")):
        y = r[i] + h / 2
        f.arrow([(A + wA, y), (B, y)], c)
    # processing -> combined HC
    f.arrow([(B + wB, r[0] + h / 2), (C, r[0] + h / 2)], "green")
    f.arrow([(B + wB, r[1] + h / 2), (C, r[1] + h / 2)], "blue")
    ym = (r[2] + r[3] + h) / 2
    f.arrow([(B + wB, ym), (C, ym)], "orange")
    # combined HC -> post-processing (one arrow, middle)
    yc = (r[0] + r[3] + h) / 2
    f.arrow([(C + wC, yc), (E, yc)], "purple")
    # post-processing -> outputs (fan to the four output rows)
    f.arrow([(E + wE, r[0] + h / 2), (D, r[0] + h / 2)], "grey")
    f.arrow([(E + wE, r[1] + h / 2), (D, r[1] + h / 2)], "grey")
    f.arrow([(E + wE, (r[2] + r[3] + h) / 2), (D, (r[2] + r[3] + h) / 2)], "grey")

    f.text(A, r[3] + h + 40,
           "Modes abinitio and proteins run step 1 (one gene finder) and step 5 only; their output stem is <tool>_ab_initio / <tool>_evidence.")
    return f


def combined_flow():
    """Zoom into step 4 (combined HC genes). Replaces the earlier drusilla_flow figure."""
    f = Fig(1200, 540, "Combined HC-gene step of Paludamentum",
            "The unified HC-gene step of Paludamentum: Tiberius and Vipsania predict genes on the genome; "
            "StringTie assembles one transcript set; TD2 (with shorten + classify + revise) and Drusilla "
            "predict ORFs on it; DIAMOND runs on all four protein sets; bedtools genomecov writes the "
            "coverage, bam2hints the introns; combine_genesets.py integrates everything into "
            "combined.gff3.")

    w, h = 170, 90
    ya, yb, yc, yd = 54, 184, 314, 410       # four rows: ab initio, StringTie, ORF callers, DIAMOND
    cols = [20, 212, 404, 596, 788, 980]
    cy = lambda y: y + h / 2

    # inputs on the left
    f.doc(cols[0], ya, w, h, "green", "Genome", "FASTA")
    f.doc(cols[0], yb, w, h, "red", "Reads", "RNA-Seq / Iso-Seq", tag="required with proteins")
    f.doc(cols[0], yc, w, h, "blue", "Proteins", "FASTA / OrthoDB v12", tag="required")

    # Row 1 — two gene finders on the genome
    f.box(cols[1], ya, w, h, "green", "Tiberius", ["ab initio prediction", "on the genome"], 1)
    f.box(cols[2], ya, w, h, "green", "Vipsania", ["--finetune then", "annotate the genome"], 1)

    # Row 2 — StringTie
    f.box(cols[1], yb, cols[2] + w - cols[1], h, "orange", "StringTie",
          ["HISAT2 + minimap2 align the reads;", "one assembly of all reads (StringTie, or --mix)"], 2)

    # Row 3 — ORF callers on the StringTie assembly
    f.box(cols[3], yb, w, h, "purple", "TD2 + revise",
          ["LongOrfs, Predict,", "shorten, revise with", "two DIAMOND runs"], 3)
    f.box(cols[4], yb, w, h, "purple", "Drusilla",
          ["filter StringTie,", "predict ORFs,", "collapse isoforms"], 3)

    # Row 4 — four DIAMOND runs
    f.box(cols[1], yc, w, h, "blue", "DIAMOND",
          ["tiberius.pep", "vs protein DB", "(14-col outfmt)"])
    f.box(cols[2], yc, w, h, "blue", "DIAMOND",
          ["vipsania.pep", "vs protein DB", "(14-col outfmt)"])
    f.box(cols[3], yc, w, h, "blue", "DIAMOND",
          ["td2.pep (revised)", "vs protein DB", "(14-col outfmt)"])
    f.box(cols[4], yc, w, h, "blue", "DIAMOND",
          ["drusilla.pep", "vs protein DB", "(14-col outfmt)"])

    # miniprothint / BAM-side evidence
    f.box(cols[3], ya, w, h, "blue", "miniprot +",
          ["miniprothint:", "scored alignments,", "hc.gff"])
    f.box(cols[4], ya, w, h, "red", "bedtools / bam2hints",
          ["rna_cov.bg.gz,", "introns.gff from the", "merged BAM"])

    # combine_genesets block
    f.box(cols[5], yb, w, h, "grey", "combine_genesets.py",
          ["4 gene sets +", "4 DIAMOND +", "miniprot/hints +", "coverage +", "introns"], 4)
    f.doc(cols[5], yd - 16, w, h, "grey", "combined.gff3", "final annotation", mono=True)

    # arrows
    for i in range(2):
        f.arrow([(cols[0] + w, cy(ya)), (cols[1], cy(ya))], "green", head=i == 0)
    f.arrow([(cols[1] + w, cy(ya)), (cols[2], cy(ya))], "green", head=False)
    f.arrow([(cols[0] + w, cy(yb)), (cols[1], cy(yb))], "red")
    # StringTie -> both ORF callers
    f.arrow([(cols[2] + w, cy(yb)), (cols[3], cy(yb))], "orange")
    f.arrow([(cols[3] + w, cy(yb)), (cols[4], cy(yb))], "orange", head=False)
    # ORF callers -> combine_genesets
    f.arrow([(cols[4] + w, cy(yb)), (cols[5], cy(yb))], "purple")
    # gene finders -> their DIAMOND runs (down)
    f.arrow([(cols[1] + w / 2, ya + h), (cols[1] + w / 2, yc)], "green")
    f.arrow([(cols[2] + w / 2, ya + h), (cols[2] + w / 2, yc)], "green")
    # ORF callers -> their DIAMOND runs (down)
    f.arrow([(cols[3] + w / 2, yb + h), (cols[3] + w / 2, yc)], "purple")
    f.arrow([(cols[4] + w / 2, yb + h), (cols[4] + w / 2, yc)], "purple")
    # proteins -> DIAMOND, miniprot (bus down the left)
    xbus = cols[0] + w + 20
    f.arrow([(cols[0] + w, cy(yc)), (xbus, cy(yc)), (xbus, cy(ya)), (cols[3], cy(ya))], "blue")
    f.dot(xbus, cy(yc), "blue")
    for x in (cols[1], cols[2], cols[3], cols[4]):
        f.arrow([(xbus, cy(yc) + 10), (x + w / 2, cy(yc) + 10), (x + w / 2, yc)], "blue", head=False)
    # DIAMOND bus into combine_genesets
    xdbus = cols[5] - 20
    f.arrow([(cols[4] + w, cy(yc)), (xdbus, cy(yc)), (xdbus, cy(yb) + 20), (cols[5], cy(yb) + 20)], "blue")
    # bam2hints / miniprothint into combine
    f.arrow([(cols[4] + w, cy(ya)), (xdbus, cy(ya)), (xdbus, cy(yb) - 20), (cols[5], cy(yb) - 20)], "red")
    f.arrow([(cols[3] + w / 2, ya + h / 2 - 2), (cols[3] + w / 2, ya + h / 2 - 2)], "blue", head=False)
    # combine -> output
    f.arrow([(cols[5] + w / 2, yb + h), (cols[5] + w / 2, yd - 16)], "grey")

    # legend bottom
    f.text(20, yd + h + 10, "Numbers correspond to the overview figure (step 1 ab initio, step 2 proteins, step 3 transcripts, step 4 combine).")
    f.text(20, yd + h + 30, "Insect clades (Insecta, Lepidoptera) add --drop-lone-single-exon tiberius to combine_genesets.py.")
    return f


def insect():
    """The unified HC flow applied to insect runs: --drop-lone-single-exon tiberius is passed."""
    f = Fig(1310, 580, "Insect run of Paludamentum",
            "The unified HC-gene step applied to an insect gene-finder model: Tiberius and Vipsania both run "
            "(Vipsania with --finetune on the insect model), both ORF callers run on the StringTie assembly, "
            "and combine_genesets receives --drop-lone-single-exon tiberius for insect clades.")
    A, wA = 30, 200                     # inputs
    X, w = [258, 440, 622, 804], 164    # processing lane
    M, wM = 1000, 162                   # combine
    O, wO = 30 + 1310 - 230, 180        # not used, kept for layout spacing

    ya, yb, yc, h, hc = 80, 220, 362, 108, 118
    cy = lambda y: y + h / 2

    for x, wid, label in ((A, wA, "INPUTS"), (X[0], M + wM - X[0], "PROCESSING"),
                          (M + wM + 20, 128, "OUTPUT")):
        f.text(x, 52, label, 13, "bold", SUB, style=' letter-spacing="1.5"')
        f.add(f'<path d="M{x},62 h{wid}" stroke="#d9dee4" stroke-width="1.5"/>')

    f.doc(A, ya, wA, h, "red", "RNA-Seq or Iso-Seq", "FASTQ, BAM, SRA, pyVARUS", tag="required")
    f.doc(A, yb, wA, h, "green", "Genome", "FASTA", tag="required")
    f.doc(A, yc, wA, h, "blue", "Proteins", "FASTA / OrthoDB v12", tag="required")

    f.box(X[0], ya, w, h, "orange", "Alignment",
          ["HISAT2, minimap2;", "BAM as given"], 3)
    f.box(X[1], ya, w, h, "orange", "StringTie",
          ["one assembly", "of all reads"], 3)
    f.box(X[2], ya, w, h, "purple", "TD2 + revise",
          ["LongOrfs, Predict,", "shorten + classify", "with two DIAMOND"], 4)
    f.box(X[3], ya, w, h, "purple", "Drusilla",
          ["model insects:", "ORFs on the", "kept StringTie"], 4)

    f.box(X[0], yb, X[3] + w - X[0], h, "green", "Both gene finders in parallel",
          ["Tiberius (insect model)   +   Vipsania Insecta, --finetune on the genome",
           "(Vipsania fine-tune is one GPU task on the whole genome)"], 1)

    f.box(X[0], yc, X[3] + w - X[0], h, "blue", "Protein evidence (required)",
          ["miniprot + miniprot-boundary-scorer + miniprothint;",
           "DIAMOND (14-col) on all four protein sets"], 2)

    f.box(M, ya, wM, yc + h - ya, "grey", "combine_genesets.py",
          ["Tiberius + Vipsania +", "revised TD2 + Drusilla,",
           "all 4 DIAMOND tables,", "miniprot / hc.gff hints,",
           "introns + rna_cov.bg.gz;", "",
           "Insecta clade adds", "--drop-lone-single-exon", "tiberius"], 4)

    f.doc(M, yc + h + 20, wM, h, "grey", "Final annotation", "combined_evidence.gff3", mono=True)

    # arrows: reads -> alignment -> stringtie -> TD2/Drusilla -> combine
    for i, c in enumerate(("red", "orange", "orange", "purple")):
        x0 = A + wA if i == 0 else X[i - 1] + w
        f.arrow([(x0, cy(ya)), (X[i], cy(ya))], c)
    f.arrow([(X[3] + w, cy(ya)), (M, cy(ya))], "purple")
    # genome -> gene finders -> combine
    f.arrow([(A + wA, cy(yb)), (X[0], cy(yb))], "green")
    f.arrow([(X[3] + w, cy(yb)), (M, cy(yb))], "green")
    # proteins -> protein evidence -> combine
    f.arrow([(A + wA, cy(yc)), (X[0], cy(yc))], "blue")
    f.arrow([(X[3] + w, cy(yc)), (M, cy(yc))], "blue")
    # combine -> final
    f.arrow([(M + wM / 2, yc + h), (M + wM / 2, yc + h + 20)], "grey")

    f.text(A, yc + h + 50 + h, "The unified HC flow is the same for every clade; insect clades add --drop-lone-single-exon tiberius to combine_genesets.")
    f.text(A, yc + h + 70 + h, "A run that differs from the clade's recommended flow in conf/hc_genes.yaml gets a warning (gene finder, fine-tune, transcripts, proteins).")
    return f


def postprocess():
    f = Fig(1160, 764, "Post-processing and quality control of Paludamentum",
            "The final annotation of step 4 (or the ab initio genes) runs through the sanity filter, "
            "gets UTRs from the StringTie assemblies when there are transcripts, and is written as GFF3, "
            "GTF, proteins and CDS. Gene set statistics always run, evidence support in evidence modes; "
            "completeness, OMArk, gffcompare, ncRNA and GO terms are switched on by options. "
            "report.html sums up all steps.")
    A, wA = 24, 180
    B, wB = 236, 220
    C, wC = 488, 200
    D, wD = 744, 356
    xl, xr = 718, 1126  # bus into and out of the QC steps
    h, r = 88, [80, 196, 312, 428]
    for x, w, s in ((A, wA, "INPUTS"), (B, wB, "POST-PROCESSING"), (C, wC, "OUTPUTS"),
                    (D, wD, "QUALITY CONTROL AND ANNOTATION")):
        f.text(x, 52, s, 13, "bold", SUB, style=' letter-spacing="1.5"')
        f.add(f'<path d="M{x},62 h{w}" stroke="#d9dee4" stroke-width="1.5"/>')

    f.doc(A, r[0], wA, h, "grey", "Annotation", "from step 4", tag="or the ab initio genes")
    f.doc(A, r[1], wA, h, "orange", "StringTie", "assembled transcripts", tag="modes with transcripts")
    f.doc(A, r[2], wA, h, "green", "Genome", "FASTA")

    f.box(B, r[0], wB, h, "grey", "Sanity filter",
          ["internal stop or broken CDS:", "transcript removed"], 1)
    f.box(B, r[1], wB, h, "orange", "UTRs",
          ["from the StringTie assemblies,", "the CDS stays as it is"], 2)
    f.box(B, r[2], wB, h, "grey", "Final files",
          ["GFF3 checked with GenomeTools;", "GTF, proteins, CDS"], 3)

    y0, hy = r[0], r[1] + h - r[0]
    f.add(f'<rect x="{C}" y="{y0}" width="{wC}" height="{hy}" rx="10" fill="#fafbfc" '
          f'stroke="#b6bec8" stroke-width="1.4" stroke-dasharray="5 4"/>')
    f.text(C + 14, y0 + 28, "Also written", 14, "bold", INK)
    for i, (name, what) in enumerate((("<stem>_with_ncRNA.gff3", "plus ncRNA genes"),
                                      ("<stem>_go.gff3", "plus GO terms"),
                                      ("qc/", "tables and plots of the QC"),
                                      ("intermediate/", "before filter and UTRs"))):
        f.text(C + 14, y0 + 54 + 42 * i, name, 11.5, fill=INK, mono=True)
        f.text(C + 14, y0 + 70 + 42 * i, what)
    f.doc(C, r[2], wC, h, "grey", "Final annotation", "<stem>.gff3", mono=True)
    f.doc(C, r[3], wC, h, "grey", "GTF and sequences", "GTF, proteins, CDS")

    # QC and annotation steps; color: what the step needs, pill: what switches it on
    steps = (("grey", "Gene set statistics", "always",
              [("genes, isoforms, exons, introns, lengths", SUB)]),
             ("blue", "Evidence support", "evidence modes",
              [("introns and CDS backed by protein or RNA-Seq hints", SUB)]),
             ("purple", "Completeness", "--busco_lineage",
              [("compleasm and BUSCO on genome and proteome", SUB)]),
             ("purple", "OMArk", "qc.omark",
              [("proteome consistency; needs an OMAmer database", SUB)]),
             ("purple", "gffcompare", "--reference_annotation",
              [("CDS compared with a reference annotation", SUB)]),
             ("purple", "ncRNA", "--ncrna",
              [("rRNA, tRNA, Rfam families; ", SUB),
               ("lncRNA with transcripts", COL["orange"][1])]),
             ("purple", "GO terms", "--fantasia",
              [("FANTASIA-Lite: GO terms by ProtT5 embeddings", SUB)]))
    hs, ys = 64, [80 + 80 * i for i in range(len(steps))]
    for (color, title, pill, parts), y in zip(steps, ys):
        fill, stroke = COL[color]
        f.add(f'<rect x="{D}" y="{y}" width="{wD}" height="{hs}" rx="10" fill="{fill}" stroke="{stroke}" stroke-width="1.6"/>')
        f.text(D + 14, y + 26, title, 15, "bold", INK)
        pw = 7.1 * len(pill) + 16
        f.add(f'<rect x="{D + wD - 10 - pw}" y="{y + 10}" width="{pw}" height="22" rx="11" fill="#ffffff" '
              f'stroke="{stroke}" stroke-width="1.2"/>')
        f.text(D + wD - 10 - pw / 2, y + 25.5, pill, 11.5, fill=INK, anchor="middle", mono=pill.startswith(("-", "qc.")))
        spans = "".join(f'<tspan fill="{c}">{italic(escape(t))}</tspan>' for t, c in parts)
        f.add(f'<text x="{D + 14}" y="{y + 50}" font-family="{SANS}" font-size="12.5" fill="{SUB}">{spans}</text>')
    yr = ys[-1] + hs + 24
    f.doc(D + (wD - 220) / 2, yr, 220, h, "grey", "report.html", "all steps in one page")

    cy = lambda y: y + h / 2
    f.arrow([(A + wA, cy(r[0])), (B, cy(r[0]))], "grey")
    f.arrow([(A + wA, cy(r[1])), (B, cy(r[1]))], "orange")
    f.arrow([(A + wA, cy(r[2])), (B, cy(r[2]))], "green")
    f.arrow([(B + wB / 2, r[0] + h), (B + wB / 2, r[1])], "grey")
    f.arrow([(B + wB / 2, r[1] + h), (B + wB / 2, r[2])], "grey")
    f.arrow([(B + wB, cy(r[2])), (C, cy(r[2]))], "grey")
    f.arrow([(C + wC / 2, r[2] + h - 4), (C + wC / 2, r[3])], "grey")
    # final annotation into the QC steps, QC steps into the report
    f.arrow([(C + wC, cy(r[2])), (xl, cy(r[2]))], "grey", head=False)
    f.dot(xl, cy(r[2]), "grey")
    f.arrow([(xl, ys[0] + hs / 2), (xl, ys[-1] + hs / 2)], "grey", head=False)
    for y in ys:
        f.arrow([(xl, y + hs / 2), (D, y + hs / 2)], "grey")
        f.arrow([(D + wD, y + hs / 2), (xr, y + hs / 2)], "grey", head=False)
    f.arrow([(xr, ys[0] + hs / 2), (xr, cy(yr)), (D + (wD + 220) / 2, cy(yr))], "grey")

    # legend: what a step needs
    yl = r[3] + 20
    f.text(A, yl, "The color of a step says what it needs:", 13, "bold", INK)
    for i, (c, s) in enumerate((("grey", "the annotation, always runs"),
                                ("orange", "transcripts: reads, Iso-Seq or --stringtie"),
                                ("blue", "hints: proteins and/or transcripts"),
                                ("purple", "the option on the step and the data it names"))):
        fill, stroke = COL[c]
        yy = yl + 22 + 24 * i
        f.add(f'<rect x="{A}" y="{yy - 12}" width="22" height="16" rx="4" fill="{fill}" stroke="{stroke}" stroke-width="1.4"/>')
        f.text(A + 32, yy + 1, s)

    f.text(A, yr + 50, "<stem> is combined_evidence in evidence modes, else <tool>_ab_initio / <tool>_evidence;", fill=SUB)
    f.text(A, yr + 71, "<tool> is tiberius or vipsania. Evidence runs without a gene finder skip post-processing.")
    return f


out = sys.argv[1]
for name, fig in (("overview", overview()), ("combined_flow", combined_flow()),
                  ("insect_flow", insect()), ("postprocess_flow", postprocess())):
    with open(f"{out}/{name}.svg", "w") as fh:
        fh.write(fig.svg())
