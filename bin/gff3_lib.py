"""
Shared GFF3 helpers of the post-processing scripts in bin/ (sanity filter,
normaliser, longest isoform, UTRs, gene support). Standard library only.

Nextflow mounts bin/ into every task, so the scripts import this module
from their own directory.

Model: an Annotation holds Genes in file order; a Gene holds its
transcript-level features (mRNA, tRNA, ...); a Transcript holds its children
(exon, CDS, UTRs). Attribute values are decoded (percent-encoding removed)
and kept as lists, because GFF3 separates multiple values by commas.

The GFF3 contract of the published files is in docs/postprocessing.md.
"""
from __future__ import annotations

import gzip
import sys
from collections import defaultdict
from typing import Dict, Iterable, Iterator, List, Optional, TextIO, Tuple
from urllib.parse import quote, unquote

TRANSCRIPT_TYPES = {
    "mRNA", "transcript", "tRNA", "rRNA", "snoRNA", "snRNA", "scaRNA", "miRNA",
    "lnc_RNA", "antisense_RNA", "ribozyme", "ncRNA", "primary_transcript",
    "pseudogenic_transcript", "SRP_RNA", "RNase_P_RNA", "RNase_MRP_RNA",
    "telomerase_RNA", "vault_RNA", "Y_RNA", "guide_RNA", "piRNA", "siRNA",
}
UTR5_TYPES = {"five_prime_UTR", "5UTR", "5'UTR"}
UTR3_TYPES = {"three_prime_UTR", "3UTR", "3'UTR"}
UTR_TYPES = UTR5_TYPES | UTR3_TYPES | {"UTR"}
# Features that the contract does not allow; they are dropped on reading
DROPPED_TYPES = {"start_codon", "stop_codon", "intron", "transcription_start_site",
                 "transcription_end_site", "tss", "tts"}
# Order of the child features of a transcript at the same start
CHILD_ORDER = {"exon": 0, "five_prime_UTR": 1, "CDS": 2, "three_prime_UTR": 3}

# Characters that GFF3 requires to be escaped in attribute values (the comma
# is the separator of multiple values and is escaped inside one value)
_SAFE = "".join(chr(c) for c in range(0x20, 0x7f) if chr(c) not in ";=&,%\t")

STOP_CODONS = {"TAA", "TAG", "TGA"}
_COMP = str.maketrans("ACGTRYKMBVDHNacgtrykmbvdhn", "TGCAYRMKVBHDNtgcayrmkvbhdn")
_BASES = "TCAG"
_AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {a + b + c: _AMINO[16 * i + 4 * j + k]
               for i, a in enumerate(_BASES) for j, b in enumerate(_BASES) for k, c in enumerate(_BASES)}


class GFF3Error(Exception):
    """A GFF3 that violates the structure the scripts need; the message names the feature."""


# ---------------------------------------------------------------- features

class Feature:
    __slots__ = ("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attrs", "line")

    def __init__(self, seqid: str, source: str, ftype: str, start: int, end: int, score: str = ".",
                 strand: str = ".", phase: str = ".", attrs: Optional[Dict[str, List[str]]] = None,
                 line: int = 0):
        self.seqid = seqid
        self.source = source
        self.type = ftype
        self.start = int(start)
        self.end = int(end)
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attrs: Dict[str, List[str]] = attrs if attrs is not None else {}
        self.line = line

    @property
    def id(self) -> Optional[str]:
        values = self.attrs.get("ID")
        return values[0] if values else None

    @property
    def parents(self) -> List[str]:
        return self.attrs.get("Parent", [])

    def get(self, key: str, default: Optional[str] = None) -> Optional[str]:
        values = self.attrs.get(key)
        return ",".join(values) if values else default

    def set(self, key: str, value: str) -> None:
        self.attrs[key] = [value]

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def copy(self, **changes) -> "Feature":
        new = Feature(self.seqid, self.source, self.type, self.start, self.end, self.score, self.strand,
                      self.phase, {k: list(v) for k, v in self.attrs.items()}, self.line)
        for key, value in changes.items():
            setattr(new, key, value)
        return new

    def describe(self) -> str:
        where = f"line {self.line}: " if self.line else ""
        return f"{where}{self.type} {self.id or '(no ID)'} {self.seqid}:{self.start}-{self.end}{self.strand}"

    def to_line(self) -> str:
        return "\t".join([self.seqid, self.source, self.type, str(self.start), str(self.end), self.score,
                          self.strand, self.phase, format_attributes(self.attrs)])


def parse_attributes(text: str) -> Dict[str, List[str]]:
    """GFF3 column 9 as {key: [decoded values]}; keys keep their order."""
    attrs: Dict[str, List[str]] = {}
    text = text.strip()
    if not text or text == ".":
        return attrs
    for part in text.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" not in part:
            # a GTF-style or broken attribute; keep it under its own text
            attrs.setdefault(unquote(part), [])
            continue
        key, value = part.split("=", 1)
        attrs[unquote(key.strip())] = [unquote(v) for v in value.split(",")]
    return attrs


def encode_value(value: str) -> str:
    return quote(value, safe=_SAFE)


def format_attributes(attrs: Dict[str, List[str]]) -> str:
    """ID first, Parent second, then the other keys in their order."""
    if not attrs:
        return "."
    keys = [k for k in ("ID", "Parent") if k in attrs] + [k for k in attrs if k not in ("ID", "Parent")]
    parts = []
    for key in keys:
        values = attrs[key]
        if not values:
            continue
        parts.append(f"{encode_value(key)}={','.join(encode_value(v) for v in values)}")
    return ";".join(parts) if parts else "."


def read_features(handle: TextIO) -> Iterator[Feature]:
    """Feature lines of a GFF3; comments, directives and a FASTA section are skipped."""
    for number, raw in enumerate(handle, 1):
        line = raw.rstrip("\n\r")
        if line.startswith("##FASTA"):
            return
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) != 9:
            raise GFF3Error(f"line {number}: {len(cols)} tab-separated columns instead of 9: {line[:200]}")
        seqid, source, ftype, start, end, score, strand, phase, attributes = cols
        try:
            begin, stop = int(start), int(end)
        except ValueError:
            raise GFF3Error(f"line {number}: start or end is not a number: {line[:200]}")
        if begin > stop:
            raise GFF3Error(f"line {number}: start {begin} > end {stop}")
        yield Feature(seqid, source, ftype, begin, stop, score.strip() or ".", strand.strip() or ".",
                      phase.strip() or ".", parse_attributes(attributes), number)


def open_text(path: str) -> TextIO:
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


# ---------------------------------------------------------------- hierarchy

class Transcript:
    def __init__(self, feature: Feature):
        self.feature = feature
        self.children: List[Feature] = []

    @property
    def id(self) -> str:
        return self.feature.id

    def of_type(self, *types: str) -> List[Feature]:
        return sorted((c for c in self.children if c.type in types), key=lambda f: (f.start, f.end))

    @property
    def exons(self) -> List[Feature]:
        return self.of_type("exon")

    @property
    def cds(self) -> List[Feature]:
        return self.of_type("CDS")

    @property
    def utrs(self) -> List[Feature]:
        return sorted((c for c in self.children if c.type in UTR_TYPES), key=lambda f: (f.start, f.end))

    @property
    def is_coding(self) -> bool:
        return any(c.type == "CDS" for c in self.children)

    @property
    def cds_length(self) -> int:
        return sum(c.length for c in self.children if c.type == "CDS")

    def cds_in_transcription_order(self) -> List[Feature]:
        cds = self.cds
        return cds if self.feature.strand != "-" else cds[::-1]


class Gene:
    def __init__(self, feature: Feature):
        self.feature = feature
        self.transcripts: List[Transcript] = []

    @property
    def id(self) -> str:
        return self.feature.id


class Annotation:
    def __init__(self):
        self.genes: List[Gene] = []

    def transcripts(self) -> Iterator[Tuple[Gene, Transcript]]:
        for gene in self.genes:
            for tx in gene.transcripts:
                yield gene, tx


def read_gff3(path: str, strict: bool = True) -> Annotation:
    """
    Read a three-level GFF3 (gene -> transcript -> exon/CDS/UTR).

    strict: a missing Parent, a duplicated ID or a feature that does not fit
    the three levels is a GFF3Error naming the feature. Not strict: such
    features are skipped with a warning on stderr. start_codon, stop_codon
    and intron features are dropped in both cases. A transcript without a
    gene gets one (ID <transcript>.gene); a CDS, exon or UTR may share its ID
    with the other parts of the same discontinuous feature.
    """
    with open_text(path) as handle:
        features = list(read_features(handle))
    return build_annotation(features, strict)


def _problem(strict: bool, message: str) -> None:
    if strict:
        raise GFF3Error(message)
    print(f"[WARN] {message}; skipped", file=sys.stderr)


def build_annotation(features: Iterable[Feature], strict: bool = True) -> Annotation:
    features = [f for f in features if f.type not in DROPPED_TYPES]
    ann = Annotation()
    genes: Dict[str, Gene] = {}
    ids: Dict[str, Feature] = {}

    def unique(feat: Feature) -> bool:
        fid = feat.id
        if fid is None:
            return True
        other = ids.get(fid)
        if other is None:
            ids[fid] = feat
            return True
        # the parts of one discontinuous feature (CDS, exon, UTR) share the ID
        if (other.type == feat.type and other.parents == feat.parents and feat.type != "gene"
                and feat.type not in TRANSCRIPT_TYPES):
            return True
        _problem(strict, f"duplicate ID {fid}: {feat.describe()} and {other.describe()}")
        return False

    for feat in features:
        if feat.type != "gene" or not unique(feat):
            continue
        if feat.id is None:
            _problem(strict, f"gene without ID: {feat.describe()}")
            continue
        if feat.parents:
            _problem(strict, f"{feat.describe()} has a Parent")
            continue
        gene = Gene(feat)
        genes[feat.id] = gene
        ann.genes.append(gene)

    child_types = set(CHILD_ORDER) | UTR_TYPES
    transcripts: Dict[str, Transcript] = {}
    children: List[Feature] = []
    for feat in features:
        if feat.type == "gene":
            continue
        parents = feat.parents
        is_tx = feat.type in TRANSCRIPT_TYPES or (
            feat.type not in child_types and parents and all(p in genes for p in parents))
        if not is_tx:
            if not parents and feat.type not in child_types:
                print(f"[WARN] {feat.describe()} is neither gene, transcript nor part of one; skipped",
                      file=sys.stderr)
            elif unique(feat):
                children.append(feat)
            continue
        if not unique(feat):
            continue
        if feat.id is None:
            _problem(strict, f"transcript without ID: {feat.describe()}")
            continue
        if len(parents) > 1:
            _problem(strict, f"{feat.describe()} has {len(parents)} parents")
            continue
        if parents:
            gene = genes.get(parents[0])
            if gene is None:
                _problem(strict, f"Parent {parents[0]} of {feat.describe()} is not a gene in the file")
                continue
        else:
            # a transcript without a gene gets one
            gid = f"{feat.id}.gene"
            gene = Gene(Feature(feat.seqid, feat.source, "gene", feat.start, feat.end, ".", feat.strand, ".",
                                {"ID": [gid]}, feat.line))
            genes[gid] = gene
            ann.genes.append(gene)
            feat.set("Parent", gid)
        tx = Transcript(feat)
        transcripts[feat.id] = tx
        gene.transcripts.append(tx)

    for feat in children:
        parents = feat.parents
        if not parents:
            _problem(strict, f"{feat.describe()} has no Parent")
            continue
        for parent in parents:
            tx = transcripts.get(parent)
            if tx is None:
                _problem(strict, f"Parent {parent} of {feat.describe()} is not a transcript in the file")
                continue
            tx.children.append(feat if len(parents) == 1 else feat.copy(attrs={**feat.attrs, "Parent": [parent]}))
    return ann


# ---------------------------------------------------------------- coordinates

def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge overlapping or adjacent 1-based inclusive intervals."""
    merged: List[List[int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]


def contained(feature: Feature, intervals: Iterable[Tuple[int, int]]) -> bool:
    return any(s <= feature.start and feature.end <= e for s, e in intervals)


def next_phase(phase: int, length: int) -> int:
    """Phase of the next CDS segment after one of this phase and length."""
    return (3 - ((length - phase) % 3)) % 3


def phase_of(feature: Feature) -> int:
    return int(feature.phase) if feature.phase in ("0", "1", "2") else 0


def set_phases(tx: Transcript) -> None:
    """Phases of all CDS segments from the first one and the segment lengths."""
    cds = tx.cds_in_transcription_order()
    if not cds:
        return
    phase = phase_of(cds[0])
    for seg in cds:
        seg.phase = str(phase)
        phase = next_phase(phase, seg.length)


def update_spans(gene: Gene) -> None:
    """Transcript span = union of its children, gene span = union of its transcripts."""
    for tx in gene.transcripts:
        if tx.children:
            tx.feature.start = min(c.start for c in tx.children)
            tx.feature.end = max(c.end for c in tx.children)
    if gene.transcripts:
        gene.feature.start = min(t.feature.start for t in gene.transcripts)
        gene.feature.end = max(t.feature.end for t in gene.transcripts)


# ---------------------------------------------------------------- sequences

def read_fasta(path: str, wanted: Optional[set] = None) -> Dict[str, str]:
    """FASTA as {first word of the header: upper-case sequence}; only `wanted` ids if given."""
    seqs: Dict[str, str] = {}
    name = None
    chunks: List[str] = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None and (wanted is None or name in wanted):
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0] if line[1:].split() else ""
                chunks = []
            elif name is not None and (wanted is None or name in wanted):
                chunks.append(line.strip())
        if name is not None and (wanted is None or name in wanted):
            seqs[name] = "".join(chunks).upper()
    return seqs


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def subseq(genome: Dict[str, str], seqid: str, start: int, end: int, strand: str) -> str:
    """Sequence of 1-based inclusive start..end, reverse complement on '-'."""
    seq = genome[seqid][start - 1:end]
    return revcomp(seq) if strand == "-" else seq


def spliced_cds(tx: Transcript, genome: Dict[str, str]) -> str:
    """Coding sequence in transcription order, starting at the phase of the first segment."""
    cds = tx.cds_in_transcription_order()
    seq = "".join(subseq(genome, c.seqid, c.start, c.end, c.strand) for c in cds)
    return seq[phase_of(cds[0]):] if cds else ""


def translate(seq: str) -> str:
    return "".join(CODON_TABLE.get(seq[i:i + 3], "X") for i in range(0, len(seq) - len(seq) % 3, 3))


# ---------------------------------------------------------------- writing

def sort_annotation(ann: Annotation) -> None:
    """Genes by seqid, start, end and ID; transcripts by start; children by start."""
    for gene in ann.genes:
        gene.transcripts.sort(key=lambda t: (t.feature.start, t.feature.end, t.id))
        for tx in gene.transcripts:
            tx.children.sort(key=lambda c: (c.start, CHILD_ORDER.get(c.type, 9), c.end))
    ann.genes.sort(key=lambda g: (g.feature.seqid, g.feature.start, g.feature.end, g.id))


def write_gff3(ann: Annotation, handle: TextIO) -> None:
    handle.write("##gff-version 3\n")
    for gene in ann.genes:
        handle.write(gene.feature.to_line() + "\n")
        for tx in gene.transcripts:
            handle.write(tx.feature.to_line() + "\n")
            for child in tx.children:
                handle.write(child.to_line() + "\n")


def gene_transcript_ids(ann: Annotation) -> Dict[str, List[str]]:
    groups: Dict[str, List[str]] = defaultdict(list)
    for gene, tx in ann.transcripts():
        groups[gene.id].append(tx.id)
    return groups
