# Post-processing and quality control

After the integration step (or after the gene finder in mode `abinitio`)
the final protein-coding annotation runs through the subworkflow
`POSTPROCESS` ([subworkflows/postprocess.nf](../subworkflows/postprocess.nf)).
It is the port of what BRAKER4 and GALBA2 do after gene prediction, without
the steps that need AUGUSTUS. Evidence runs without a gene finder have no
final annotation and skip it.

`<stem>` below is `<tool>_evidence`, or `<tool>_ab_initio` in mode
`abinitio`; `<tool>` is `tiberius` or `vipsania`.

| Step | Runs | Parameters |
| --- | --- | --- |
| [Sanity filter](#sanity-filter) | always (report only with `postprocess.sanity_filter: false`) | `postprocess.sanity_filter` |
| [UTRs](#utrs) | modes `rnaseq`, `isoseq`, `mixed` | `postprocess.utr`, `postprocess.max_utr_extension` |
| [Final files](#final-files-and-the-gff3-contract) | always | |
| [Completeness](#completeness) | with `qc.busco_lineage` | `qc.busco_lineage`, `qc.busco_download_path`, `qc.compleasm`, `qc.busco` |
| [Evidence support](#evidence-support-and-statistics) | evidence modes | `qc.gene_support` |
| [Gene set statistics](#evidence-support-and-statistics) | always | `qc.statistics` |
| [OMArk](#omark) | `qc.omark: true` | `qc.omamer_db`, `qc.ete_taxa_path` |
| [gffcompare](#gffcompare) | with `qc.reference_annotation` | `qc.reference_annotation` |
| [ncRNA](#ncrna) | `ncrna.run: true` | `ncrna.rfam_dir`, `ncrna.trnascan_high_confidence`, `ncrna.lncrna` |
| [GO terms](#go-terms-with-fantasia-lite) | `fantasia.run: true` (GPU) | `fantasia.*` |
| [Report](#report) | always | `qc.report` |

All parameters are in [parameters.md](parameters.md#post-processing-and-quality-control).
On the command line: `--busco_lineage`, `--reference_annotation`, `--ncrna`,
`--fantasia`.

## Sanity filter

`bin/sanity_filter_gff3.py` checks every transcript of the merged annotation
(`intermediate/<tool>_merged.gff3`, or the *ab initio* file in mode
`abinitio`) and writes `qc/sanity_filter.tsv` (`gene_id transcript_id action
reason`, counts in `#` lines on top):

| Reason | Action |
| --- | --- |
| `no_cds`: a transcript without CDS | removed |
| `cds_length_mod3`: CDS length, less the phase of the first CDS (a 5' partial transcript starts mid-codon), not a multiple of 3 | removed |
| `cds_structure`: CDS segments overlap or lie on different strands or sequences | removed |
| `cds_outside_exon`: a CDS segment outside the exons of its transcript | removed |
| `internal_stop`: a stop codon before the last codon | removed |
| `extended_stop`: the CDS ends without a stop codon and the next three bases are one; the CDS (and its exon) is extended by 3 bp | kept |
| `no_stop`: the CDS ends without a stop codon (3' partial; a stop codon is only added when all three bases lie in the exon of the CDS end, unless no exon follows; one that reaches into the intron is not, as the GT donor site fakes TA\|G) | kept |
| `non_atg_start`: the first codon is not ATG | kept |

A gene without transcripts is removed. Unlike BRAKER4, broken transcripts are
removed, not predicted again (there are no AUGUSTUS parameters).
`postprocess.sanity_filter: false` writes the report (action `flagged`) and
changes nothing.

**Stop codons.** The CDS of the final annotation includes the stop codon, as
in AUGUSTUS output. This was checked on the final *T. rubripes* annotations of
the TransDecoder flow (2026-10-07): 41,297 of 41,299 Tiberius transcripts and
41,125 of 41,150 Vipsania transcripts end with a stop codon inside the CDS.
The sanity filter on these files removed 0 Tiberius transcripts and 17
Vipsania transcripts (all `cds_length_mod3`), extended none and noted
`no_stop` for 2 and 8 transcripts. (`gffread -y` drops the terminal stop
codon from the protein sequences, so the proteins end without `*`.)

## UTRs

`bin/add_utrs_from_stringtie.py`, a GFF3 port of BRAKER4's
`stringtie2utr.py`: coding transcripts without UTRs get UTRs from the
StringTie assemblies of the run (short reads, and the one assembly of the
Drusilla flow or of `stringtie`, as `--stringtie`; Iso-Seq assemblies of the
TransDecoder flow as `--longread`, whose match wins over a short-read match).
A multi-exon transcript matches a StringTie transcript that contains all its
CDS introns; a single-exon transcript matches an overlapping StringTie
transcript without an intron inside the CDS; of several matches the longest
(spliced length) is used. The parts of the StringTie exons outside the CDS
become `five_prime_UTR`/`three_prime_UTR` (source `stringtie2utr`), and the
exons become the union of CDS and UTRs. The CDS never changes. UTRs reach at
most `postprocess.max_utr_extension` bp (5000) beyond the CDS and stop at the
nearest same-strand gene that is a barrier: one that has its own StringTie
match, or one that already has UTRs (the HC genes).
`qc/utr_report.tsv` lists per coding transcript `matched_by` (`short`, `long`,
`none`, `has_utr`) and the UTR lengths. The annotation before the UTRs is
`intermediate/<stem>_sanity_filtered.gff3`.

Limitations: a neighbour without a StringTie match and without UTRs does not
stop a UTR, so a read-through StringTie transcript can extend a UTR over it.
Differences to BRAKER4: the HC genes stop a neighbour's UTR (the original
only stops at a matched neighbour), all StringTie transcripts of an intron are
used (the original kept one per intron), StringTie transcript IDs are unique
per sequence and strand, the standard library replaces intervaltree, and a
UTR is added on a side of the CDS only if a StringTie exon contains that end
of the CDS (the original can join a UTR to the CDS by an intron StringTie
does not have).

## Final files and the GFF3 contract

`FINALIZE_ANNOTATION` publishes `<stem>.gff3`, `<stem>.gtf` (`gffread -T`;
transcript, exon and CDS lines, no gene lines, no UTR lines), the proteins of
all isoforms `<stem>_proteins.fa` (`gffread -y -S`) and the coding sequences
`<stem>_cds.fa` (`gffread -x`).

Every GFF3 that the pipeline publishes (`<stem>.gff3`,
`<stem>_with_ncRNA.gff3`, `<stem>_go.gff3`, `<stem>_with_ncRNA_go.gff3`,
`qc/<stem>_longest_isoform.gff3`) satisfies this contract; `bin/normalize_gff3.py`
writes it, `bin/validate_gff3.sh` checks it, and a file that fails the check
fails the run.

1. First line `##gff-version 3`; nine tab-separated columns; attributes
   `key=value` separated by `;`, values percent-encoded where GFF3 requires it.
2. Three levels: `gene` -> transcript -> `exon`, `CDS`, `five_prime_UTR`,
   `three_prime_UTR`. Transcripts are `mRNA` in the coding file; `tRNA`,
   `rRNA`, `snoRNA`, `lnc_RNA`, ... in the ncRNA file. No `start_codon`,
   `stop_codon` or `intron` features.
3. Unique IDs; every feature but a gene has exactly one `Parent`, written
   before it. Gene IDs `gene_000001`, transcripts `gene_000001.t1`, exons
   `<tx>.exon<N>`, CDS `<tx>.cds` (one ID for all segments of a transcript),
   UTRs `<tx>.utr5p<N>`/`<tx>.utr3p<N>`, numbered in genome order. ncRNA IDs
   are `<stem>-rRNA_<n>_<name>`, `<stem>-<seq>.trna<n>`, `<stem>-ncRNA_<i>`,
   `<stem>-lncRNA_<n>` with genes `gene-<rnaID>`.
4. Coding genes have `gene_biotype=protein_coding`; ncRNA genes their biotype.
5. The gene span is the union of its transcripts, a transcript span the union
   of its exons; CDS and UTR segments lie in exons; exons of a transcript do
   not overlap; the CDS includes the stop codon; CDS phases follow from the
   first segment; features are sorted by sequence, gene start, transcript and
   start.
6. `gt gff3validator` (GenomeTools 1.6.5) accepts the file, and `gffread -E`
   prints no warning.

The real *T. rubripes* annotations of both gene finders pass after the
sanity filter.

## Completeness

With `qc.busco_lineage` (e.g. `eukaryota_odb12`, `vertebrata_odb12`; a name
without `_odbNN` gets `_odb12`), compleasm and BUSCO assess the genome and the
proteome:

| File | Content |
| --- | --- |
| `qc/compleasm_genome/summary.txt`, `qc/compleasm_proteins/summary.txt` | compleasm 0.2.9 |
| `qc/busco_genome_short_summary.txt`, `qc/busco_proteins_short_summary.txt` | BUSCO 6.1.0 |
| `qc/completeness.tsv` | all four as one table (`assessment genome_or_proteome complete single duplicated fragmented missing n`); compleasm's F and I count as fragmented |
| `qc/<stem>_longest_isoform.gff3`, `qc/<stem>_longest_isoform_proteins.fa` | the proteome that is assessed |
| `qc/completeness.png` | bar chart, drawn by the report |

The proteome is the longest coding isoform of each gene
(`bin/longest_isoform.py`), not all isoforms as in BRAKER4, so that
alternative isoforms do not count as duplicated BUSCOs. Proteins of 100,000 aa
or more are left out (hmmsearch refuses them).

The lineage is read from `qc.busco_download_path` (default
`~/.cache/paludamentum/busco`), layout `<dir>/lineages/<lineage>/`. A missing
lineage is downloaded once on the submitting host (label `download`) with
`compleasm.py download`, which also fetches what compleasm checks before
every run (`file_versions.tsv`, the `eukaryota` lineage, the placement files
and their `.done` markers); the lineage archive is then extracted completely,
so that BUSCO (`--offline`) finds `dataset.cfg`. A cache made with
`busco --download` alone lacks the compleasm markers and is completed by the
download step. compleasm runs through `bin/compleasm_wrapper.py` (BRAKER4),
which lets compleasm use the cached files when the compute node has no
internet (compleasm's own `URLError` class otherwise shadows urllib's), and
leaves out over-long proteins. compleasm 0.2.7 (the version of the plan) fails
on the current `file_versions.tsv` of busco-data.ezlab.org; 0.2.9 works.

Without `qc.busco_lineage` there is no completeness assessment; the
pipeline and the launcher say so in one line.

## Evidence support and statistics

`qc/gene_support.tsv` (`bin/gene_support.py`, port of BRAKER4's
`gene_support_summary.py`): per coding transcript, the CDS introns and CDS
segments supported by RNA-Seq hints, protein hints and any hint of
`hintsfile.gff`. An intron is supported by an `intron` hint with the same
coordinates and strand, a CDS segment by an overlapping `CDSpart`, `exonpart`,
`CDS` or `exon` hint. Hint sources of Paludamentum:

| `src=` | Evidence class | Written by | Feature types |
| --- | --- | --- | --- |
| `E` | RNA-Seq (short reads and Iso-Seq) | `bam2hints --intronsonly` and `filterIntronsFindStrand.pl` on the merged short-read and Iso-Seq BAMs; pyVARUS `hints.gff`; summed by `merge_intron_hints.py` | `intron` |
| `P` | protein | `aln2hints.pl --prg=miniprot` on miniprothint's scored miniprot alignments | `CDSpart`, `intron` |
| `C` | both | not written by Paludamentum; counted as both, as in BRAKER4 | |

miniprothint's high-confidence hints (`hc.gff`) are not part of
`hintsfile.gff`; the Drusilla flow uses them. Short reads and Iso-Seq both
write `src=E` and cannot be told apart. Not in mode `abinitio` (no hints).

`qc/gene_set_statistics.txt` and the plots `qc/isoform_and_exon_structure.png`,
`qc/transcript_lengths.png`, `qc/introns_per_gene.png` and, with the support
table, `qc/evidence_support.png` (`bin/gene_set_statistics.py`, port of
BRAKER4's `gene_set_statistics.py`, PNG only).

## OMArk

`qc.omark: true` with `qc.omamer_db` (the OMAmer database, e.g. `LUCA.h5` of
<https://omabrowser.org/All/>, about 15 GB, downloaded by you): OMAmer places
the proteins of all isoforms, OMArk reads them grouped per gene and writes
`qc/omark_summary.txt`.

OMArk reads the NCBI taxonomy through ete3, which needs it as a database
(`taxa.sqlite`). The pipeline builds it once on the submitting host (label
`download`, image `omark`) in `qc.ete_taxa_path`, by default
`~/.cache/paludamentum/ncbi_taxonomy`, from `taxdump.tar.gz` of
<https://ftp.ncbi.nih.gov/pub/taxonomy/>, which is downloaded there (md5
checked) unless the directory already holds it: for a submitting host without
internet, put the tarball there. The OMArk task gets the database staged
(`omark -e`) and needs no internet. The build takes about six minutes and
3 GB of memory (`taxa.sqlite` is about 0.8 GB). The database stays as built;
delete `taxa.sqlite` to rebuild it from a newer tarball. (ete3 itself never reads a
tarball from `~/.etetoolkit`: without a database it fetches the md5 and the
tarball from NCBI, which is why the task cannot build it on a node without
internet.)

## gffcompare

With `qc.reference_annotation` (GFF3 or GTF, gzipped or not), the CDS of the final annotation
are compared with the CDS of the reference, as the benchmarks of this pipeline
and BRAKER4 do: `gffcompare --strict-match -e 3 -T` on CDS lines.
`qc/gffcompare.stats` has sensitivity and precision at base, exon, intron,
transcript and gene (locus) level.

## ncRNA

`ncrna.run: true` (commands of BRAKER4 `rules/ncrna/`):

| File | Tool |
| --- | --- |
| `ncrna/rRNA.gff3` | pybarrnap 0.5.1 (`--kingdom euk`) |
| `ncrna/tRNAs.gff3`, `ncrna/tRNAs.txt` | tRNAscan-SE 2.0.12 (`-E`; with `ncrna.trnascan_high_confidence` the EukHighConfidenceFilter) |
| `ncrna/infernal.tblout`, `ncrna/ncRNAs_infernal.gff3` | Infernal 1.1.5 `cmscan --cut_ga --rfam --nohmmonly` against Rfam 15.1, one task per genome chunk; `bin/infernal_to_gff3.py` keeps the hits above the gathering threshold and, of a set of overlapping hits, the best-scoring one (`olp` other than `=`), and types each by its Rfam family |
| `ncrna/lncRNAs.gff3`, `ncrna/feelnc_classifier.txt` | FEELnc 0.2 on the merged StringTie assemblies (modes with transcripts, `ncrna.lncrna`; see below) |
| `<stem>_with_ncRNA.gff3` | the final annotation plus the ncRNA genes (`bin/merge_ncrna_gff3.py`; priority rRNA > tRNA > Rfam > lncRNA; an ncRNA that overlaps coding exons is dropped) |

Rfam 15.1 (`Rfam.cm`, `Rfam.clanin`, the `cmpress` index) is read from
`ncrna.rfam_dir`, else downloaded once into `~/.cache/paludamentum/rfam` on
the submitting host and checked against its SHA-256. A `rfam_dir` without the
`cmpress` index is pressed in each task. cmscan is the slowest step on large
genomes; it runs on the genome chunks of the gene finder split (20 Mb, at
most 20 chunks).

FEELnc (`bin/feelnc_to_gff3.py` for the GFF3): `FEELnc_filter.pl` keeps the
assembled transcripts of at least 200 bp with more than one exon that do not
overlap the final annotation (the candidates); `FEELnc_codpot.pl --mode=shuffle`
trains a random forest on the annotated mRNAs and on shuffled copies of them
and keeps the candidates without coding potential; `FEELnc_classifier.pl`
lists the coding genes next to each lncRNA (`feelnc_classifier.txt`, by
StringTie transcript ID, which the `lnc_RNA` features carry as `Alias`).
FEELnc cannot train on fewer than 100 candidates or fewer than 100 annotated
transcripts: `lncRNAs.gff3` is then empty but for a comment line that says so
(also in the task log), as when no candidate is without coding potential. Any
other FEELnc error fails the task and the run; `ncrna.lncrna: false` skips
the step.

## GO terms with FANTASIA-Lite

`fantasia.run: true` (GPU, at least 15 GB free GPU memory; commands of
BRAKER4 `rules/postprocessing/run_fantasia.smk`): ProtT5 embeddings of all
proteins and GO terms by similarity to the FANTASIA lookup table.

| Parameter | Content |
| --- | --- |
| `fantasia.hf_cache_dir` | Hugging Face cache with `Rostlab/prot_t5_xl_uniref50` (the snapshot with `pytorch_model.bin` is used, as in BRAKER4) |
| `fantasia.lookup_dir` | `lookup_table.npz`, `annotations.json`, `accessions.json` (Zenodo record 17720428) |
| `fantasia.min_score` | minimum score of a GO term (0.5) |
| `fantasia.additional_params` | appended to `fantasia_pipeline.py` |

Outputs: `<stem>_go.gff3` (and `<stem>_with_ncRNA_go.gff3` with `ncrna.run`)
with `Ontology_term` on mRNAs and genes, `qc/fantasia/results.csv`,
`qc/fantasia/failed_sequences.csv`, `qc/fantasia/fantasia_summary.txt`,
`qc/fantasia/fantasia_go_terms.tsv`, `qc/fantasia/fantasia_go_categories.png`
(a placeholder figure when no GO term reaches `fantasia.min_score`).
The FANTASIA-Lite image is run by Nextflow; `conf/base.config` binds
`bin/fantasia_generate_embeddings.py` (BRAKER4's patched embedding script)
over the image's copy. FANTASIA-Lite was validated by BRAKER4 on A100 GPUs
only. Give the process a suitable queue in your config (see
[conf/user_hpc_template.config](../conf/user_hpc_template.config)).

`FANTASIA_GPU_CHECK` (`bin/fantasia_gpu_check.sh`, same labels and so the
same image and queue as `FANTASIA_ANNOTATE`) probes for a CUDA GPU with at
least 15 GB free memory at the start of the run: a run whose `gpu` label
gives no usable GPU stops in the first minutes, instead of failing in
`FANTASIA_ANNOTATE` after hours, with the report never written. The same
probe runs again in `FANTASIA_ANNOTATE` before the model loads.

## Report

`report.html` (`bin/paludamentum_report.py`): run summary, output files, gene
set statistics, completeness, evidence support, sanity filter, UTRs, ncRNA,
OMArk, gffcompare, FANTASIA, software versions (`qc/software_versions.tsv`)
and the references of `citations.md`. One self-contained file; sections
without data are left out.

## Outputs

```
<outdir>/
  <stem>.gff3, <stem>.gtf, <stem>_proteins.fa, <stem>_cds.fa
  <stem>_with_ncRNA.gff3                       ncrna.run
  <stem>_go.gff3, <stem>_with_ncRNA_go.gff3    fantasia.run
  report.html
  intermediate/
    <tool>_merged.gff3                         evidence modes: ab initio + HC genes, before the sanity filter
    <tool>_ab_initio.gff3                      ab initio predictions, before the sanity filter
    <stem>_sanity_filtered.gff3                after the sanity filter, before the UTRs
  qc/
    sanity_filter.tsv, utr_report.tsv
    gene_support.tsv, gene_set_statistics.txt, *.png
    software_versions.tsv
    compleasm_genome/summary.txt, compleasm_proteins/summary.txt,
    busco_genome_short_summary.txt, busco_proteins_short_summary.txt,
    completeness.tsv, completeness.png,
    <stem>_longest_isoform.gff3, <stem>_longest_isoform_proteins.fa     qc.busco_lineage
    omark_summary.txt                          qc.omark
    gffcompare.stats                           qc.reference_annotation
    fantasia/                                  fantasia.run
  ncrna/                                       ncrna.run
```

## Checks on real data

On the final *T. rubripes* annotations of the TransDecoder flow (Tiberius
vertebrates and Vipsania etb1go6q, runs of 2026-09-26) with the reference
annotation of the benchmark (`annot_cds.gff`) (gffcompare `--strict-match -e 3 -T` on CDS lines,
sensitivity/precision):

| Annotation | Transcript level | Gene (locus) level |
| --- | --- | --- |
| Tiberius, before the sanity filter | 44.3 / 43.3 | 72.2 / 67.3 |
| Tiberius, after the sanity filter | 44.3 / 43.3 | 72.2 / 67.3 |
| Tiberius, after the UTRs | 44.3 / 43.3 | 72.2 / 67.3 |
| Vipsania, before the sanity filter | 39.6 / 38.8 | 64.6 / 58.9 |
| Vipsania, after the sanity filter | 39.6 / 38.8 | 64.6 / 58.9 |

The UTR step changed no CDS (compared line by line) and added UTRs to 14,421
of the 23,946 Tiberius transcripts without UTRs (the 17,353 HC transcripts
had UTRs already). compleasm 0.2.9 on the longest isoforms (eukaryota_odb12):
S 79.84 %, D 15.50 %, F 2.33 %, M 2.33 %.

## Not ported

From BRAKER4 and GALBA2: re-prediction of transcripts with in-frame stop
codons (needs AUGUSTUS parameters; the sanity filter removes them instead),
TSEBRA (`bin/merge_annotations.py` is the merge), AGAT (the normaliser
replaces it), the BUSCO rescue `best_by_compleasm` (deferred: only the
Drusilla flow drops *ab initio* genes; to be decided with the completeness
numbers of the Drusilla and TransDecoder flows), repeat masking (a step before
the prediction), translation tables other than 1.
