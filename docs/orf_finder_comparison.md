# ORF finder of the HC gene step: TransDecoder, TD2 and Drusilla

The high-confidence (HC) gene step gets its ORFs from TransDecoder 5.7.1. This
page compares it with two alternatives for a Vipsania run on *Takifugu
rubripes*:

- **TD2** ([TransDecoder2](https://github.com/Markusjsommer/TD2), version
  1.1.0), which scores ORFs with PSAURON. Set `transdecoder: td2` in the
  params file.
- **Drusilla**, the ORF flow of branch `drusilla`, which is selected there for
  vertebrate models.

TD2 is slightly better than TransDecoder (+0.8 gene F1). The Drusilla flow is
about 11 points ahead of both.

## Result

Final gene set, accuracy in %:

| ORF flow | Gene Sn | Gene Pr | Gene F1 | Transcript Sn | Transcript Pr | Exon Sn | Exon Pr | Transcripts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TransDecoder 5.7.1 | 69.3 | 63.3 | 66.16 | 42.7 | 41.4 | 88.1 | 88.3 | 41 147 |
| TD2 1.1.0 | 70.3 | 64.0 | 67.00 | 43.4 | 41.8 | 88.3 | 88.4 | 41 489 |
| Drusilla | 74.9 | 81.2 | 77.92 | 48.6 | 57.7 | 86.4 | 92.3 | 33 613 |

The Vipsania *ab initio* prediction is the same in all three runs: gene Sn 46.7,
Pr 43.7, F1 45.15 (23 788 transcripts).

ORF set alone, before the merge with the Vipsania genes:

| ORF set | Gene Sn | Gene Pr | Gene F1 | Transcript Sn | Transcript Pr | Exon Sn | Exon Pr | Transcripts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HC genes, TransDecoder 5.7.1 | 53.6 | 89.1 | 66.93 | 32.9 | 75.7 | 60.7 | 98.0 | 17 359 |
| HC genes, TD2 1.1.0 | 55.3 | 90.2 | 68.56 | 34.0 | 76.8 | 62.1 | 98.1 | 17 701 |
| Drusilla ORFs, stop/start fixed | 67.1 | 82.7 | 74.09 | 43.0 | 67.5 | 78.3 | 95.1 | 25 434 |

## What the numbers mean

- **TD2 compared to TransDecoder.** The HC set gains 1.7 points in sensitivity
  and 1.1 in precision. After the merge with the Vipsania genes 0.8 points of
  gene F1 remain.
- **Drusilla compared to both.** The Drusilla ORFs alone are more sensitive
  than the HC genes (67.1 against 53.6 and 55.3) and less precise (82.7
  against 89.1 and 90.2).
- **Precision of the final set.** The TransDecoder flow merges the HC genes
  with all Vipsania genes and ends at 63 to 64 % gene precision. The Drusilla
  flow keeps only the Vipsania genes that pass its LightGBM filter and ends at
  81 %. A different ORF finder in the TransDecoder flow does not change that.

## TD2 and TransDecoder on the same transcripts

| | TransDecoder 5.7.1 | TD2 1.1.0 |
| --- | --- | --- |
| ORFs | 28 638 | 24 833 |
| complete | 22 939 | 22 651 |
| 5' partial | 3 734 | 374 |
| 3' partial | 1 515 | 1 157 |
| internal | 450 | 651 |
| on the minus strand of the transcript | 2 027 | 613 |
| HC transcripts after the DIAMOND filter | 17 359 | 17 701 |
| Wall time of `TD_ALL` | 11 min 24 s | 6 min 14 s |
| CPU use | 1 core | 9 cores |
| Peak memory | 0.6 GB | 3.2 GB |

Both ran with their defaults and on both strands of the transcripts: TD2 with
a minimal ORF length of 90 aa and a PSAURON cutoff of 0.5, TransDecoder with a
minimal ORF length of 100 aa.

The FASTA headers of TD2 follow the TransDecoder format with two exceptions.
Minus-strand ORFs have their coordinates as `<start>-<end>` with start > end,
which `TD_ALL` rewrites to the TransDecoder order. `len:` counts the stop
codon; the HC steps do not use it.

## Setup

| Item | Value |
| --- | --- |
| Genome, reference | *T. rubripes* (GCF_901000725.2) of the Tiberius vertebrate test set, reference CDS annotation of the same assembly |
| Evidence | one VARUS RNA-Seq BAM (mode `rnaseq`), proteins of four related species |
| Gene finder | Vipsania, model `etb1go6q` (Vertebrata) |
| TransDecoder and TD2 runs | Paludamentum aab81d8, image `tiberius:2.0.8`, 2026-09-30. Both reuse one Vipsania prediction (`vipsania.result`) |
| Drusilla run | branch `drusilla` with hint rescue, same inputs, 2026-09-27 |
| Evaluation | `gffcompare --strict-match -e 3 -T` (0.12.6) on the CDS lines; gene level is the locus level of gffcompare |

TD2 is not in the `tiberius:2.0.8` image. For this comparison TD2 1.1.0,
psauron 1.1.3 and the CPU build of PyTorch were installed with `pip` into a
directory that the site config adds to `PYTHONPATH` and `PATH` of `TD_ALL`.


## Limits

- One species and one run per flow.
- TransDecoder is not deterministic. An earlier run with the same inputs gave
  gene F1 66.11 for the final set and 66.83 for the HC genes.
- TD2 options such as `--precise` (`td2_predict_args`) and the strand-specific
  mode were not tested.
