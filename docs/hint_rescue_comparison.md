# Hint rescue: Paludamentum compared with the original scripts

The hint rescue predicts partial gene finder genes again with Tiberius, using
the protein hints of the best alignment chain at each locus. Paludamentum
(processes `HINT_RESCUE_LOCI` and `HINT_RESCUE_TIBERIUS` in
`modules/drusilla.nf`) uses the scripts of the original integration
(`prepare_hint_rescue_loci.py`, `filter_and_merge_rescue_gtf.py`,
`chainedHints.py`) unchanged. The inputs and the Tiberius call differ in
several places from the original SLURM scripts.

## Differences

| | Original scripts | Paludamentum |
|---|---|---|
| `seq_len` of Tiberius | default of the model config (400 050): each locus in one piece | 99 990 (`drusilla.rescue_seq_len`): loci longer than 100 kb are split into windows |
| Tiberius model | always `vertebrates` | the model of the run (*Bos taurus*: `mammalia_softmasking_v2`); `vertebrates` for Vipsania runs (`drusilla.rescue_model_cfg`) |
| ORFs passed as `--orfs_gtf` | raw Drusilla ORFs, neither fixed nor collapsed | final Drusilla ORFs, stop/start fixed and subsequence-collapsed (`drusilla_orfs.gtf`). The agreement check that would read them is switched off in the script, so the choice has no effect |
| Partial and correct genes | `tiberius_lgb_partial.gtf` and `tiberius_lgb_correct.gtf` from the raw Tiberius GTF | from the merged *ab initio* GFF3, converted to transcript and CDS lines. Partial: most likely class `partial` among the transcripts with P(partial) + P(correct) ≥ 0.5; how the original partial file was built is not recorded |
| Final gene set | `cat` of LightGBM "correct" + rescue + ORFs | `merge_annotations.py --mode full` of the same three sets |

## What the differences mean

- **`seq_len`.** With 2 305 loci (*T. rubripes*), each padded to 400 kb,
  Tiberius needed more than 90 GB of memory, so Paludamentum uses 99 990.
  Windows can cut genes at their borders in the few loci longer than 100 kb
  (median 66 kb, maximum 399 kb).
- **Model.** On *Bos taurus* the rescue uses the mammalia model of the run,
  not `vertebrates`.
- **ORFs.** With the agreement check switched off, the ORF file does not
  change the loci. On *T. rubripes* the rescue was byte-identical for two
  different RNA-Seq BAMs.
- **Final gene set.** `cat` and the merge give the same gene-level accuracy;
  the merge removes duplicate transcripts, which raises transcript precision.

## Effect of the rescue in Paludamentum (gene F1)

| | without rescue | with rescue | effect |
|---|---|---|---|
| *T. rubripes*, Tiberius vertebrates (poster BAM) | 80.44 | 80.26 | −0.18 |
| *T. rubripes*, Vipsania Vertebrata (poster BAM) | 77.75 | 77.92 | +0.17 |
| *Bos taurus*, Tiberius mammalia (poster BAM) | 79.40 | 79.46 | +0.06 |
| *Bos taurus*, 6 new pyVARUS BAMs (with and without Logan) | 78.09–79.30 | 78.32–79.62 | +0.14 to +0.33 |
