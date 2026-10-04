# Drusilla flow

For vertebrate gene finder models, Paludamentum replaces the TransDecoder
high-confidence genes with Drusilla ORFs and filters the *ab initio*
predictions. The flow is derived from the one behind the Tiberius evidence
results on the GCB 2026 poster
([doi:10.13140/RG.2.2.24444.91521](https://doi.org/10.13140/RG.2.2.24444.91521)),
which used an earlier Drusilla model and had no stop/start codon fix, no hint
rescue and no merge of overlapping genes. With this flow, the gene-level F1
is 80.46 on *Takifugu rubripes* with the Tiberius `vertebrates` model and
79.71 on *Bos taurus* with `mammalia_softmasking_v2` (Paludamentum 0.4.0).
It runs when all of these hold:

- the run has transcripts (mode `rnaseq`, `isoseq` or `mixed`),
- the gene finder model is a Tiberius model whose `target_species` is
  `Vertebrata` or `Mammalia` (`vertebrates`, `mammalia*`), or Vipsania
  `Vertebrata` (`etb1go6q`); `drusilla.run: true` forces the flow for other
  vertebrate models,
- `drusilla.lgb_model` is not `null` (the default is the released model).

<p align="center">
  <img src="../figures/drusilla_flow.svg" alt="Drusilla flow: StringTie assembly, transcript filter, Drusilla ORFs and codon fix give the high-confidence genes; LightGBM filter and hint rescue treat the ab initio genes; both are merged into the final annotation" width="100%">
</p>

Steps (numbered as in the figure):

1. One StringTie assembly of all reads (short reads, Iso-Seq with `-L`, or
   both with `--mix`). An assembly made outside the pipeline can be given
   with `stringtie`; it must be that one assembly, with the `cov` and `TPM`
   attributes of StringTie (not the output of `stringtie --merge`).
2. Transcripts are kept if length >= 300, coverage >= 3 and TPM >= 1 (TPM >= 0.5
   for transcripts of 3000 nt or longer).
3. `drusilla annotate` predicts the ORFs of the kept transcripts, including
   ORFs truncated at the 3' or 5' end of a transcript. The released model
   (`vertebrates`: weights of the training run `cnn_lstm_vertebrates_run009`)
   is downloaded once on the submitting host and checked against the
   `weights_sha256` of its manifest in `drusilla/model_cfg/`, so the GPU nodes
   need no internet for it. With `cache_dir`, Drusilla takes the model from
   that cache instead (`drusilla models download vertebrates` with
   `DRUSILLA_CACHE_DIR` set to it fills it once).
4. Stop codon fix: ORFs with an early stop and truncated ORFs are extended to
   a stop codon supported by a miniprot alignment. Start codon fix: ORFs
   without an upstream in-frame stop are extended to a start codon hint of
   miniprothint. Isoforms whose CDS is a subsequence of another one are then
   collapsed.
5. A LightGBM model classifies each *ab initio* transcript as wrong, partial or
   correct, from miniprot alignment and miniprothint hint features. Transcripts
   with P(partial) + P(correct) >= 0.5 whose most likely class is `correct`
   are kept. The model (`drusilla_lgb_3class_v1`, a LightGBM text model) is
   downloaded by Nextflow at the start of the run and checked against its
   sha256, so the compute nodes need no internet for it.
6. Hint rescue: loci of `partial` transcripts without a kept transcript are
   predicted again by Tiberius with the hints of the best protein chain of the
   locus; loci without a protein chain are predicted *ab initio*, without hints.
   With `rescue_orf_filter: true`, loci where a Drusilla ORF already has all
   introns of the chain are skipped; the filter is off by default, and the
   F1 values above were measured without it. This needs a Tiberius with
   `--hints` (branch `hint_integration`, with bricks2marble from its branch
   `intron_hints`). The step therefore runs in its own image (see
   [Containers](containers.md)). `rescue_tiberius`
   is another `tiberius.py` to use instead. If the Tiberius has no `--hints`,
   the step is skipped with a warning in the Nextflow log, and no GPU task is
   started.
7. The kept and rescued *ab initio* transcripts and the Drusilla ORFs are merged
   into `<tool>_evidence.gff3`.

```yaml
drusilla:
  run: auto              # auto: vertebrate models only; true: always; false: never
  model: vertebrates     # released Drusilla model
  weights: null          # a local .weights.h5 instead (needs config)
  config: null
  cache_dir: null        # model cache (DRUSILLA_CACHE_DIR); null = downloaded once on the submitting host
  batch_size: null       # null = Drusilla sizes the batch from the GPU memory
  shards: 1              # parallel Drusilla processes; e.g. 24 on a 48-core CPU node
  min_coding_length: 200 # minimal CDS length of a Drusilla ORF
  lgb_model: https://bioinf.uni-greifswald.de/bioinf/drusilla/models/drusilla_lgb_3class_v1.tar.gz
                         # LightGBM model of the ab initio filter (required): the
                         # archive (URL or file), its unpacked directory, or a .txt with its .json
  lgb_model_sha256: d5bf3a9774914c6ca6421507c4cb956ef813b0bd2305bf45717af5f5eb42f1da
                         # sha256 of the archive; null for another model
  lgb_threshold: 0.5
  lgb_keep: correct      # classes kept after the threshold, e.g. "correct,partial"
  fix_stop: true         # stop codon fix of the ORFs
  fix_start: true        # start codon fix of the ORFs (needs fix_stop)
  rescue: true           # hint rescue of partial ab initio genes
  rescue_tiberius: null  # tiberius.py with --hints; null = the one of the hint rescue image
  rescue_model_cfg: null # model file or name; null = the Tiberius model of the run, vertebrates for Vipsania
  rescue_flank: 25000
  rescue_hint_weight: 2.5
  rescue_orf_filter: false # skip loci where a Drusilla ORF has all introns of the chain
  rescue_seq_len: 99990   # Tiberius seq_len of the rescue loci
  min_length: 300        # StringTie pre-filter
  min_cov: 3
  min_tpm: 1
  long_length: 3000
  min_tpm_long: 0.5
```

Drusilla runs on a GPU (label `gpu`) in the Drusilla image, which also runs
the LightGBM filter (`lightgbm` 4.7.0, `pyfaidx`, `pandas`). Without a
GPU, send `DRUSILLA_ANNOTATE` to CPU nodes in your site config and set
`shards`: one Drusilla process uses only one to two cores, so `shards`
splits the transcripts by gene into parts that run in parallel. On Bos taurus
24 shards used 26 of 48 cores and 40 GB.

All parameters of the flow are described in
[parameters.md](parameters.md#drusilla-flow).

## Known issues

- The LightGBM model was trained on Drusilla ORFs of 48 vertebrates (listed
  in the `.json` of the archive) and is applied to gene finder predictions.
- The hint rescue needs the Tiberius branch `hint_integration` and the
  bricks2marble branch `intron_hints`, which are not part of a release. It
  runs in a separate image until they are. With a site config that runs the
  rescue in the released Tiberius image, it is skipped with a warning.
- The scripts of the flow (`bin/filter_stringtie_gtf.py`,
  `compute_orf_features.py`, `apply_lgb_model_gtf.py`, `fix_stop_by_miniprot.py`,
  `prepare_hint_rescue_loci.py`, `filter_and_merge_rescue_gtf.py`,
  `chainedHints.py`) are copies from tiberius_orf_finder and Tiberius.
