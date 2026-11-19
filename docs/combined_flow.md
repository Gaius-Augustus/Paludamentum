# Evidence and integration steps (combined HC flow)

This page describes steps 1 to 4 of the [overview](../README.md#what-the-pipeline-does)
in technical detail: the tools and options of each step, the integration with
`bin/combine_genesets.py`, the selection models of the modes and the clade
table. Step 5 is described in [postprocessing.md](postprocessing.md).

<p align="center">
  <img src="../figures/combined_flow.svg" alt="Combined HC-gene step: Tiberius and Vipsania predict genes on the genome; StringTie assembles one transcript set from the reads; TD2 (with shorten + classify + revise) and Drusilla predict ORFs on it; DIAMOND runs on all four protein sets; bedtools genomecov writes the coverage, bam2hints the introns; combine_genesets.py integrates everything into combined.gff3" width="100%">
</p>

## Steps of a run with evidence

1. ***Ab initio* predictions.** Both gene finders run in parallel: Tiberius
   on genome chunks, and Vipsania with `--finetune` on the whole genome
   (default on; one GPU task). An existing Tiberius or Vipsania prediction
   (`--result`) is used instead of running that gene finder, see
   [Existing predictions](../README.md#existing-predictions).
2. **Protein evidence.** Proteins (FASTA files and/or OrthoDB v12
   partitions) are aligned with miniprot, scored with
   miniprot-boundary-scorer, and converted to hints with miniprothint. For
   very large protein databases the Tiberius *ab initio* proteins are used to
   select the most relevant source species first (DIAMOND).
3. **Transcript evidence.** Short reads are aligned with HISAT2, Iso-Seq reads
   with minimap2. Libraries with a low alignment rate are dropped. The
   library strandedness of every short-read BAM is recovered before the
   assembly (`strand_fix`, on by default). One StringTie assembly of all
   reads is made (`stringtie`, `stringtie -L` or `stringtie --mix`). An
   existing StringTie assembly (`--stringtie`) is used instead, see
   [Existing StringTie assemblies](../README.md#existing-stringtie-assemblies).
   Each read type can also come as BAM files (`rnaseq_bam`, `isoseq_bam`)
   or as its own assembly (`rnaseq_stringtie`, `isoseq_stringtie`); a mixed
   run with a given assembly uses the short-read and the Iso-Seq assembly
   side by side, see [parameters.md](parameters.md#stringtie-assemblies).
4. **Combined high-confidence genes.** Two ORF callers run on the StringTie
   assembly: TD2 (with `TD2.LongOrfs`, `TD2.Predict`, then two auxiliary
   DIAMOND passes and `revise_pep.py` for the partial-ORF rescue; the
   revised ORFs are mapped back to the genome with
   `cdna_alignment_orf_to_genome_orf.pl`) and Drusilla (StringTie pre-filter
   and `drusilla annotate` with the Drusilla model of all clades, see
   [drusilla_flow.md](drusilla_flow.md)). DIAMOND runs on the proteins of
   all four gene sets (Tiberius, Vipsania, revised TD2, Drusilla).
   `bedtools genomecov -bg -split` writes a read-depth bedGraph, and
   `bam2hints --intronsonly` plus `filterIntronsFindStrand.pl` writes intron
   hints. `bin/combine_genesets.py` integrates the four gene sets, the four
   DIAMOND tables, miniprot alignments and miniprothint hints, the bedGraph
   and the introns into one GFF3 (`intermediate/combined.gff3`), which the
   post-processing turns into `paludamentum.gff3`.

The same flow runs for every clade. The clade of the gene finder model only
decides two switches of `combine_genesets.py`, see
[Clade table](#clade-table).

## Integration step by step

Steps of the combined HC flow (the default run; the other modes leave out the inputs they do not have):

1. One StringTie assembly of all reads (short, Iso-Seq with `-L`, or both
   with `--mix`). A `--stringtie` assembly made outside the pipeline is
   used as it is; it must carry the `cov` and `TPM` attributes. The one
   assembly is not passed through `stringtie --merge` (only the GFF3 that
   TD2 needs is written from it); several assemblies would be merged.
2. **Drusilla** filters the StringTie assembly (length, coverage, TPM) and
   predicts the ORFs of the kept transcripts with the one Drusilla model of
   all clades (`drusilla_model`, bundled with Drusilla). Isoforms whose CDS
   is a subsequence of another one are collapsed.
3. **TD2** predicts ORFs on the full StringTie assembly. ORFs are then
   shortened to their first methionine where incomplete,
   `CDS_CLASSIFY_AND_REVISE` runs two auxiliary DIAMOND passes to decide
   whether to keep the original or the trimmed start, and the revised ORFs
   are mapped back to the genome with
   `cdna_alignment_orf_to_genome_orf.pl`.
4. For each of the four gene sets (Tiberius, Vipsania, revised TD2,
   Drusilla), `PROTEIN_FROM_GFF` extracts proteins with gffread (`-S`:
   internal stops are written as `*`, which DIAMOND accepts).
   `bin/check_internal_stops.py` writes a `WARN` line to the run log when
   more than 5% of a gene set's proteins have an internal stop. The usual
   cause is wrong CDS phases in a reused `<tool>.result`, which
   `bin/fix_cds_phases.py` repairs; the run continues.
   `DIAMOND_BLASTP_FULL` searches the proteins against the protein database with
   the 14-column outfmt (`qseqid sseqid pident length mismatch gapopen
   qstart qend sstart send evalue bitscore qlen slen`). It reads GFF3 and
   GTF (the Drusilla ORFs) and keeps at most 10 isoforms per gene; an
   isoform beyond the cap has no DIAMOND hit. The same step drops malformed
   transcripts (CDS parts of one transcript that overlap, or a feature with
   start beyond end) and names them in a warning on stderr, because gffread
   0.12.7 corrupts its heap on them. This affects only the protein extraction:
   the gene set handed to the combiner is unchanged, so such a transcript has
   no protein and no DIAMOND hit.
5. `bedtools genomecov -ibam BAM -bg -split | gzip > rna_cov.bg.gz` writes
   the read-depth bedGraph, and `bam2hints --intronsonly` plus
   `filterIntronsFindStrand.pl` writes the intron hints. In mixed mode both
   intron hint files (short reads and Iso-Seq) go to `combine_genesets.py`,
   which sums them per intron. pyVARUS directories and a `--stringtie`
   assembly without reads have no BAM: `combine_genesets.py` then runs
   without `--coverage` (pyVARUS hints are still used). Transcripts of
   full-length reads (`chains`, or in mixed mode the `stringtie -L` assembly
   of the Iso-Seq reads alone) go to `combine_genesets.py --chains`; the ORF
   callers do not see them.
6. `combine_genesets.py` integrates all of the above. By default
   (`combine.classifier = "${projectDir}/conf/combiner_classifier_v1.json"`)
   it runs in classifier mode: the candidate pool is built as with
   `--keep-all` (every unique CDS structure, all tier rules off), and a
   gradient-boosted tree model (400 trees of depth 3, exported to JSON,
   evaluated with numpy; scikit-learn is not needed at run time) gives each
   candidate P(structure = reference CDS) from 47 features of the report
   (tools, structure, locus context, RNA-Seq intron and coverage evidence,
   DIAMOND and miniprot protein evidence; no clade input). Per locus the
   best candidate is the primary if its P reaches the locus threshold
   (0.25, or 0.10 for a fragmented assembly, see
   [Fragmented assemblies](#fragmented-assemblies)); a locus below the
   threshold is still kept if a multi-exon candidate has a full-length
   protein hit and protein hints for its introns (see
   [Protein rescue](#protein-rescue)); further candidates are
   isoforms if P reaches 0.30. The `reason`
   column of the report is `primary`, `primary_protein_rescue`,
   `isoform_classifier`,
   `isoform_rejected`, `locus_below_min` or `dropped_lone_single_exon`, and
   the last column is `p_classifier`. With `combine.classifier = false` the
   tier rules apply (merging, clustering, tiering with the four tool sets,
   RNA / protein / DIAMOND adjustments; `combine.min_locus` is the minimum
   tier of a locus); see [drusilla_flow.md](drusilla_flow.md) and
   `combine_genesets.py --help`. The full-length transcript evidence
   (`--chains`, mixed and Iso-Seq runs) is not a model feature (no training
   species has long reads); in classifier mode it raises the log-odds of a
   structure whose CDS chain lies in one full-length transcript by 1.5
   (`--chain-logit-full`).

## Selection models of the modes

All selection models are gradient-boosted trees of the same recipe, trained
on the keep-all candidate pools of the same 15 species (4 insects, 7 fungi,
2 vertebrates, 2 plants) with gffcompare labels. They differ in the features
that the inputs of the mode provide. The models are clade-agnostic.

| Run | Parameter | File | Features | Thresholds (locus / isoform) |
| --- | --- | --- | --- | --- |
| default run and modes `isoseq`, `mixed` | `combine.classifier` | `conf/combiner_classifier_v1.json` | 47: tools, structure, locus context, RNA-Seq intron and coverage evidence, DIAMOND and miniprot protein evidence | 0.25 / 0.30 (locus 0.30 for a clade outside the clade table; 0.10 on a fragmented assembly, see below); protein rescue of loci below the locus threshold, see below |
| mode `proteins` | `combine.classifier_proteins` | `conf/combiner_classifier_protonly_v1.json` | 41: no RNA-Seq features | 0.20 / 0.30 (locus 0.30 for a clade outside the clade table) |
| evidence mode without a Tiberius model | `combine.classifier_vipsania` | `conf/combiner_classifier_vipsania_v1.json` | as the default, without the Tiberius gene set | 0.225 / 0.30 |
| mode `abinitio` with both models | `combine.classifier_abinitio` | `conf/combiner_classifier_abinitio_v1.json` | 18: tools, structure and locus context only | 0.30 / 0.25 |

- **Mode `proteins`.** `combine_genesets.py` receives the Tiberius and
  Vipsania gene sets, their two DIAMOND tables and the miniprot and
  miniprothint evidence: no StringTie assembly, no ORF callers, no coverage
  track and no intron hints. `combine.classifier_proteins = false` selects
  with the tier rules, whose RNA rules stay off. `combine.t_locus` sets
  another locus threshold.
- **Evidence mode without a Tiberius model.** Vipsania (fine-tuned on the
  genome by default) is the only *ab initio* gene set, its proteins select
  the species of a large protein database, and `combine_genesets.py`
  integrates three gene sets (Vipsania, revised TD2, Drusilla). The RNA-Seq
  multiplicity is normalised by the median over the structures that all
  tools of the run predicted (three here, four with Tiberius).
- **Mode `abinitio` with both models.** `combine_genesets.py` integrates the
  two predictions without any evidence: no protein database, no DIAMOND, no
  miniprot, no transcripts.

## Fragmented assemblies

In the evidence-mode classifier run with a Tiberius model,
`combine_genesets.py` reads the genome (`--genome`), computes the number of
sequences, length, N50 and L50 (the number of the longest sequences that hold
half of the genome) and, if the L50 exceeds `combine.fragmented_l50` (100),
uses `combine.t_locus_fragmented` (0.10) as locus threshold. Many genes of a
fragmented assembly have only partial models whose low P would cost the whole
locus. Precedence: `combine.t_locus`, then the fragmented threshold, then the
0.30 of a clade new to the classifier (`--unknown-clade`), then the model's.
`combine.fragmented_l50 = null` switches the rule off. The isoform threshold
is unchanged. The other modes and models do not use the rule.

Evidence (Nov 2026, replay of the selection against NCBI RefSeq): two
fragmented assemblies (L50 639 and 1939) gained +0.11 and +0.08 gene-level F1
at 0.10, with exact-CDS F1 unchanged (+0.008, +0.005). Fifteen contiguous
assemblies (L50 <= 33) gain little from a lower threshold (at most +0.03 at
0.15) and lose exact-CDS F1 at 0.10 (-0.01 to -0.07). No benchmark assembly
has an L50 above 56.

## Protein rescue

The classifier's P estimates an exact match to the reference CDS. A conserved
gene whose only model is imperfect (typically predicted by one gene finder
alone, with partial RNA-Seq support) can therefore stay below the locus
threshold and be missing from the gene set, which shows as lower BUSCO
completeness. In classifier mode, a locus whose best candidate stays below the
locus threshold is kept all the same if one of its candidates has

- at least two CDS exons,
- a best DIAMOND hit that covers at least `combine.rescue_protein_cov` (0.9) of
  the candidate and of the protein (report columns `dmd_qcov`, `dmd_tcov`), and
- protein hints for at least `combine.rescue_protein_introns` (0.5) of its
  introns (report column `prot_intron_frac`).

The candidate with the highest P among those becomes the primary transcript
(`reason` = `primary_protein_rescue`); the other candidates of the locus follow
the isoform threshold as usual. Missing evidence never qualifies a candidate.
Genomes with few introns per gene gain little, because the rule needs a
multi-exon candidate.

Scope: the same as the fragmented-assembly threshold, the classifier of the
evidence modes with Tiberius. Not used in mode `proteins`, with the
Vipsania-only classifier, in the ab initio combination or with the tier rules.
The rule acts in addition to whichever locus threshold applies (`combine.t_locus`,
fragmented, unknown clade or the model's). `combine.rescue_protein_cov = null`
switches it off. The script `combine_genesets.py` has it off by default
(`--rescue-protein-cov C`, `--rescue-protein-introns F`); the pipeline passes
both options.

Evidence (Nov 2026, 36 benchmark species, C = 0.9, F = 0.5): BUSCO
completeness +1.1 points on average (vertebrates +1.5); exact-match gene-level
F1 -0.3 (sensitivity +0.5, precision -1.2); exon-level F1 +1.7. Without the
intron condition the rule cost more precision (gene F1 -0.8), because
single-exon candidates with a full-length hit that overlap no reference gene
were let in. An intron share of 0.75 kept less of the BUSCO gain (+0.7).

`intermediate/combined.settings.tsv` (`key`, `value`) records what a run
used: `selection`, `assembly_sequences`, `assembly_length`, `assembly_n50`,
`assembly_l50`, `fragmented_l50`, `fragmented` (0/1) and, in classifier mode,
`t_locus`, `t_locus_source` (`t_locus`, `fragmented_assembly`,
`model_unknown_clade` or `model`) and `t_iso`, then `rescue_protein_cov` and
`rescue_protein_introns` (empty when the rescue is off) and `rescued_loci` (the
number of loci with `primary_protein_rescue`); with the tier rules
`min_locus` instead.

## Clade table

`conf/hc_genes.yaml` (`hc_table`) lists the known clades of the gene finder
models. It does not select a workflow: the flow is the same for every clade.
It decides two things.

- `classifier_clade` names the training clade of the selection model that
  covers the clade: `Vertebrata` and `Mammalia` use `vertebrates`,
  `Mesangiospermae` (the Tiberius angiosperms model) uses `embryophyta`,
  `Fungi` uses `fungi`, `Insecta` and `Lepidoptera` use `insects`. A clade
  with neither `classifier_clade` nor `standard_thresholds`, or one not in the
  table, is new to the model and the pipeline passes `--unknown-clade`, which
  raises the locus threshold of the default model to 0.30.
- `standard_thresholds: true` marks a clade the classifier has not seen but
  that gets the thresholds of the training clades (0.25 with the evidence-mode
  model) instead of `--unknown-clade`. The Tiberius chlorophyta and diatoms
  models (Chlorophyta, Bacillariophyta) use it: on two chlorophytes and two
  diatoms it was as accurate as 0.30 or slightly better. A Vipsania-only
  diatom run (Vipsania model Stramenopiles, no Tiberius model) still counts
  as an unknown clade, since Stramenopiles covers more than diatoms.
- For `Insecta` and `Lepidoptera` the pipeline passes
  `--drop-lone-single-exon tiberius` to `combine_genesets.py`, because
  Tiberius overpredicts single-exon genes in insects. In classifier mode
  lone single-exon Tiberius structures without a miniprot match are never
  selected (`dropped_lone_single_exon`); in rules mode they are demoted.

A new clade is added to the table with its `classifier_clade` (or with
`standard_thresholds: true`). The Drusilla
model does not depend on the clade (`drusilla.model`, default
`drusilla_model`).
