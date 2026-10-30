# Parameters

All inputs and settings of the pipeline go into one YAML file. A commented
template is [conf/parameters.yaml](../conf/parameters.yaml). Pass the file to
the launcher:

```bash
paludamentum --params_yaml parameters.yaml --nf_config slurm_generic
```

`--nf_config` takes a path or the name of a config in `conf/`, see
[hpc.md](hpc.md). The launcher also accepts most parameters on the command
line (`paludamentum --help`); command line values override the file. The
launcher writes the merged parameters of a run to `<outdir>/params.yaml`.

## Input Data
### Genomic Sequences (Required Input)
Path to the genome FASTA file.

Example:
```yaml
genome: /path/to/genome.fa
```

### Protein Sequences (required for every mode except `abinitio`)
Path to one or more FASTA files with protein sequences aligned to the target genome.
Multiple files are concatenated in input order. Evidence runs (`rnaseq`,
`isoseq`, `mixed`) always need protein evidence: `combine_genesets.py`
reads four DIAMOND tables, miniprot alignments and miniprothint hints.

Example:
```yaml
proteins: /path/to/proteins.faa
```

Or a list:
```yaml
proteins: [
  /path/to/proteins1.faa,
  /path/to/proteins2.faa,
]
```

### ODB12 Partitions (Optional)
Provide one or more named partitions to download and merge into the protein evidence.
Available partitions: `Metazoa`, `Vertebrata`, `Viridiplantae`, `Arthropoda`, `Fungi`,
`Alveolata`, `Stramenopiles`, `Amoebozoa`, `Euglenozoa`, `Eukaryota`.

Example:
```yaml
odb12Partitions: [
  Metazoa,
  Fungi,
]
```

### RNA-Seq/Iso-Seq (Local)
The pipeline accepts local long and short read data. Add the absolute paths of your files to:
`rnaseq_single`, `rnaseq_paired`, `isoseq`.

The launcher expands `~` and environment variables and resolves relative
paths against the directory where it is launched, not the directory of
`params.yaml`; Nextflow itself does neither, so a params file that is passed
to `nextflow run` directly needs absolute paths. Paths must not contain spaces.

#### `rnaseq_paired`
Three forms are accepted:

1. **Single glob string** (most concise — matches one or many libraries):
   ```yaml
   rnaseq_paired: "/abs/path/to/rnaseq/*_{1,2}.fastq.gz"
   ```
   A list of glob strings is **not** supported; use one glob that covers
   all libraries, or switch to the explicit list form below.

2. **List of explicit `[r1, r2]` pairs** (one entry per library):
   ```yaml
   rnaseq_paired:
     - ["/abs/path/RNA/SRR37921248_1.fastq", "/abs/path/RNA/SRR37921248_2.fastq"]
     - ["/abs/path/RNA/SRR37921249_1.fastq", "/abs/path/RNA/SRR37921249_2.fastq"]
     - ["/abs/path/RNA/SRR37921250_1.fastq", "/abs/path/RNA/SRR37921250_2.fastq"]
     - ["/abs/path/RNA/SRR37921251_1.fastq", "/abs/path/RNA/SRR37921251_2.fastq"]
   ```

3. **Flat list of exactly two FASTQ files** (treated as a single library):
   ```yaml
   rnaseq_paired:
     - /abs/path/sample_R1.fastq.gz
     - /abs/path/sample_R2.fastq.gz
   ```

#### `rnaseq_single` and `isoseq`
Either a glob string or a YAML list of FASTQ files:
```yaml
rnaseq_single:
  - /abs/path/sample_a.fastq.gz
  - /abs/path/sample_b.fastq.gz
```


### RNA-Seq/Iso-Seq (SRA-download)
RNA-Seq and Iso-Seq libraries can also be automatically downloaded from the Sequence-Read-Archive by specifying their SRA IDs at
`rnaseq_sra_paired`, `rnaseq_sra_single`, `isoseq_sra`.  
An accession whose download fails twice (a network error, or a run that `fasterq-dump` cannot read) is dropped
with an `Error is ignored` note, as a library below `min_alignment_rate` is dropped; the other libraries go on,
and `-resume` tries the dropped accession again. The reads stay in the Nextflow work directory; with
`keep_downloads: true` they are also copied to `sra_downloads/` in `outdir`.


## Parameters
You can also set parameters of the pipeline within the file, default parameters that you can overwrite are:

| Parameter                  | Default Value                          | Description                                                                                                  |
| -------------------------- | -------------------------------------- | ------------------------------------------------------------------------------------------------------------ |
| `threads`                  | `48`                                   | Default number of CPUs reserved per task. Every tool runs with exactly the CPUs reserved for its task (`task.cpus`); a site config that sets `process.cpus` overrides this default. |
| `outdir`                   | `"results"` (the launcher sets `<genefinder>_results`) | Directory where all final results are written.                                                               |
| `scoring_matrix`           | `"conf/blosum62.csv"` | Amino acid substitution scoring matrix used by homology-based tools.                                         |
| `mode`                     | inferred                               | Pipeline mode, see [Mode](#mode).                                                                            |
| `keep_downloads`           | `false`                                | Copy the reads downloaded from SRA to `sra_downloads/` in `outdir`. Off, they are kept in the work directory only, which saves a second copy of every library. |
| `min_alignment_rate`       | `80`                                   | RNA-Seq and Iso-Seq libraries whose alignment rate (`samtools flagstat`, percent mapped) is below this value are dropped. If every RNA-Seq (or every Iso-Seq) library is dropped, the run stops with an error. |
| `strand_fix`               | `true`                                 | Recover the library strandedness of every short-read BAM before merging and StringTie (`STRAND_FIX_BAM`, `drusilla strand-fix-bam`): the protocol of each run is inferred from its spliced reads, and the `XS` tag is set on all reads of stranded runs. Reads of unstranded runs are unchanged. Not applied to Iso-Seq BAMs; pyVARUS directories bring no BAM. The per-BAM report is `intermediate/strand_fix/<bam>.strand_runs.tsv`. Set `false` to skip it. |
| `strand_fix_run_regex`     | none                                   | Regex whose group 1 is the run id in the read names of a user BAM (`rnaseq_bam`) that pools several runs, e.g. from VARUS. Without it Drusilla's default is used (leading letters and digits, e.g. `SRR123456`). A HISAT2 BAM of the pipeline is always one run. |
| `hc_table`                 | `"conf/hc_genes.yaml"`                 | Clade table for the combined HC flow: `drusilla_model` for each clade, insect clades trigger `--drop-lone-single-exon tiberius`, `recommended` of a clade warns a run that differs from it. Today Vertebrata, Mammalia (`drusilla_model: vertebrates`), Insecta, Lepidoptera (`drusilla_model: insects`). |
| `td2_predict_args`         | none                                   | Options appended to `TD2.Predict`, for example `"--precise"`.                                                |
| `combine.min_locus`        | `4`                                    | Minimum tier of the best model of a locus for `combine_genesets.py` to keep the locus. Everything else uses the script defaults; see `combine_genesets.py --help`. |

The location of the required executables is set by default so that they are available in your path, unless you are using the Singularity container:
| Tool Parameter                       | Default Value                           | Description                                                         |
| ------------------------------------ | --------------------------------------- | ------------------------------------------------------------------- |
| `tools.hisat2`                       | `"hisat2"`                              | Path or module name for HISAT2 aligner.                             |
| `tools.hisat2_build`                 | `"hisat2-build"`                        | Path to HISAT2 genome index builder.                                |
| `tools.minimap2`                     | `"minimap2"`                            | Path to Minimap2 aligner (used for Iso-Seq).                        |
| `tools.stringtie`                    | `"stringtie"`                           | Path to StringTie3 for transcript assembly.                          |
| `tools.samtools`                     | `"samtools"`                            | Path to Samtools for BAM/FASTQ operations.                          |
| `tools.transdecoder_util_gtf2fa`     | `"gtf_genome_to_cdna_fasta.pl"`         | Utility for converting GTF + genome → cDNA FASTA.                   |
| `tools.transdecoder_util_orf2genome` | `"cdna_alignment_orf_to_genome_orf.pl"` | Maps the revised TD2 ORFs from transcript to genome coordinates.    |
| `tools.transdecoder_gtf2gff`         | `"gtf_to_alignment_gff3.pl"`            | Converts GTF to alignment-style GFF3 (writes the GFF3 that TD2 needs from the StringTie assembly). |
| `tools.td2_longorfs`                 | `"TD2.LongOrfs"`                        | Long ORF detection of TD2.                                          |
| `tools.td2_predict`                  | `"TD2.Predict"`                         | ORF prediction of TD2.                                              |
| `tools.diamond`                      | `"diamond"`                             | DIAMOND, used both for species selection (12-col) and the four combine_genesets runs (14-col). |
| `tools.miniprot`                     | `"miniprot"`                            | miniprot: protein-to-genome aligner.                                |
| `tools.miniprot_boundary_scorer`     | `"miniprot_boundary_scorer"`            | miniprot boundary scorer for boundary refinement.                   |
| `tools.miniprothint`                 | `"miniprothint.py"`                     | miniprothint wrapper for protein hints.                             |
| `tools.bam2hints`                    | `"bam2hints"`                           | Converts RNA-Seq BAM alignments into intron hints.                  |
| `tools.bedtools`                     | `"bedtools"`                            | `bedtools genomecov -bg -split` writes `rna_cov.bg.gz` for combine_genesets. |



### Aligned RNA-Seq reads (BAM)
Already aligned short reads, for example from VARUS, go to `rnaseq_bam` as a
list of BAM files. They skip the HISAT2 step. A BAM whose header does not say
`SO:coordinate` is sorted first; a coordinate-sorted one is used as it is.
Like the HISAT2 BAMs, user BAMs go through the strand fix (`strand_fix`); the
run ids are taken from the read names (`strand_fix_run_regex`).

### pyVARUS output directories
[pyVARUS](https://github.com/Gaius-Augustus/pyVARUS) writes `stringtie.gtf`
and `hints.gff` next to `VARUS.bam`, made as the pipeline makes them from a
BAM. Pass the output directory instead of the BAM; this also works after
`varus run --drop-bam`. The directories skip mapping, `bam2hints` and
StringTie. A pyVARUS directory brings no BAM: with pyVARUS input only,
there is no coverage track, and `combine_genesets.py` runs without
`--coverage` (the intron hints are used).

```yaml
rnaseq_varus:                # short-read runs
  - /path/to/varus_short
isoseq_varus:                # --longreads runs
  - /path/to/varus_long
```

`rnaseq_varus` counts as short reads and `isoseq_varus` as Iso-Seq when the
mode is inferred. The pipeline checks that the run used the same genome (MD5
in `VARUS.manifest.tsv`) and the right read type. A directory of an older
pyVARUS without `stringtie.gtf` needs
`varus assemble genome.fa --short DIR/VARUS.bam --outdir DIR` first (`--long`
for long reads).

The hints of several sources of one read type (several directories, or
directories and reads) are summed as `bam2hints` would count them on all
reads. In mixed mode both hint files (short reads and Iso-Seq) go to
`combine_genesets.py`, which sums them per intron. The combined HC flow
needs one assembly of all reads of a read type: in `rnaseq` and `isoseq`
mode pass exactly one directory and no other reads of that type. In `mixed`
mode it needs one `stringtie --mix` assembly of both read types:

```sh
varus assemble genome.fa --short varus_short/VARUS.bam --long varus_long/VARUS.bam --outdir varus_mixed
```

```yaml
rnaseq_varus: /path/to/varus_short
isoseq_varus: /path/to/varus_long
mixed_varus:  /path/to/varus_mixed
```

Then pass no other short reads or Iso-Seq. `mixed_varus` is used only in
mixed mode. `varus replay` rebuilds a BAM deleted by
`--drop-bam`.

### StringTie assemblies

A StringTie assembly made outside the pipeline goes to `stringtie`
(`--stringtie`): exactly one GTF; several files are an error. It replaces the
assembly. It can be combined with reads, BAM files and Iso-Seq reads: these
are mapped and merged into one BAM only, nothing is assembled. The BAM gives
the coverage track and the intron hints.

```yaml
stringtie: /path/to/stringtie.gtf
```

`stringtie` counts as short reads when the mode is inferred, so proteins and
an assembly give the mode `rnaseq`; a forced `isoseq` or `mixed` accepts it
as well. The combined HC flow takes exactly one assembly and uses it as it
is (no `stringtie --merge`). It cannot be combined with pyVARUS directories
(`rnaseq_varus`, `isoseq_varus`, `mixed_varus`): that is an error. The GTF
needs the `cov` and `TPM` attributes that StringTie writes when it assembles
reads; an assembly without them (the output of `stringtie --merge`) stops
the run. A missing file is an error. Without reads there is no BAM, so no
coverage track and no intron hints: `combine_genesets.py` runs without
`--coverage` and `--introns`, its RNA rules stay off, and `hintsfile.gff`
holds protein hints only.

## Gene finder: Tiberius

In evidence modes Tiberius always runs; set `tiberius.model_cfg` to a
model configuration. In `abinitio` / `proteins` mode set `tiberius.run:
true` as well (one gene finder runs there). The available configurations
are the YAML files in `tiberius/model_cfg/` of the Tiberius submodule (also
listed in
[Tiberius/model_cfg](https://github.com/Gaius-Augustus/Tiberius/tree/main/model_cfg)).

| Parameter | Default | Description |
| --- | --- | --- |
| `tiberius.model_cfg` | none | Name of a Tiberius model configuration, or path to a configuration file. Required in evidence modes (unless `tiberius.result` is set). |
| `tiberius.run` | `false` | Only consulted in `abinitio` / `proteins` modes: `true` selects Tiberius as the single gene finder (equivalent to `genefinder: tiberius`). Ignored in evidence modes, which always run both finders. |
| `tiberius.model_dir` | none | Directory that holds the extracted weights directory (the archive name of `weights_url` without `.tar.gz`, e.g. `vertebrates_weights`), instead of the download at the start of the run. Tiberius then runs with `--model` and never downloads. |
| `tiberius.result` | none | Existing Tiberius prediction (GFF3 or GTF). Reused instead of running Tiberius. In evidence modes `tiberius.model_cfg` can then be omitted; the clade comes from `tiberius.model_cfg` (if given) via `target_species`. |
| `tiberius.min_split_size` | `20000000` | Minimal size in bp of a genome chunk. |
| `tiberius.max_files` | `20` | Maximal number of genome chunks, which is the upper limit of parallel Tiberius tasks. |
| `tiberius.max_parallel` | unlimited | Cap of concurrently running Tiberius tasks, e.g. `1` on a single-GPU workstation. |
| `tiberius.batch_size` | automatic | Forwarded to `tiberius.py --batch_size`. |
| `tiberius.seq_len` | automatic | Forwarded to `tiberius.py --seq_len`. |

## Gene finder: Vipsania

In evidence modes Vipsania always runs; set `vipsania.model` and (default)
leave `vipsania.finetune: true`. In `abinitio` / `proteins` mode set
`vipsania.run: true` to pick Vipsania as the single gene finder (or use
`genefinder: vipsania`). All Vipsania-specific parameters are in
[vipsania.md](vipsania.md).

## Drusilla

Drusilla is one of the two ORF callers of the combined HC flow (step 4),
alongside TD2 (see [drusilla_flow.md](drusilla_flow.md)). The parameters
below configure the Drusilla model, the StringTie pre-filter and the
Drusilla GPU task.

| Parameter | Default | Description |
| --- | --- | --- |
| `drusilla.model` | `drusilla_model` of the clade in `hc_table` | Released Drusilla model (`drusilla models list`): `vertebrates` or `insects`. |
| `drusilla.result` | none | Existing Drusilla ORFs GTF: used as the Drusilla gene set instead of running `FILTER_STRINGTIE`, `DRUSILLA_ANNOTATE` and `COLLAPSE_ORFS`. The file is passed to `combine_genesets.py` as the `drusilla` input. |
| `drusilla.weights`, `drusilla.config` | none | A local `.weights.h5` file and its architecture YAML instead of a released model. |
| `drusilla.cache_dir` | none | Drusilla model cache (`DRUSILLA_CACHE_DIR`), e.g. filled by `drusilla models download vertebrates`. Without it the archive of the model's `weights_url` in `drusilla/model_cfg/` is downloaded once on the submitting host (label `download`), checked against its `weights_sha256` and staged into the Drusilla task. |
| `drusilla.batch_size` | automatic | Drusilla batch size. |
| `drusilla.shards` | `1` | Parallel Drusilla processes, each on a part of the transcripts. For CPU nodes, where one process uses only one to two cores (e.g. 24 with 48 CPUs). |
| `drusilla.min_coding_length` | `200` | Minimal CDS length of a Drusilla ORF. |
| `drusilla.min_length`, `min_cov`, `min_tpm` | `300`, `3`, `1` | StringTie pre-filter: transcript length, coverage and TPM. |
| `drusilla.long_length`, `min_tpm_long` | `3000`, `0.5` | Relaxed TPM for transcripts of at least `long_length` nt. |
| `drusilla.drop_unstranded` | `true` | StringTie pre-filter: drop transcripts with strand `.`. |
| `drusilla.unstranded_both_strands` | `false` | `--unstranded-both-strands`: also scan the reverse strand of transcripts with strand `.`. Needs `drop_unstranded: false`; otherwise the run warns and the option has no effect. |
| `drusilla.flank_bp` | none | Genomic context in bp on both sides of every transcript, with ORF completion into it (`--flank-bp N --flank-3p-bp M --flank-keep`). `null` or `0` = off. |
| `drusilla.flank_3p_bp` | `flank_bp` | Context on the 3' side alone (`--flank-3p-bp`); enough on its own to turn the completion on. |
| `drusilla.flank_min_tx_frac` | Drusilla default | `--flank-min-tx-frac`: minimal fraction of an ORF inside the transcript. Only passed with `flank_bp` or `flank_3p_bp`. |
| `drusilla.orfs_per_tx` | Drusilla default (`1`) | `--orfs-per-tx`: ORFs per transcript, longest first; `0` = all. |

## Post-processing and quality control

The steps after the integration, see [postprocessing.md](postprocessing.md).
Command line: `--busco_lineage`, `--reference_annotation`, `--ncrna`,
`--fantasia`; the other keys go into the params file.

| Parameter | Default | Description |
| --- | --- | --- |
| `postprocess.sanity_filter` | `true` | Remove transcripts with an internal stop codon, a CDS length that is not a multiple of 3, a broken CDS structure or no CDS; `false`: report only (`qc/sanity_filter.tsv`). |
| `postprocess.utr` | `true` | UTRs from the StringTie assemblies for transcripts without UTRs (modes `rnaseq`, `isoseq`, `mixed`). |
| `postprocess.max_utr_extension` | `5000` | Maximal UTR extension in bp beyond the CDS. |
| `qc.busco_lineage` | none | BUSCO/compleasm lineage, e.g. `eukaryota_odb12`, `vertebrata_odb12` (without `_odbNN`: odb12). Without it there is no completeness assessment. |
| `qc.busco_download_path` | `~/.cache/paludamentum/busco` | Lineage cache, layout `<dir>/lineages/<lineage>/`; a missing lineage is downloaded once on the submitting host. |
| `qc.compleasm`, `qc.busco` | `true` | compleasm and BUSCO on genome and proteome (with `qc.busco_lineage`). |
| `qc.omark` | `false` | OMArk on the proteome; needs `qc.omamer_db`. |
| `qc.omamer_db` | none | OMAmer database (e.g. `LUCA.h5`, about 15 GB). |
| `qc.ete_taxa_path` | `~/.cache/paludamentum/ncbi_taxonomy` | NCBI taxonomy for OMArk: ete3's `taxa.sqlite`, built once on the submitting host from `taxdump.tar.gz`, which is downloaded there unless the directory holds it. |
| `qc.reference_annotation` | none | Reference annotation (GFF3 or GTF, gzipped or not): `gffcompare` at CDS level, `qc/gffcompare.stats`. |
| `qc.gene_support` | `true` | Hint support per transcript, `qc/gene_support.tsv` (every mode but `abinitio`). |
| `qc.statistics` | `true` | Gene set statistics and plots. |
| `qc.report` | `true` | `report.html`. |
| `ncrna.run` | `false` | rRNA, tRNA, Rfam and lncRNA genes, `<stem>_with_ncRNA.gff3`. |
| `ncrna.rfam_dir` | `~/.cache/paludamentum/rfam` | Directory with `Rfam.cm` and `Rfam.clanin` (Rfam 15.1); unset: downloaded once on the submitting host. |
| `ncrna.trnascan_high_confidence` | `false` | tRNAscan-SE high-confidence filter (EukHighConfidenceFilter). |
| `ncrna.lncrna` | `true` | lncRNA with FEELnc (modes with transcripts). |
| `fantasia.run` | `false` | GO terms with FANTASIA-Lite (GPU), `<stem>_go.gff3`. |
| `fantasia.hf_cache_dir` | none | Hugging Face cache with `Rostlab/prot_t5_xl_uniref50` (required with `fantasia.run`). |
| `fantasia.lookup_dir` | none | `lookup_table.npz`, `annotations.json`, `accessions.json` of Zenodo record 17720428 (required with `fantasia.run`). |
| `fantasia.min_score` | `0.5` | Minimum score of a GO term. |
| `fantasia.device` | `cuda` | `cuda`, or `cpu`: the ProtT5 embeddings on the CPU, no GPU and no GPU probe (much slower; route `FANTASIA_ANNOTATE` to a CPU queue, see [postprocessing.md](postprocessing.md#go-terms-with-fantasia-lite)). |
| `fantasia.proteins_per_task` | none | Proteins per `FANTASIA_ANNOTATE` task; the tasks run in parallel and `FANTASIA_MERGE` joins their results. None: 2000 with `fantasia.device: cpu`, all proteins in one task with `cuda`; `0`: one task. |
| `fantasia.additional_params` | none | Appended to `fantasia_pipeline.py`. |
| `tools.gffread`, `tools.compleasm`, `tools.busco`, `tools.gffcompare`, `tools.trnascan`, `tools.cmscan`, `tools.cmpress`, `tools.barrnap`, `tools.omamer`, `tools.omark` | the names on `PATH` in the images | Executables of the post-processing (`--check_tools` checks those of the steps that are on). GenomeTools (`gt`, GFF3 validation of the final files) has no override and is called by name. |

## Mode

The pipeline infers its mode from the inputs, see the table in the
[README](../README.md#inputs-and-modes). Set `mode` to force one of
`abinitio`, `proteins`, `rnaseq`, `isoseq`, `mixed`. `tiberius` is accepted as
the historic name of `abinitio`; any other value is an error, as is a forced
mode whose inputs are missing (for example `rnaseq` without reads or a StringTie assembly).
