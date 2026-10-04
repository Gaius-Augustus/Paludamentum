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
Multiple files are concatenated in input order. RNA-Seq and Iso-Seq evidence
need protein evidence, because the high-confidence genes are selected by
protein homology.

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


## Parameters
You can also set parameters of the pipeline within the file, default parameters that you can overwrite are:

| Parameter                  | Default Value                          | Description                                                                                                  |
| -------------------------- | -------------------------------------- | ------------------------------------------------------------------------------------------------------------ |
| `threads`                  | `48`                                   | Default number of CPUs reserved per task. Every tool runs with exactly the CPUs reserved for its task (`task.cpus`); a site config that sets `process.cpus` overrides this default. |
| `outdir`                   | `"results"` (the launcher sets `<genefinder>_results`) | Directory where all final results are written.                                                               |
| `scoring_matrix`           | `"conf/blosum62.csv"` | Amino acid substitution scoring matrix used by homology-based tools.                                         |
| `mode`                     | inferred                               | Pipeline mode, see [Mode](#mode).                                                                            |
| `min_alignment_rate`       | `80`                                   | RNA-Seq and Iso-Seq libraries whose alignment rate (`samtools flagstat`, percent mapped) is below this value are dropped. |
| `transdecoder`             | `"td1"`                                | ORF finder of the HC gene step: `td1` (TransDecoder 5.7.1) or `td2` ([TD2](https://github.com/Markusjsommer/TD2), experimental). TD2 is not in the container image and must be on the `PATH` of the task. See [orf_finder_comparison.md](orf_finder_comparison.md). |
| `td2_predict_args`         | none                                   | Options appended to `TD2.Predict`, for example `"--precise"`.                                                |

The location of the required executables is set by default so that they are available in your path, unless you are using the Singularity container:
| Tool Parameter                       | Default Value                           | Description                                                         |
| ------------------------------------ | --------------------------------------- | ------------------------------------------------------------------- |
| `tools.hisat2`                       | `"hisat2"`                              | Path or module name for HISAT2 aligner.                             |
| `tools.hisat2_build`                 | `"hisat2-build"`                        | Path to HISAT2 genome index builder.                                |
| `tools.minimap2`                     | `"minimap2"`                            | Path to Minimap2 aligner (used for Iso-Seq).                        |
| `tools.stringtie`                    | `"stringtie"`                           | Path to StringTie3 for transcript assembly.                          |
| `tools.samtools`                     | `"samtools"`                            | Path to Samtools for BAM/FASTQ operations.                          |
| `tools.transdecoder_longorfs`        | `"TransDecoder.LongOrfs"`               | Tool to detect long ORFs for protein prediction.                    |
| `tools.transdecoder_predict`         | `"TransDecoder.Predict"`                | Main TransDecoder protein prediction step.                          |
| `tools.transdecoder_util_gtf2fa`     | `"gtf_genome_to_cdna_fasta.pl"`         | Utility for converting GTF + genome → cDNA FASTA.                   |
| `tools.transdecoder_util_orf2genome` | `"cdna_alignment_orf_to_genome_orf.pl"` | Maps ORF predictions from cDNA-level to genome coordinates.         |
| `tools.transdecoder_gtf2gff`         | `"gtf_to_alignment_gff3.pl"`            | Converts GTF to alignment-style GFF3 for TransDecoder.              |
| `tools.td2_longorfs`                 | `"TD2.LongOrfs"`                        | Long ORF detection of TD2 (`transdecoder: td2`).                    |
| `tools.td2_predict`                  | `"TD2.Predict"`                         | ORF prediction of TD2 (`transdecoder: td2`).                        |
| `tools.diamond`                      | `"diamond"`                             | Diamond alignment tool for protein-to-genome searches.              |
| `tools.bedtools`                     | `"bedtools"`                            | Bedtools for genomic interval operations.                           |
| `tools.miniprot`                     | `"miniprot"`                            | MiniProt — protein-to-genome aligner.                               |
| `tools.miniprot_boundary_scorer`     | `"miniprot_boundary_scorer"`            | Additional MiniProt scoring utility for boundary refinement.        |
| `tools.miniprothint`                 | `"miniprothint.py"`                     | MiniProtHint wrapper script for generating protein hints.           |
| `tools.bam2hints`                    | `"bam2hints"`                           | Converts RNA-Seq BAM alignments into genome annotation hint format. |



### Aligned RNA-Seq reads (BAM)
Already aligned short reads, for example from VARUS, go to `rnaseq_bam` as a
list of BAM files. They skip the HISAT2 step.

### pyVARUS output directories
[pyVARUS](https://github.com/Gaius-Augustus/pyVARUS) writes `stringtie.gtf`
and `hints.gff` next to `VARUS.bam`, made as the pipeline makes them from a
BAM. Pass the output directory instead of the BAM; this also works after
`varus run --drop-bam`. The directories skip mapping, `bam2hints` and
StringTie.

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
reads. TransDecoder flow: each directory adds its assembly to the
per-library assemblies. The Drusilla flow needs one assembly of all reads of
a read type: in `rnaseq` and `isoseq` mode pass exactly one directory and no
other reads of that type. In `mixed` mode it needs one `stringtie --mix`
assembly of both read types:

```sh
varus assemble genome.fa --short varus_short/VARUS.bam --long varus_long/VARUS.bam --outdir varus_mixed
```

```yaml
rnaseq_varus: /path/to/varus_short
isoseq_varus: /path/to/varus_long
mixed_varus:  /path/to/varus_mixed
```

Then pass no other short reads or Iso-Seq. `mixed_varus` is used only by the
Drusilla flow in mixed mode. `varus replay` rebuilds a BAM deleted by
`--drop-bam`.

### StringTie assemblies

A StringTie assembly made outside the pipeline goes to `stringtie`
(`--stringtie`): one GTF, a list, or a glob. It replaces mapping and
assembly, and it can be combined with reads and BAM files.

```yaml
stringtie: /path/to/stringtie.gtf
```

`stringtie` counts as short reads when the mode is inferred, so proteins and
an assembly give the mode `rnaseq`; a forced `isoseq` or `mixed` accepts it
as well. In the TransDecoder flow all given assemblies are merged with those
of the reads (`stringtie --merge`). The Drusilla flow takes exactly one
assembly and no other RNA-Seq or Iso-Seq input, and needs the `cov` and
`TPM` attributes that StringTie writes when it assembles reads; an assembly
without them (the output of `stringtie --merge`) stops the run. A missing
file is an error. No intron hints are made for a given assembly.

## Gene finder: Tiberius

To run Tiberius set `tiberius.run: true` and choose the model configuration
for your clade. The available configurations are the YAML files in
`tiberius/model_cfg/` of the Tiberius submodule (also listed in
[Tiberius/model_cfg](https://github.com/Gaius-Augustus/Tiberius/tree/main/model_cfg)).
The launcher resolves a name such as `diatoms` to that file.

| Parameter | Default | Description |
| --- | --- | --- |
| `tiberius.run` | `false` (the launcher sets `true` for the selected gene finder) | Run Tiberius and merge its predictions with the HC genes. |
| `tiberius.model_cfg` | none | Name of a Tiberius model configuration, or path to a configuration file. |
| `tiberius.model_dir` | none | Directory that holds the extracted weights directory (the archive name of `weights_url` in the model configuration without `.tar.gz`, e.g. `vertebrates_weights`), instead of the download at the start of the run. Tiberius then runs with `--model` and never downloads; the task fails if the directory is missing or empty. Without `model_dir` the archive of `weights_url` is downloaded once on the submitting host (label `download`) and staged into the tasks, so the GPU nodes need no internet. |
| `tiberius.result` | none | Existing Tiberius prediction (`--result`): the GFF3 or GTF of `tiberius.py`, or the `tiberius_ab_initio.gff3` of an earlier run. It is used instead of running Tiberius, in every mode; a missing file is an error. `tiberius.model_cfg` is then optional: it is not run, but the model of a vertebrate prediction selects the [Drusilla flow](drusilla_flow.md) and is the model of its hint rescue. |
| `tiberius.min_split_size` | `20000000` | Minimal size in bp of a genome chunk. |
| `tiberius.max_files` | `20` | Maximal number of genome chunks, which is the upper limit of parallel Tiberius tasks. |
| `tiberius.max_parallel` | unlimited | Cap of concurrently running Tiberius tasks, for example `1` on a single-GPU workstation. |
| `tiberius.batch_size` | automatic | Forwarded to `tiberius.py --batch_size`. Without it, Tiberius sizes the batch from the GPU memory; on GPUs above 80 GB (e.g. 96 GB RTX PRO 6000) that choice overflows int32 in TensorFlow, so the pipeline caps it (`bin/tiberius_batch_size.py`). |
| `tiberius.seq_len` | automatic | Forwarded to `tiberius.py --seq_len`. |

For the prediction the genome is split into smaller FASTA files, so that
Tiberius can run on several GPUs in parallel.

The pipeline downloads the weights of the model configuration
(`weights_url` in the YAML) once at the start of the run, on the submitting
host like the other downloads, and stages them into every Tiberius task, so
the GPU nodes need no internet. The same holds for the Tiberius model of the
hint rescue. To skip the download, for example when the submitting host has
no internet either, download and extract the archive once and set
`model_dir` to the directory that holds the extracted weights directory: the
archive name without `.tar.gz`, for example `vertebrates_weights` or `fungi`.
Tiberius runs with `--model` on these weights and never downloads; a task
fails if the directory is missing or empty.

## Gene finder: Vipsania

Set `vipsania.run: true` and `vipsania.model`. All parameters are described in
[vipsania.md](vipsania.md). Only one gene finder can run. `genefinder: tiberius|vipsania`
selects one explicitly; without it the block with `run: true` is used.

## Drusilla flow

For runs with transcripts and a vertebrate gene finder model, Drusilla ORFs
replace the TransDecoder HC genes and a LightGBM model filters the *ab initio*
predictions. The steps are described in
[drusilla_flow.md](drusilla_flow.md).

| Parameter | Default | Description |
| --- | --- | --- |
| `drusilla.run` | `auto` | `auto`: on for Tiberius models whose `target_species` is `Vertebrata` or `Mammalia` (`vertebrates`, `mammalia*`) and Vipsania `Vertebrata`, if the run has transcripts; `true`: always (an error without transcripts), e.g. for a custom vertebrate model of another clade; `false`: never. |
| `drusilla.model` | `vertebrates` | Released Drusilla model (`drusilla models list`). |
| `drusilla.weights`, `drusilla.config` | none | A local `.weights.h5` file and its architecture YAML instead of a released model. |
| `drusilla.cache_dir` | none | Drusilla model cache (`DRUSILLA_CACHE_DIR`), e.g. filled by `drusilla models download vertebrates`, for a submitting host without internet. Without it the archive of the model's `weights_url` in `drusilla/model_cfg/` is downloaded once on the submitting host (label `download`), checked against its `weights_sha256` and staged into the Drusilla task, so the GPU nodes need no internet. A model name that is not in `drusilla/model_cfg/` is downloaded by Drusilla in the task. |
| `drusilla.batch_size` | automatic | Drusilla batch size; without it Drusilla sizes the batch from the GPU memory. |
| `drusilla.shards` | `1` | Parallel Drusilla processes, each on a part of the transcripts. For CPU nodes, where one process uses only one to two cores (e.g. 24 with 48 CPUs). |
| `drusilla.min_coding_length` | `200` | Minimal CDS length of a Drusilla ORF. |
| `drusilla.fix_stop`, `drusilla.fix_start` | `true` | Stop codon fix of the ORFs with miniprot alignments; start codon fix with miniprothint start hints (needs `fix_stop`). |
| `drusilla.min_length`, `min_cov`, `min_tpm` | `300`, `3`, `1` | StringTie pre-filter: transcript length, coverage and TPM. |
| `drusilla.long_length`, `min_tpm_long` | `3000`, `0.5` | Relaxed TPM for transcripts of at least `long_length` nt. |
| `drusilla.lgb_model` | released `drusilla_lgb_3class_v1` | LightGBM model of the *ab initio* filter: the archive (URL or file), its unpacked directory, or a `.txt` text model with its `.json`. `null` switches the flow off. |
| `drusilla.lgb_model_sha256` | sha256 of the released archive | Checked after the download; set it to `null` for another model. |
| `drusilla.lgb_threshold` | `0.5` | Transcripts with P(partial) + P(correct) at or above it are candidates. |
| `drusilla.lgb_keep` | `correct` | Classes kept among the candidates, e.g. `correct,partial`. |
| `drusilla.rescue` | `true` | Hint rescue: loci of partial *ab initio* genes are predicted again by Tiberius with hints. |
| `drusilla.rescue_tiberius` | the one of the hint rescue image | Another `tiberius.py` with `--hints`. Without `--hints` the rescue is skipped with a warning. |
| `drusilla.rescue_model_cfg` | the Tiberius model of the run | Tiberius model of the rescue: a model configuration file, or the name of a model in `model_cfg/` of the hint rescue image. `vertebrates` for Vipsania runs and for `tiberius.result` without `tiberius.model_cfg`. The weights in `tiberius.model_dir` are used only for the model of the run; otherwise the weights of the rescue model are downloaded once on the submitting host (for a name, `weights_url` is read from `model_cfg/` of the Tiberius submodule; a name that is not there is downloaded by Tiberius in the task). |
| `drusilla.rescue_flank` | `25000` | Flank in bp around each rescue locus. |
| `drusilla.rescue_hint_weight` | `2.5` | `tiberius.py --hint_weight` of the rescue. |
| `drusilla.rescue_orf_filter` | `false` | ORF-agreement filter: skip rescue loci where a Drusilla ORF already has all introns of the best protein chain. Off by default, as benchmarked. |
| `drusilla.rescue_seq_len` | `99990` | `tiberius.py --seq_len` of the rescue (divisible by 18). |

## Mode

The pipeline infers its mode from the inputs, see the table in the
[README](../README.md#inputs-and-modes). Set `mode` to force one of
`abinitio`, `proteins`, `rnaseq`, `isoseq`, `mixed`. `tiberius` is accepted as
the historic name of `abinitio`; any other value is an error, as is a forced
mode whose inputs are missing (for example `rnaseq` without reads or a StringTie assembly).
