# Parameters

All inputs and settings of the pipeline go into one YAML file. A commented
template is [conf/parameters.yaml](../conf/parameters.yaml). Pass the file to
the launcher of your gene finder:

```bash
python tiberius.py --params_yaml parameters.yaml --nf_config slurm_generic
```

`--nf_config` takes a path or the name of a config in `conf/`, see
[hpc.md](hpc.md). The Tiberius launcher also accepts most parameters on the
command line; command line values override the file. The launcher writes the
merged parameters of a run to `<outdir>/params.yaml`.

## Input Data
### Genomic Sequences (Required Input)
Path to the genome FASTA file.

Example:
```bash
genome: /path/to/genome.fa
```

### Protein Sequences (Required Input)
Path to one or more FASTA files with protein sequences aligned to the target genome.
Multiple files are concatenated in input order.

Example:
```bash
proteins: /path/to/proteins.faa
```

Or a list:
```bash
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
```bash
odb12Partitions: [
  Metazoa,
  Fungi,
]
```

### RNA-Seq/Iso-Seq (Local)
The pipeline accepts local long and short read data. Add the absolute paths of your files to:
`rnaseq_single`, `rnaseq_paired`, `isoseq`.

Use absolute paths. Nextflow does **not** expand `~`, and relative paths are
resolved against the directory where Nextflow is launched, not the directory
of `params.yaml`.

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
| `threads`                  | `48`                                   | Global default number of CPU threads used by many processes.                                                 |
| `outdir`                   | `"results"`                            | Directory where all final results are written.                                                               |
| `scoring_matrix`           | `"conf/blosum62.csv"` | Amino acid substitution scoring matrix used by homology-based tools.                                         |

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
| `tools.diamond`                      | `"diamond"`                             | Diamond alignment tool for protein-to-genome searches.              |
| `tools.bedtools`                     | `"bedtools"`                            | Bedtools for genomic interval operations.                           |
| `tools.miniprot`                     | `"miniprot"`                            | MiniProt — protein-to-genome aligner.                               |
| `tools.miniprot_boundary_scorer`     | `"miniprot_boundary_scorer"`            | Additional MiniProt scoring utility for boundary refinement.        |
| `tools.miniprothint`                 | `"miniprothint.py"`                     | MiniProtHint wrapper script for generating protein hints.           |
| `tools.bam2hints`                    | `"bam2hints"`                           | Converts RNA-Seq BAM alignments into genome annotation hint format. |



### Aligned RNA-Seq reads (BAM)
Already aligned short reads, for example from VARUS, go to `rnaseq_bam` as a
list of BAM files. They skip the HISAT2 step.

## Gene finder: Tiberius

To run Tiberius set `tiberius.run: true` and choose the model configuration
for your clade. The available configurations are listed in
[Tiberius/model_cfg](https://github.com/Gaius-Augustus/Tiberius/tree/main/model_cfg)
and by `python tiberius.py --list_cfg`.

| Parameter | Default | Description |
| --- | --- | --- |
| `tiberius.run` | `false` (the Tiberius launcher sets `true`) | Run Tiberius and merge its predictions with the HC genes. |
| `tiberius.model_cfg` | none | Name of a Tiberius model configuration, or path to a configuration file. |
| `tiberius.result` | none | Existing Tiberius prediction (GTF/GFF3). It is used instead of running Tiberius. |
| `tiberius.min_split_size` | `20000000` | Minimal size in bp of a genome chunk. |
| `tiberius.max_files` | `20` | Maximal number of genome chunks, which is the upper limit of parallel Tiberius tasks. |
| `tiberius.max_parallel` | unlimited | Cap of concurrently running Tiberius tasks, for example `1` on a single-GPU workstation. |
| `tiberius.batch_size` | automatic | Forwarded to `tiberius.py --batch_size`. |
| `tiberius.seq_len` | automatic | Forwarded to `tiberius.py --seq_len`. |

For the prediction the genome is split into smaller FASTA files, so that
Tiberius can run on several GPUs in parallel.

## Gene finder: Vipsania

Set `vipsania.run: true` and `vipsania.model`. All parameters are described in
[vipsania.md](vipsania.md). Only one gene finder can run. `genefinder: tiberius|vipsania`
selects one explicitly; without it the block with `run: true` is used.

## Mode

The pipeline infers its mode from the inputs, see the table in the
[README](../README.md#inputs-and-modes). Set `mode` to force one of
`abinitio`, `proteins`, `rnaseq`, `isoseq`, `mixed`. `tiberius` is accepted as
the historic name of `abinitio`.
