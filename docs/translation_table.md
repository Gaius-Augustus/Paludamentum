# Genomes with another genetic code than the standard one

`translation_table` (`--translation_table`) is the NCBI translation table of
the nuclear genome. The default is 1, the standard code. Supported are the
nuclear codes that turn one or two stop codons of the standard code into an
amino acid:

| Table | Organisms (examples) | Reassigned | Stop codons |
| --- | --- | --- | --- |
| 1 | standard code | – | TAA, TAG, TGA |
| 6 | ciliates such as *Tetrahymena* and *Paramecium*, Dasycladaceae, *Hexamita* | TAA, TAG = Gln | TGA |
| 10 | *Euplotes* | TGA = Cys | TAA, TAG |
| 15 | *Blepharisma* | TAG = Gln | TAA, TGA |
| 29 | *Mesodinium* | TAA, TAG = Tyr | TGA |
| 30 | peritrich ciliates (*Vorticella*, *Zoothamnium*) | TAA, TAG = Glu | TGA |

Not supported: tables whose stop codon depends on its context (27, 28, 31),
table 12 and table 26 (a sense codon with another amino acid), and the codes
of organelles and bacteria.

```yaml
genome: genome.fa
proteins: [proteins.faa]
rnaseq_paired: "RNA/*_{1,2}.fastq.gz"
translation_table: 6
```

No gene finder model is given. The modes with transcript evidence (`rnaseq`,
`isoseq`, `mixed`) and mode `abinitio` are supported; mode `proteins` needs
Tiberius and is not.

## What changes with another table than 1

| Step | Table 1 | Another table |
| --- | --- | --- |
| Tiberius | runs when `tiberius.model_cfg` is given | does not run: its models are table-1 models. `tiberius.model_cfg` stops the run; an existing prediction (`tiberius.result`) is accepted |
| Vipsania | pretrained model (`vipsania.model`), finetuned on the genome | trained from scratch on the genome, then annotates it; see below |
| TD2 | `TD2.LongOrfs`, `TD2.Predict` | `-G <table>` for both. The ORF search reads the table; `TD2.Predict` 1.1.0 takes the option, and PSAURON scores with the standard code |
| Drusilla | `drusilla annotate` | `--translation-table <table>`: frame recoding, see below |
| miniprot | `miniprot --aln` | `-T <table>` |
| Proteins of gene sets and of the final annotation | `gffread -y -S` | `bin/translate_gff.py`, the same records translated with the table |
| Sanity filter | internal TAA, TAG, TGA remove a transcript | the stop codons of the table do; a CDS without stop codon is extended by a stop codon of the table |

Not adapted: BUSCO and compleasm in genome mode search the genome with the
standard code, so their genome completeness is too low for such a genome. The
completeness of the predicted proteins is not affected. The selection models
of `combine_genesets.py` were trained on genomes of the standard code.

## Vipsania

The genetic code is part of a Vipsania model, and the released models are
table-1 models. Vipsania trains without annotation, on the genome sequence
alone, so the pipeline trains a model for the genome (`TRAIN_VIPSANIA`, one
GPU task, `vipsania train --translation_table <table>`) and annotates the
genome with it (`RUN_VIPSANIA_TRAINED`). `vipsania.model` is not needed. If
it is given, it is not run; its name selects the clade settings of
`conf/hc_genes.yaml` only. `vipsania.result` skips both tasks.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `vipsania.train_config` | `configs/base_10M.json` | Model configuration; a relative path is one of the Vipsania image (`/opt/vipsania/configs`) |
| `vipsania.train_epochs` | `300` | Epochs |
| `vipsania.train_steps` | `800` | Steps per epoch |
| `vipsania.train_min_contig` | `20000` | Sequences of the genome below this length are left out of the training set (Vipsania skips sequences below 20 kb) |

`intermediate/vipsania/train_set.tsv` holds the number of sequences and bases
of the training set, `intermediate/vipsania/vipsania_train.log` the training
log. A genome without a sequence of `train_min_contig` bp stops the run.

## Drusilla

The Drusilla model was trained on genomes of the standard code and labels a
CDS with an in-frame TAA, TAG or TGA as non-coding. With
`--translation-table`, Drusilla scans every transcript three times, once per reading frame, with the reassigned stop codons of that frame
replaced by a standard codon of the same amino acid (table 6: TAA to CAA, TAG
to CAG). An ORF is kept when it lies in the frame its copy was recoded for;
ORFs of the three copies that overlap on the genome are resolved in favour of
the longest. The model is not retrained. Run time and memory of
`DRUSILLA_ANNOTATE` are three times those of a table-1 run.

The option is in the branch `translation-table` of Drusilla. The image
`gaiusaugustus/drusilla:0.3.0` that `conf/base.config` pins does not have it:
a run with another table than 1 needs an image built from that branch
(`process.withLabel: drusilla` in your config) until the next Drusilla
release is pinned.

Measured on *Tetrahymena thermophila* (table 6, GCF_000189635.1, 26,996
reference genes) with a StringTie assembly of three RNA-Seq libraries (35,375
transcripts), exact match of the CDS structure:

| ORF caller | ORFs | Exact | Precision | Sensitivity |
| --- | ---: | ---: | ---: | ---: |
| Drusilla, `--translation-table 6` | 24,513 | 12,912 | 0.527 | 0.478 |
| TD2, `-G 6` | 33,542 | 11,711 | 0.349 | 0.434 |
| Longest ORF from ATG to TGA | 26,443 | 11,721 | 0.443 | 0.434 |
| Drusilla without `--translation-table` (reference transcripts) | 1,401 | 31 | 0.022 | 0.001 |
| TD2 without `-G` | 3,705 | 20 | 0.005 | 0.001 |
