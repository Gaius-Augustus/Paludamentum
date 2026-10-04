nextflow.enable.dsl=2

process MINIPROT_ALIGN {
  label 'container', 'bigmem'

  input: path genome; path proteins

  output:
    path "miniprot/miniprot.aln", emit: aln

  script: """
  mkdir -p miniprot
  ${params.tools.miniprot} -t ${task.cpus} --aln ${genome} ${proteins} > miniprot/miniprot.aln
  """

  stub:
  """
  mkdir -p miniprot
  touch miniprot/miniprot.aln
  """
}

process MINIPROT_BOUNDARY_SCORE {
  label 'container'

  input:
    path aln
    path score_matrix

  output:
    path "miniprot/miniprot_parsed.gff", emit: gff

  script:
  """
  mkdir -p miniprot
  ${params.tools.miniprot_boundary_scorer} \
    -s ${score_matrix} \
    -o miniprot/miniprot_parsed.gff \
    < ${aln}
  """

  stub:
  """
  mkdir -p miniprot
  touch miniprot/miniprot_parsed.gff
  """
}

process MINIPROTHINT_CONVERT {
  label 'container'

  input: path gff

  output:
    path "miniprot/miniprot.gtf", emit: gtf
    path "miniprot/miniprot_trainingGenes.gff", emit: traingff
    path "miniprot/hc.gff", emit: hc_hints

  script: """
  mkdir -p miniprot
  ${params.tools.miniprothint} ${gff} --workdir miniprot --ignoreCoverage --topNperSeed 10 --minScoreFraction 0.5
  """

  stub:
  """
  mkdir -p miniprot
  touch miniprot/miniprot.gtf miniprot/miniprot_trainingGenes.gff miniprot/hc.gff
  """
}

process ALN2HINTS {
  label 'container'

  input: path gtf
  output: path "hints_protein.gff", emit: hints

  script: """
  aln2hints.pl --in=${gtf} --out=prot_hintsfile.aln2hints.temp.gff --prg=miniprot --priority=4
  cp prot_hintsfile.aln2hints.temp.gff hints_protein.gff
  """

  stub:
  """
  touch hints_protein.gff
  """
}

process PREPROCESS_PROTEINDB {
  label 'container', 'bigmem'

  input:
    path proteinDB
    path tiberius_prot

  output: path "protein_preprocessed.fa"

  script: """
    # Count protein sequences (FASTA headers start with '>'); grep -c prints 0
    # itself for an empty database, its exit code 1 must not stop the task
    N_PROT=\$(grep -c '^>' ${proteinDB} || true)
    # The species selection below needs OrthoDB-style ids (<taxid>_<n>:<hex>);
    # other databases (UniProt, NCBI) carry no species in the id.
    N_ODB=\$(grep -c -E '^>[0-9]+_[0-9]+:' ${proteinDB} || true)

    echo "[PREPROCESS_PROTEINDB] Number of proteins in input: \$N_PROT (OrthoDB-style ids: \$N_ODB)" >&2

    if [[ "\$N_PROT" -le 1000000 || "\$N_ODB" -lt "\$N_PROT" ]]; then
        if [[ "\$N_PROT" -le 1000000 ]]; then
            echo "[PREPROCESS_PROTEINDB] <= 1,000,000 proteins – using full DB." >&2
        else
            echo "[PREPROCESS_PROTEINDB] > 1,000,000 proteins, but the ids are not OrthoDB-style, so no species can be selected – using full DB." >&2
        fi
        # Sanitize headers: keep only the first whitespace-delimited token after '>'.
        # An embedded tab (e.g. OrthoDB-style headers) survives into DIAMOND's
        # sseqid as the literal two characters '\\t', breaking downstream
        # id matching in hc_module.getting_hc_supported_by_proteins.
        awk '/^>/ { split(substr(\$0,2), a, /[ \\t]/); print ">" a[1]; next } { print }' \
            ${proteinDB} > protein_preprocessed.fa
    else
        echo "[PREPROCESS_PROTEINDB] > 1,000,000 proteins – running DIAMOND soft filter." >&2

        # Sanitize headers on the way into the database (same as above, streamed,
        # no copy of the FASTA). OrthoDB headers carry a tab, which DIAMOND
        # reports once per sequence ("Tabulator character in sequence title"):
        # gigabytes of log for millions of proteins. The per-sequence warnings
        # are dropped from stderr; everything else DIAMOND reports is kept.
        awk '/^>/ { split(substr(\$0,2), a, /[ \\t]/); print ">" a[1]; next } { print }' \
            ${proteinDB} \
          | ${params.tools.diamond} makedb --db prot_db --threads ${task.cpus} \
              2> >(grep -v -e 'Tabulator character in sequence title' >&2)

        # gffread emits '.' for internal stop codons in malformed CDS predictions;
        # DIAMOND rejects them. Replace with '*' (canonical stop) so the query
        # is valid and the species-ranking pass can still run.
        awk '/^>/{print; next} {gsub(/\\./, "*"); print}' \
            ${tiberius_prot} > tiberius_proteins.clean.fa

        ${params.tools.diamond} blastp \
          --query tiberius_proteins.clean.fa \
          --db prot_db \
          --out diamond_hits.tsv \
          --outfmt 6 qseqid sseqid pident length evalue bitscore qlen slen \
          --evalue 1e-5 \
          --max-target-seqs 200 \
          --very-sensitive \
          --threads ${task.cpus}

        rank_species_from_diamond.py diamond_hits.tsv 13 > species_rank.tsv

        if [[ ! -s top_species.txt ]]; then
            echo "[PREPROCESS_PROTEINDB] No species could be ranked (no DIAMOND hits): using the full DB." >&2
            awk '/^>/ { split(substr(\$0,2), a, /[ \\t]/); print ">" a[1]; next } { print }' \
                ${proteinDB} > protein_preprocessed.fa
            exit 0
        fi

        awk '
        BEGIN {
          while ((getline < "top_species.txt") > 0) {
            wanted[\$1] = 1
          }
        }
        /^>/ {
          hdr = substr(\$0, 2)
          split(hdr, a, /[ \t]/)
          id = a[1]                 # e.g. 101020_0:000003
          species = id
          sub(/_.*/, "", species)   # species = 101020
          keep = (species in wanted)
        }
        keep && /^>/ { print ">" id; next }
        keep { print }
      ' ${proteinDB} > protein_preprocessed.fa
    fi
  """

  stub:
  """
  touch protein_preprocessed.fa
  """
}
