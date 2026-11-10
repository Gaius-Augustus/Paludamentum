nextflow.enable.dsl=2

// TD2: StringTie GTF -> transcripts.fasta, LongOrfs, Predict. Produces
// transcript-coordinate ORFs that the TD2 partial-ORF rescue (SHORTEN /
// CDS_CLASSIFY_AND_REVISE) refines, then TD_REVISED_TO_GENOME maps the
// result to the genome.
process TD_ALL {
  label 'container'
  input:
    path gtf
    path genome

  output:
    path "transdecoder/transcripts.fasta",                    emit: cdna
    path "transdecoder/transcripts.fasta.transdecoder.pep",   emit: pep
    path "transdecoder/*.transdecoder_dir",                   emit: longdir

  script:
  """
  mkdir -p transdecoder

  # 1) GTF -> transcript FASTA
  ${params.tools.transdecoder_util_gtf2fa} ${gtf} ${genome} > transdecoder/transcripts.fasta

  # 2) + 3) TD2 writes the final files to the working directory; renamed to
  # the TransDecoder names for the downstream steps.
  cd transdecoder
  ${params.tools.td2_longorfs} -t transcripts.fasta -O transcripts.fasta.transdecoder_dir -@ ${task.cpus}
  ${params.tools.td2_predict}  -t transcripts.fasta -O transcripts.fasta.transdecoder_dir ${params.td2_predict_args ?: ''}
  for ext in pep cds gff3 bed; do mv transcripts.fasta.TD2.\$ext transcripts.fasta.transdecoder.\$ext; done

  # TD2 writes minus-strand ORF ranges as <start>-<end> with start > end;
  # TransDecoder convention is <left>-<right>. The downstream orf-to-genome
  # helpers expect the TransDecoder form.
  perl -i -pe 's/:(\\d+)-(\\d+)\\(-\\)\$/":" . (\$1 < \$2 ? "\$1-\$2" : "\$2-\$1") . "(-)"/e if /^>/' \\
      transcripts.fasta.transdecoder.pep transcripts.fasta.transdecoder.cds
  """

  stub:
  """
  mkdir -p transdecoder/transcripts.fasta.transdecoder_dir
  touch transdecoder/transcripts.fasta transdecoder/transcripts.fasta.transdecoder.pep
  """
}

process SHORTEN_INCOMPLETE_ORFS {
  label 'container'
  input:
    path pep
  output:
    path "shortened_candidates.pep", emit: pep_short
  script:
  """
  shorten_incomplete.py ${pep} -o shortened_candidates.pep
  """

  stub:
  """
  touch shortened_candidates.pep
  """
}

process CDS_CLASSIFY_AND_REVISE {
  label 'container'
  input:
    path diamond_normal,  stageAs: 'diamond_normal.tsv'
    path diamond_short,   stageAs: 'diamond_short.tsv'
    path transdecoder_pep
    path shortened_pep

  output:
    path "revised_candidates.pep", emit: revised_pep
    path "classifications.json",   emit: classes

  script:
  """
  revise_pep.py \
    --diamond_normal diamond_normal.tsv \
    --diamond_short  diamond_short.tsv  \
    --transdecoder_pep ${transdecoder_pep} \
    --shortened_pep    ${shortened_pep} \
    --revised_pep      revised_candidates.pep \
    --classifications_json classifications.json
  """

  stub:
  """
  touch revised_candidates.pep classifications.json
  """
}

// The revised TD2 peptides are in transcript coordinates (one entry per ORF
// of a StringTie transcript). cdna_alignment_orf_to_genome_orf.pl of
// TransDecoder maps them back onto the genome using the StringTie assembly;
// combine_genesets.py reads the result as the TD2 gene set. All isoforms are
// kept (combine_genesets does its own ranking).
process TD_REVISED_TO_GENOME {
  label 'container'
  publishDir "${params.outdir}/intermediate", pattern: "td2.gff3", mode: 'copy', overwrite: true
  input:
    path revised_pep
    path transcripts_fasta
    path stringtie_gtf
    path stringtie_gff3

  output:
    path "td2.gff3", emit: gff3

  script:
  """
  td_revised_to_genome.py \
      --revised-pep       ${revised_pep} \
      --stringtie-gtf     ${stringtie_gtf} \
      --stringtie-gff3    ${stringtie_gff3} \
      --transcripts-fasta ${transcripts_fasta} \
      --transdecoder-util ${params.tools.transdecoder_util_orf2genome} \
      -o td2.gff3
  """

  stub:
  """
  touch td2.gff3
  """
}
