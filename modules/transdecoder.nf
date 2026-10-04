nextflow.enable.dsl=2

include { orfFinder } from '../lib_nf/functions.nf'

//  TransDecoder (td1) or TD2 (td2): gtf->fasta, LongOrfs, Predict
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
  // params.transdecoder, else orf_finder of the HC table (conf/hc_genes.yaml)
  def td = orfFinder(params)
  if( td == 'td2' )
  """
  mkdir -p transdecoder

  # 1) GTF -> transcript FASTA
  ${params.tools.transdecoder_util_gtf2fa} ${gtf} ${genome} > transdecoder/transcripts.fasta

  # 2) + 3) TD2 writes the final files to the working directory; renamed to
  # the TransDecoder names for the HC steps
  cd transdecoder
  ${params.tools.td2_longorfs} -t transcripts.fasta -O transcripts.fasta.transdecoder_dir -@ ${task.cpus}
  ${params.tools.td2_predict}  -t transcripts.fasta -O transcripts.fasta.transdecoder_dir ${params.td2_predict_args ?: ''}
  for ext in pep cds gff3 bed; do mv transcripts.fasta.TD2.\$ext transcripts.fasta.transdecoder.\$ext; done

  # The FASTA headers follow TransDecoder except for minus-strand ORFs: TD2
  # writes <start>-<end> with start > end, TransDecoder <left>-<right>
  perl -i -pe 's/:(\\d+)-(\\d+)\\(-\\)\$/":" . (\$1 < \$2 ? "\$1-\$2" : "\$2-\$1") . "(-)"/e if /^>/' \\
      transcripts.fasta.transdecoder.pep transcripts.fasta.transdecoder.cds
  """
  else
  """
  mkdir -p transdecoder

  # 1) GTF -> transcript FASTA
  ${params.tools.transdecoder_util_gtf2fa} ${gtf} ${genome} > transdecoder/transcripts.fasta

  # 2) Long ORFs
  ${params.tools.transdecoder_longorfs} -O transdecoder -t transdecoder/transcripts.fasta

  # 3) Predict CDS/peptides
  ${params.tools.transdecoder_predict}  -O transdecoder -t transdecoder/transcripts.fasta
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
