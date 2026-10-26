nextflow.enable.dsl=2

// Integrate the four gene sets (Tiberius, Vipsania, revised TD2, Drusilla)
// with the DIAMOND tables, miniprot/miniprothint protein hints and RNA-Seq
// coverage/intron hints into one annotation (combined.gff3). The insect flag
// (passed from main.nf for Insecta/Lepidoptera gene finder models) adds
// --drop-lone-single-exon tiberius. All combine_genesets.py thresholds other
// than --min-locus use the script defaults.
process COMBINE_GENESETS {
  label 'container'
  publishDir "${params.outdir}/intermediate", pattern: "combined.{gff3,report.tsv}", mode: 'copy', overwrite: true

  input:
    path tib_gff
    path vip_gff
    path dru_gtf
    path td_gff
    path introns
    path cov
    path dia_tib
    path dia_vip
    path dia_dru
    path dia_td2
    path mp_train
    path mp_rep
    path hc_hints
    val  insect
    val  min_locus

  output:
    path "combined.gff3",       emit: gff3
    path "combined.report.tsv", emit: report

  script:
  def drop = insect ? "--drop-lone-single-exon tiberius" : ""
  """
  combine_genesets.py \\
      --tiberius ${tib_gff} \\
      --vipsania ${vip_gff} \\
      --drusilla ${dru_gtf} \\
      --td2      ${td_gff} \\
      --introns  ${introns} \\
      --coverage ${cov} \\
      --diamond-tiberius ${dia_tib} \\
      --diamond-vipsania ${dia_vip} \\
      --diamond-drusilla ${dia_dru} \\
      --diamond-td2      ${dia_td2} \\
      --diamond-format full \\
      --miniprot-train ${mp_train} \\
      --miniprot-rep   ${mp_rep} \\
      --protein-hints  ${hc_hints} \\
      --min-locus ${min_locus} ${drop} \\
      -o combined.gff3 \\
      --report combined.report.tsv
  """

  stub:
  """
  touch combined.gff3 combined.report.tsv
  """
}
