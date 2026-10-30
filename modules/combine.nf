nextflow.enable.dsl=2

// Integrate the four gene sets (Tiberius, Vipsania, revised TD2, Drusilla)
// with the DIAMOND tables, miniprot/miniprothint protein hints and RNA-Seq
// coverage/intron hints into one annotation (combined.gff3). The insect flag
// (passed from main.nf for Insecta/Lepidoptera gene finder models) adds
// --drop-lone-single-exon tiberius. All combine_genesets.py thresholds other
// than --min-locus use the script defaults.
//
// Inputs that may share a file name are staged apart (Nextflow refuses two
// input files of one name): the four DIAMOND tables are all diamond.tsv, and
// in mixed mode both intron hint files are named after the merged BAM.
process COMBINE_GENESETS {
  label 'container'
  publishDir "${params.outdir}/intermediate", pattern: "combined.{gff3,report.tsv}", mode: 'copy', overwrite: true

  input:
    path tib_gff
    path vip_gff
    path dru_gtf, stageAs: 'drusilla/*'        // drusilla_orfs.gtf, or the file of params.drusilla.result
    path td_gff
    // the intron hints per read type and the coverage track; an empty
    // placeholder (no BAM: pyVARUS directories, a given assembly without
    // reads) is left out of the command below
    path introns, stageAs: 'introns/?/*'
    path cov, stageAs: 'coverage/*'
    path dia_tib, stageAs: 'diamond/tiberius.tsv'
    path dia_vip, stageAs: 'diamond/vipsania.tsv'
    path dia_dru, stageAs: 'diamond/drusilla.tsv'
    path dia_td2, stageAs: 'diamond/td2.tsv'
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
  # RNA-Seq evidence only where there is some: without --introns and
  # --coverage the RNA rules of combine_genesets.py stay off, whereas an
  # empty hints file would count every intron as unsupported
  rna_args=""
  for f in ${introns}; do
      if [ -s "\$f" ]; then rna_args="\$rna_args \$f"; fi
  done
  if [ -n "\$rna_args" ]; then rna_args="--introns\$rna_args"; fi
  if [ -s ${cov} ]; then rna_args="\$rna_args --coverage ${cov}"; fi
  combine_genesets.py \\
      --tiberius ${tib_gff} \\
      --vipsania ${vip_gff} \\
      --drusilla ${dru_gtf} \\
      --td2      ${td_gff} \\
      --diamond-tiberius ${dia_tib} \\
      --diamond-vipsania ${dia_vip} \\
      --diamond-drusilla ${dia_dru} \\
      --diamond-td2      ${dia_td2} \\
      --diamond-format full \\
      --miniprot-train ${mp_train} \\
      --miniprot-rep   ${mp_rep} \\
      --protein-hints  ${hc_hints} \\
      --min-locus ${min_locus} ${drop} \\
      \$rna_args \\
      -o combined.gff3 \\
      --report combined.report.tsv
  """

  stub:
  """
  touch combined.gff3 combined.report.tsv
  """
}
