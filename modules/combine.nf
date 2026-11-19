nextflow.enable.dsl=2

// Integrate the four gene sets (Tiberius, Vipsania, revised TD2, Drusilla)
// with the DIAMOND tables, miniprot/miniprothint protein hints, RNA-Seq
// coverage/intron hints and, when there are any, full-length transcripts
// (--chains: the StringTie -L assembly of long reads, or params.chains) into
// one annotation (combined.gff3). The insect flag (passed from main.nf for
// Insecta/Lepidoptera gene finder models) adds --drop-lone-single-exon
// tiberius. Selection: the learned classifier (params.combine.classifier,
// a JSON model, with the options of clf_opts: --unknown-clade when the clade
// has neither classifier_clade nor standard_thresholds in the HC table, --t-locus / --t-iso from
// params.combine.t_locus / t_iso) or, with the empty placeholder instead of
// a model, the tier rules with --min-locus. The genome goes in for the
// assembly L50: in frag_opts (--fragmented-l50 N --t-locus-fragmented T, from
// params.combine) a fragmented assembly gets the lower locus threshold. The
// same string carries the protein rescue (--rescue-protein-cov C
// --rescue-protein-introns F, from params.combine.rescue_protein_*): a locus
// below the locus threshold is kept if a multi-exon candidate has a
// full-length protein hit and protein hints for its introns. The
// assembly statistics and the thresholds used are written to
// combined.settings.tsv. All other combine_genesets.py thresholds use the
// script defaults.
//
// Inputs that may share a file name are staged apart (Nextflow refuses two
// input files of one name): the four DIAMOND tables are all diamond.tsv, and
// in mixed mode both intron hint files are named after the merged BAM.
//
// Mode proteins (HC_COMBINED with transcripts = false) passes the empty
// placeholder for the Drusilla and TD2 gene sets, their DIAMOND tables, the
// intron hints, the coverage track and the chains: combine_genesets.py reads
// an empty gene set and an empty DIAMOND table as no structures and no hits,
// and the RNA rules stay off as without reads. The classifier is then the
// protein-only model (params.combine.classifier_proteins). The ab initio
// combination (mode abinitio with both gene finders) also passes the
// placeholder for the miniprot models, the miniprothint hints and the DIAMOND
// tables of the gene finders, and selects with params.combine.classifier_abinitio.
process COMBINE_GENESETS {
  label 'container'
  publishDir "${params.outdir}/intermediate", pattern: "combined.{gff3,report.tsv,settings.tsv}", mode: 'copy', overwrite: true

  input:
    path tib_gff
    path vip_gff
    path dru_gtf, stageAs: 'drusilla/*'        // drusilla_orfs.gtf, or the file of params.drusilla.result
    path td_gff
    // the intron hints per read type and the coverage track (bam2hints and
    // genomecov of the BAM, or both derived from a given assembly by
    // ASSEMBLY_HINTS); an empty placeholder (no BAM: pyVARUS directories)
    // is left out of the command below
    path introns, stageAs: 'introns/?/*'
    path cov, stageAs: 'coverage/*'
    // full-length transcripts (GTF/GFF3 with exon lines) with the placeholder,
    // or the placeholder alone: then combine_genesets.py runs without --chains
    path chains, stageAs: 'chains/?/*'
    path dia_tib, stageAs: 'diamond/tiberius.tsv'
    path dia_vip, stageAs: 'diamond/vipsania.tsv'
    path dia_dru, stageAs: 'diamond/drusilla.tsv'
    path dia_td2, stageAs: 'diamond/td2.tsv'
    // the miniprot models and miniprothint hints, staged apart because all
    // three are the empty placeholder in the ab initio combination (no
    // protein evidence) and are then left out of the command below
    path mp_train, stageAs: 'miniprot_train/*'
    path mp_rep, stageAs: 'miniprot_rep/*'
    path hc_hints, stageAs: 'miniprot_hints/*'
    // the classifier model (JSON), or the empty placeholder for the tier rules
    path classifier, stageAs: 'classifier/*'
    // the genome, read only for its L50 (fragmented assemblies, frag_opts)
    path genome, stageAs: 'genome/*'
    val  insect
    val  min_locus
    // further classifier options as one string: --unknown-clade, --t-locus X, --t-iso Y (may be empty)
    val  clf_opts
    // --fragmented-l50 N --t-locus-fragmented T (fragmented-assembly threshold) and
    // --rescue-protein-cov C --rescue-protein-introns F (protein rescue), each only where it applies; may be empty
    val  frag_opts

  output:
    path "combined.gff3",         emit: gff3
    path "combined.report.tsv",   emit: report
    path "combined.settings.tsv", emit: settings

  script:
  def drop = insect ? "--drop-lone-single-exon tiberius" : ""
  """
  # learned selection with the shipped model; the empty placeholder means tier rules
  clf_args=""
  if [ -s ${classifier} ]; then clf_args="--classifier ${classifier} ${clf_opts} ${frag_opts}"; fi
  # RNA-Seq evidence only where there is some: without --introns and
  # --coverage the RNA rules of combine_genesets.py stay off, whereas an
  # empty hints file would count every intron as unsupported
  rna_args=""
  for f in ${introns}; do
      if [ -s "\$f" ]; then rna_args="\$rna_args \$f"; fi
  done
  if [ -n "\$rna_args" ]; then rna_args="--introns\$rna_args"; fi
  if [ -s ${cov} ]; then rna_args="\$rna_args --coverage ${cov}"; fi
  chain_args=""
  for f in ${chains}; do
      if [ -s "\$f" ]; then chain_args="\$chain_args \$f"; fi
  done
  if [ -n "\$chain_args" ]; then chain_args="--chains\$chain_args"; fi
  # protein evidence only where there is some (not in the ab initio combination)
  prot_args=""
  if [ -s ${mp_train} ]; then prot_args="\$prot_args --miniprot-train ${mp_train}"; fi
  if [ -s ${mp_rep} ];   then prot_args="\$prot_args --miniprot-rep ${mp_rep}"; fi
  if [ -s ${hc_hints} ]; then prot_args="\$prot_args --protein-hints ${hc_hints}"; fi
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
      --min-locus ${min_locus} ${drop} \\
      \$prot_args \$clf_args \$rna_args \$chain_args \\
      --genome ${genome} \\
      -o combined.gff3 \\
      --report combined.report.tsv \\
      --settings combined.settings.tsv
  """

  stub:
  """
  touch combined.gff3 combined.report.tsv combined.settings.tsv
  """
}
