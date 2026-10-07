nextflow.enable.dsl=2

include { truthy } from '../lib_nf/functions.nf'

// Non-coding RNA annotation (params.ncrna.run, docs/postprocessing.md): rRNA
// (pybarrnap), tRNA (tRNAscan-SE), Rfam families (Infernal cmscan, Rfam 15.1),
// lncRNA (FEELnc, modes with transcripts), merged with the protein-coding
// annotation into <stem>_with_ncRNA.gff3. The commands follow BRAKER4
// rules/ncrna/*.smk (commit 3535ed3). IDs carry the stem as prefix.

// Rfam 15.1 (Rfam.cm, Rfam.clanin, cmpress index), once, on the submitting
// host, written straight into the cache directory (staged as `rfam`). The
// version matches bin/rfam_family_types.tsv.
process DOWNLOAD_RFAM {
  label 'download', 'postprocess'

  input:
    path rfam

  output:
    path "downloaded.txt", emit: done

  script:
  """
  set -euo pipefail
  base=https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1
  curl -fsSL -o ${rfam}/Rfam.cm.gz \$base/Rfam.cm.gz
  curl -fsSL -o ${rfam}/Rfam.clanin \$base/Rfam.clanin
  ( cd ${rfam} && printf '%s\\n' \\
      'e2e636e9ce138dae506769aa74d067ce35466e94002818cfe125ea2a3cab8168  Rfam.cm.gz' \\
      '8055b3aec9be36004663b7b7643ea128b906e6c9656266189afd30fc0da34a1d  Rfam.clanin' \\
      | sha256sum -c - )
  gunzip -f ${rfam}/Rfam.cm.gz
  rm -f ${rfam}/Rfam.cm.i1?
  ${params.tools.cmpress} ${rfam}/Rfam.cm
  echo "Rfam 15.1" > downloaded.txt
  """

  stub:
  """
  echo "Rfam 15.1" > downloaded.txt
  """
}

process BARRNAP {
  label 'postprocess'
  publishDir "${params.outdir}/ncrna", mode:'copy', overwrite: true, pattern: 'rRNA.gff3'

  input:
    val stem
    path genome

  output:
    path "rRNA.gff3", emit: gff
    path "versions.tsv", emit: versions

  script:
  """
  ${params.tools.barrnap} --kingdom euk --threads ${task.cpus} ${genome} > rRNA.gff3.tmp
  if [ -s rRNA.gff3.tmp ] && grep -qv '^#' rRNA.gff3.tmp; then
      awk -F'\\t' -v OFS='\\t' -v p="${stem}" '
          BEGIN { n = 1 }
          /^#/ { print; next }
          {
              oldname = ""
              if( match(\$9, /Name=[^;]+/) ) oldname = substr(\$9, RSTART + 5, RLENGTH - 5)
              newname = p "-rRNA_" n "_" oldname
              gsub(/Name=[^;]+/, "Name=" newname, \$9)
              if( \$9 ~ /ID=/ ) gsub(/ID=[^;]+/, "ID=" newname, \$9)
              else \$9 = "ID=" newname ";" \$9
              n++
              print
          }' rRNA.gff3.tmp > rRNA.gff3
  else
      echo "##gff-version 3" > rRNA.gff3
  fi
  rm -f rRNA.gff3.tmp
  printf 'pybarrnap\\t%s\\t%s\\n' "\$(${params.tools.barrnap} --version 2>&1 | awk '{print \$NF}')" "${task.container ?: 'none'}" > versions.tsv
  """

  stub:
  """
  echo "##gff-version 3" > rRNA.gff3
  touch versions.tsv
  """
}

process TRNASCAN {
  label 'postprocess'
  publishDir "${params.outdir}/ncrna", mode:'copy', overwrite: true, pattern: 'tRNAs.{gff3,txt}'

  input:
    val stem
    path genome

  output:
    path "tRNAs.gff3", emit: gff
    path "tRNAs.txt", emit: txt
    path "versions.tsv", emit: versions

  script:
  def highconf = truthy(params.ncrna.trnascan_high_confidence)
  """
  LC_ALL=C ${params.tools.trnascan} -E ${highconf ? '-H --detail -f tRNAs.ss' : ''} --thread ${task.cpus} -q --forceow \\
      -o tRNAs.txt --gff tRNAs.gff3.tmp ${genome}
  touch tRNAs.txt tRNAs.gff3.tmp
  if ${highconf} && [ -s tRNAs.gff3.tmp ] && grep -qv '^#' tRNAs.gff3.tmp; then
      LC_ALL=C EukHighConfidenceFilter --result tRNAs.txt --ss tRNAs.ss --output . \\
          --prefix tRNAs.highconf --remove
      # keep the GFF lines of the high-confidence tRNAs (<seq>.trna<N>)
      awk -F'\\t' -v OFS='\\t' '
          FNR == NR {
              if( FNR > 3 ) { seq = \$1; gsub(/[[:space:]]/, "", seq); n = \$2; gsub(/[[:space:]]/, "", n); keep[seq ".trna" n] = 1 }
              next
          }
          /^#/ { print; next }
          {
              key = ""
              if( \$3 == "exon" ) { if( match(\$9, /Parent=[^;]+/) ) key = substr(\$9, RSTART + 7, RLENGTH - 7) }
              else if( match(\$9, /ID=[^;]+/) ) key = substr(\$9, RSTART + 3, RLENGTH - 3)
              if( key in keep ) print
          }' tRNAs.highconf.out tRNAs.gff3.tmp > tRNAs.gff3.hc
      mv tRNAs.gff3.hc tRNAs.gff3.tmp
  fi
  if [ -s tRNAs.gff3.tmp ] && grep -qv '^#' tRNAs.gff3.tmp; then
      prefix_trna_ids.sh "${stem}" tRNAs.gff3.tmp > tRNAs.gff3
  else
      echo "##gff-version 3" > tRNAs.gff3
  fi
  rm -f tRNAs.gff3.tmp
  printf 'tRNAscan-SE\\t%s\\t%s\\n' "\$(LC_ALL=C ${params.tools.trnascan} --help 2>&1 | sed -n 's/.*tRNAscan-SE \\([0-9.][0-9.]*\\).*/\\1/p' | head -n 1)" "${task.container ?: 'none'}" > versions.tsv
  """

  stub:
  """
  echo "##gff-version 3" > tRNAs.gff3
  touch tRNAs.txt versions.tsv
  """
}

// One task per genome chunk (SPLIT_GENOME); chunks hold whole sequences, so
// the --cut_ga hits and the clan overlap filter are those of the whole genome.
process CMSCAN {
  tag "${chunk.baseName}"
  label 'postprocess', 'bigmem'

  input:
    path chunk
    path rfam

  output:
    path "${chunk.baseName}.tblout", emit: tblout

  script:
  """
  # a user's rfam_dir without the cmpress index: press a copy in the task
  cm=${rfam}/Rfam.cm
  if [ ! -f "\$cm.i1m" ]; then
      cp ${rfam}/Rfam.cm Rfam.cm
      ${params.tools.cmpress} Rfam.cm
      cm=Rfam.cm
  fi
  ${params.tools.cmscan} --cut_ga --rfam --nohmmonly --clanin ${rfam}/Rfam.clanin --fmt 2 \\
      --cpu ${task.cpus} --tblout ${chunk.baseName}.tblout "\$cm" ${chunk} > /dev/null
  """

  stub:
  """
  touch ${chunk.baseName}.tblout
  """
}

process INFERNAL_TO_GFF3 {
  label 'container'
  publishDir "${params.outdir}/ncrna", mode:'copy', overwrite: true

  input:
    val stem
    path tblouts, stageAs: 'tblout/?/*'

  output:
    path "infernal.tblout", emit: tblout
    path "ncRNAs_infernal.gff3", emit: gff

  script:
  """
  cat ${tblouts} > infernal.tblout
  infernal_to_gff3.py -i infernal.tblout -o ncRNAs_infernal.gff3 -p ${stem} \\
      --family-types \$(dirname \$(command -v infernal_to_gff3.py))/rfam_family_types.tsv
  """

  stub:
  """
  touch infernal.tblout
  echo "##gff-version 3" > ncRNAs_infernal.gff3
  """
}

// lncRNA candidates of the StringTie assemblies (merged) that do not overlap
// the final annotation and have no coding potential (FEELnc)
process FEELNC {
  label 'feelnc'
  publishDir "${params.outdir}/ncrna", mode:'copy', overwrite: true, pattern: '{lncRNAs.gff3,feelnc_classifier.txt}'

  input:
    val stem
    path assembly
    path gtf
    path genome

  output:
    path "lncRNAs.gff3", emit: gff
    path "feelnc_classifier.txt", emit: classifier

  script:
  """
  FEELnc_filter.pl -i ${assembly} -a ${gtf} --monoex=-1 --size=200 -p ${task.cpus} > candidate_lncrna.gtf || true
  n=\$(awk -F'\\t' '\$3 == "transcript"' candidate_lncrna.gtf | wc -l)
  echo "##gff-version 3" > lncRNAs.gff3
  echo "# No lncRNA candidates" > feelnc_classifier.txt
  if [ "\$n" -gt 0 ]; then
      FEELnc_codpot.pl -i candidate_lncrna.gtf -a ${gtf} -g ${genome} --mode=shuffle --outdir=codpot_out -p ${task.cpus} || true
      lnc=\$(find codpot_out -name '*_lncRNA.gtf' 2>/dev/null | head -n 1)
      if [ -n "\$lnc" ] && [ -s "\$lnc" ]; then
          FEELnc_classifier.pl -i "\$lnc" -a ${gtf} > feelnc_classifier.txt || true
          # lnc_RNA + exons, IDs <stem>-lncRNA_<n>
          awk -F'\\t' -v OFS='\\t' -v p="${stem}" '
              function tx_id(attrs) {
                  if( match(attrs, /transcript_id "[^"]+"/) ) return substr(attrs, RSTART + 15, RLENGTH - 16)
                  return ""
              }
              function new_id(tid) {
                  if( !(tid in id) ) { id[tid] = p "-lncRNA_" n; n++ }
                  return id[tid]
              }
              BEGIN { print "##gff-version 3"; n = 1 }
              /^#/ { next }
              \$3 == "transcript" { lid = new_id(tx_id(\$9)); \$3 = "lnc_RNA"; \$9 = "ID=" lid ";Name=" lid ";biotype=lncRNA"; print }
              \$3 == "exon" { lid = new_id(tx_id(\$9)); k[lid]++; \$9 = "ID=" lid ".exon" k[lid] ";Parent=" lid; print }
          ' "\$lnc" > lncRNAs.gff3
      fi
  fi
  """

  stub:
  """
  echo "##gff-version 3" > lncRNAs.gff3
  touch feelnc_classifier.txt
  """
}

// <stem>_with_ncRNA.gff3: the final annotation and the ncRNA genes (priority
// rRNA > tRNA > Rfam > lncRNA), validated against the GFF3 contract
process MERGE_NCRNA {
  label 'container'
  publishDir "${params.outdir}/", mode:'copy', overwrite: true

  input:
    val stem
    path coding
    path rrna
    path trna
    path infernal
    path lncrna, stageAs: 'lnc/*'

  output:
    path "${stem}_with_ncRNA.gff3", emit: gff3

  script:
  """
  merge_ncrna_gff3.py --coding ${coding} --ncrna ${rrna} ${trna} ${infernal} ${lncrna} -o ${stem}_with_ncRNA.gff3
  validate_gff3.sh ${stem}_with_ncRNA.gff3
  """

  stub:
  """
  touch ${stem}_with_ncRNA.gff3
  """
}
