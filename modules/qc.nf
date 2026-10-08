nextflow.enable.dsl=2

// Quality control of the final annotation (docs/postprocessing.md): hint
// support per transcript, gene set statistics, OMArk, gffcompare against a
// reference, the software versions and report.html.

process GENE_SUPPORT {
  label 'container'
  publishDir "${params.outdir}/qc", mode:'copy', overwrite: true

  input:
    path gff3
    path hints

  output:
    path "gene_support.tsv", emit: tsv

  script:
  """
  gene_support.py --gff3 ${gff3} --hints ${hints} --out gene_support.tsv
  """

  stub:
  """
  touch gene_support.tsv
  """
}

// gene_set_statistics.txt and its plots; the evidence plot with a support table
process GENE_SET_STATISTICS {
  label 'postprocess'
  publishDir "${params.outdir}/qc", mode:'copy', overwrite: true

  input:
    path gff3
    path support, stageAs: 'support/*'

  output:
    path "gene_set_statistics.txt", emit: text
    path "*.png", emit: plots

  script:
  def supportArg = support ? "--support ${support}" : ''
  """
  gene_set_statistics.py --gff3 ${gff3} --out-dir . ${supportArg}
  """

  stub:
  """
  touch gene_set_statistics.txt isoform_and_exon_structure.png transcript_lengths.png introns_per_gene.png
  ${support ? 'touch evidence_support.png' : ''}
  """
}

// The NCBI taxonomy for OMArk, once, on the submitting host (label
// 'download'): ete3's taxa.sqlite in the cache directory
// (params.qc.ete_taxa_path, else ~/.cache/paludamentum/ncbi_taxonomy), built
// from taxdump.tar.gz, which is downloaded there unless the directory holds
// it. ete3 (3.1.3) never reads a tarball from ~/.etetoolkit: without a
// database it fetches the md5 and the tarball from NCBI, so the OMArk task
// gets the built database (omark -e) and needs no internet.
process DOWNLOAD_NCBI_TAXONOMY {
  label 'download', 'omark'

  input:
    path cache

  output:
    path "downloaded.txt", emit: done

  script:
  """
  set -euo pipefail
  # one build at a time per cache directory (runs that start together); a run
  # that waited here finds taxa.sqlite complete
  exec 9> ${cache}/.paludamentum_download.lock
  flock 9
  if [ -s ${cache}/taxa.sqlite ]; then
      echo "taxa.sqlite is already in ${cache}"
      echo "NCBI taxonomy" > downloaded.txt
      exit 0
  fi
  if [ ! -f ${cache}/taxdump.tar.gz ]; then
      python - <<'END_DOWNLOAD'
import hashlib, os, shutil, urllib.request
url = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
part = "${cache}/taxdump.tar.gz.part"
with urllib.request.urlopen(url + ".md5") as response:
    md5 = response.read().decode().split()[0]
with urllib.request.urlopen(url) as response, open(part, "wb") as out:
    shutil.copyfileobj(response, out)
with open(part, "rb") as tarball:
    if hashlib.md5(tarball.read()).hexdigest() != md5:
        raise SystemExit("taxdump.tar.gz: md5 mismatch after the download")
os.replace(part, "${cache}/taxdump.tar.gz")
END_DOWNLOAD
  fi
  # ete3 builds the database (and <db>.traverse.pkl next to it) from the
  # tarball and writes its temporary tables into the working directory. Built
  # under a temporary name: a killed build must not leave a taxa.sqlite behind.
  rm -f ${cache}/taxa.sqlite.part ${cache}/taxa.sqlite.part.traverse.pkl
  python -c "from ete3 import NCBITaxa; NCBITaxa(dbfile='${cache}/taxa.sqlite.part', taxdump_file='${cache}/taxdump.tar.gz')"
  [ -s ${cache}/taxa.sqlite.part ] || { echo "ete3 built no ${cache}/taxa.sqlite" >&2; exit 1; }
  mv ${cache}/taxa.sqlite.part.traverse.pkl ${cache}/taxa.sqlite.traverse.pkl
  mv ${cache}/taxa.sqlite.part ${cache}/taxa.sqlite
  echo "NCBI taxonomy" > downloaded.txt
  """

  stub:
  """
  echo "NCBI taxonomy" > downloaded.txt
  """
}

// OMArk on the proteins of all isoforms, grouped per gene (isoforms.splice)
process OMARK {
  label 'omark'
  publishDir "${params.outdir}/qc", mode:'copy', overwrite: true, pattern: 'omark_summary.txt'

  input:
    path proteins
    path gff3
    path omamer_db
    path taxa          // directory with taxa.sqlite (DOWNLOAD_NCBI_TAXONOMY)

  output:
    path "omark_summary.txt", emit: summary
    path "versions.tsv", emit: versions

  script:
  """
  filter_proteins.py --in ${proteins} --out proteins.fa --strip-stop
  ${params.tools.omamer} search --db ${omamer_db} --query proteins.fa --nthreads ${task.cpus} --out proteome.omamer
  # one line per gene: its transcript IDs, separated by ';'
  awk -F'\\t' '\$3 == "mRNA" { id = ""; parent = "";
      n = split(\$9, a, ";");
      for( i = 1; i <= n; i++ ) { if( a[i] ~ /^ID=/ ) id = substr(a[i], 4); if( a[i] ~ /^Parent=/ ) parent = substr(a[i], 8) }
      if( !(parent in tx) ) order[++k] = parent;
      tx[parent] = (parent in tx) ? tx[parent] ";" id : id }
      END { for( i = 1; i <= k; i++ ) print tx[order[i]] }' ${gff3} > isoforms.splice
  # HOME is the task directory in the container, so that the caches of the
  # tools (matplotlib) stay there. The NCBI taxonomy is the staged taxa.sqlite
  # (-e); ete3 does not touch ~/.etetoolkit or the internet with it.
  export HOME=\$PWD
  ${params.tools.omark} -f proteome.omamer -d ${omamer_db} -i isoforms.splice -o omark -e ${taxa}/taxa.sqlite
  summary=\$(find omark -name '*_detailed_summary.txt' | head -n 1)
  [ -n "\$summary" ] || summary=\$(find omark -name '*.sum' | head -n 1)
  [ -n "\$summary" ] || { echo "OMArk wrote no summary" >&2; exit 1; }
  cp "\$summary" omark_summary.txt
  # omark has no --version (it prints its usage); the version of the package
  printf 'OMArk\\t%s\\t%s\\n' "\$(python -c 'import importlib.metadata as m; print(m.version("omark"))' 2>/dev/null || echo unknown)" "${task.container ?: 'none'}" > versions.tsv
  printf 'OMAmer\\t%s\\t%s\\n' "\$(${params.tools.omamer} --version 2>&1 | awk '{print \$NF}')" "${task.container ?: 'none'}" >> versions.tsv
  """

  stub:
  """
  touch omark_summary.txt versions.tsv
  """
}

// CDS-level comparison with a reference annotation, as BRAKER4 and the
// benchmarks of this pipeline evaluate (gffcompare --strict-match -e 3 -T)
process GFFCOMPARE {
  label 'postprocess'
  publishDir "${params.outdir}/qc", mode:'copy', overwrite: true, pattern: 'gffcompare.stats'

  input:
    path gtf
    path reference, stageAs: 'reference/*'

  output:
    path "gffcompare.stats", emit: stats
    path "versions.tsv", emit: versions

  script:
  """
  # decompress first (a gzipped GFF3 must not reach gffread), then decide the
  # format on the name without .gz: GTF is used as is, anything else is converted
  ref="${reference}"
  name="\${ref%.gz}"
  if [[ "\$ref" == *.gz ]]; then
      gunzip -c "\$ref" > "reference.\${name##*.}"
      ref="reference.\${name##*.}"
  fi
  case "\$name" in
      *.gtf) ;;
      *) ${params.tools.gffread} "\$ref" -T -o reference.gtf; ref=reference.gtf ;;
  esac
  awk -F'\\t' '\$3 == "CDS"' "\$ref" > reference.CDS.gtf
  awk -F'\\t' '\$3 == "CDS"' ${gtf} > prediction.CDS.gtf
  ${params.tools.gffcompare} --strict-match -e 3 -T -r reference.CDS.gtf -o gffcompare prediction.CDS.gtf
  printf 'GffCompare\\t%s\\t%s\\n' "\$(${params.tools.gffcompare} --version 2>&1 | awk '{print \$NF}')" "${task.container ?: 'none'}" > versions.tsv
  """

  stub:
  """
  touch gffcompare.stats versions.tsv
  """
}

// qc/software_versions.tsv: the version lines of the post-processing tasks,
// sorted and without duplicates, after the pipeline line
process SOFTWARE_VERSIONS {
  label 'container'
  publishDir "${params.outdir}/qc", mode:'copy', overwrite: true

  input:
    val pipeline_lines
    path versions, stageAs: 'versions/?/*'

  output:
    path "software_versions.tsv", emit: tsv

  script:
  def head = pipeline_lines.join('\n')
  """
  cat > software_versions.tsv <<'END_VERSIONS'
tool	version	image
${head}
END_VERSIONS
  cat versions/*/* 2>/dev/null | sort -u >> software_versions.tsv || true
  """

  stub:
  """
  touch software_versions.tsv
  """
}

// report.html from the published files, staged as they are published:
// top-level files (with hintsfile.gff and the launcher's params.yaml), qc/
// and ncrna/. completeness.png, drawn by the report, is published to qc/.
process REPORT {
  label 'postprocess'
  publishDir "${params.outdir}/", mode:'copy', overwrite: true, pattern: 'report.html'
  publishDir "${params.outdir}/qc", mode:'copy', overwrite: true, pattern: 'completeness.png'

  input:
    path top, stageAs: 'staged/*'
    path qc, stageAs: 'staged/qc/*'
    path ncrna, stageAs: 'staged/ncrna/*'
    path fantasia, stageAs: 'staged/qc/fantasia/*'
    path citations, stageAs: 'staged/citations.md'
    path hints, stageAs: 'staged/hintsfile.gff'       // [] in mode abinitio
    path params_yaml, stageAs: 'staged/params.yaml'   // [] when Nextflow was run without the launcher
    val run_info

  output:
    path "report.html", emit: html
    path "completeness.png", optional: true, emit: completeness

  script:
  def json = groovy.json.JsonOutput.toJson(run_info)
  """
  cat > run_info.json <<'END_RUN_INFO'
${json}
END_RUN_INFO
  paludamentum_report.py --dir staged --out report.html --run-info run_info.json --completeness-png completeness.png
  """

  stub:
  """
  touch report.html
  """
}
