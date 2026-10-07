nextflow.enable.dsl=2

// Post-processing of the final protein-coding annotation (docs/postprocessing.md):
//   1. sanity filter (internal stop codons, CDS length, structure),
//   2. UTRs from the StringTie assemblies (modes with transcripts),
//   3. the published <stem>.gff3, .gtf, _proteins.fa, _cds.fa,
//   4. completeness (compleasm, BUSCO) with qc.busco_lineage,
//   5. hint support, gene set statistics, OMArk, gffcompare,
//   6. ncRNA genes (ncrna.run): <stem>_with_ncRNA.gff3,
//   7. GO terms (fantasia.run): <stem>_go.gff3,
//   8. software versions and report.html.
// stem: <tool>_evidence, or <tool>_ab_initio in mode abinitio.

include { truthy } from '../lib_nf/functions.nf'
include { EMPTY_FILE as EMPTY_PLACEHOLDER } from '../modules/util.nf'
include { SANITY_FILTER; ADD_UTRS; FINALIZE_ANNOTATION; LONGEST_ISOFORM } from '../modules/finalize.nf'
include { DOWNLOAD_BUSCO_LINEAGE; FILTER_BUSCO_PROTEINS; COMPLEASM_GENOME; COMPLEASM_PROTEINS;
          BUSCO_GENOME; BUSCO_PROTEINS; COMPLETENESS_SUMMARY } from '../modules/completeness.nf'
include { GENE_SUPPORT; GENE_SET_STATISTICS; DOWNLOAD_NCBI_TAXONOMY; OMARK; GFFCOMPARE; SOFTWARE_VERSIONS;
          REPORT } from '../modules/qc.nf'
include { DOWNLOAD_RFAM; BARRNAP; TRNASCAN; CMSCAN; INFERNAL_TO_GFF3; FEELNC; FEELNC_TO_GFF3; MERGE_NCRNA } from '../modules/ncrna.nf'
include { STRINGTIE_MERGE as STRINGTIE_MERGE_ALL } from '../modules/assembly.nf'
include { SPLIT_GENOME as SPLIT_GENOME_NCRNA } from '../modules/genefinder.nf'
include { FANTASIA_ANNOTATE; FANTASIA_SUMMARY; FANTASIA_DECORATE as FANTASIA_DECORATE_CODING;
          FANTASIA_DECORATE as FANTASIA_DECORATE_NCRNA } from '../modules/fantasia.nf'

// BUSCO lineage with its OrthoDB version (eukaryota -> eukaryota_odb12), or null
def buscoLineage(value) {
    def lineage = value?.toString()?.trim()
    if( !lineage ) return null
    return lineage.contains('_odb') ? lineage : "${lineage}_odb12"
}

// A cache directory of downloaded databases: the given path, else
// ~/.cache/paludamentum/<name>. Created, so that tasks can stage it.
def cacheDir(value, String name) {
    def dir = file(value ? value.toString() : "${System.getenv('HOME')}/.cache/paludamentum/${name}")
    dir.mkdirs()
    return dir
}

workflow POSTPROCESS {

    take:
    gff         // channel: final protein-coding GFF3 (0 or 1 item)
    genome      // channel: genome FASTA (value)
    prefix      // value: name of the gene finder, 'tiberius' or 'vipsania'
    mode        // value: pipeline mode
    asm_short   // channel: StringTie assemblies of short reads (may be empty)
    asm_long    // channel: StringTie assemblies of Iso-Seq reads (may be empty)
    hints       // channel: hintsfile.gff (empty in mode abinitio)
    run_info    // value: map for the report (version, mode, genefinder, model, hc)

    main:
    def stem = mode == 'abinitio' ? "${prefix}_ab_initio" : "${prefix}_evidence"
    def transcripts = mode in ['rnaseq', 'isoseq', 'mixed']
    def empty = EMPTY_PLACEHOLDER()
    def versions = channel.empty()
    def qc_files = channel.empty()

    // ---- the final annotation
    def filtered = SANITY_FILTER(stem, gff, genome)
    qc_files = qc_files.mix(filtered.report)
    def pre_final = filtered.gff
    if( truthy(params.postprocess.utr) && transcripts ) {
        def utr = ADD_UTRS(stem, pre_final, asm_short.collect().ifEmpty([]), asm_long.collect().ifEmpty([]))
        pre_final = utr.gff
        qc_files = qc_files.mix(utr.report)
    }
    def fin = FINALIZE_ANNOTATION(stem, pre_final, genome)
    versions = versions.mix(fin.versions)
    def final_gff3 = fin.gff3.first()
    def final_gtf  = fin.gtf.first()
    def proteins   = fin.proteins.first()

    // ---- completeness
    def lineage = buscoLineage(params.qc.busco_lineage)
    def doCompleasm = lineage && truthy(params.qc.compleasm)
    def doBusco     = lineage && truthy(params.qc.busco)
    if( !lineage ) log.info "Completeness assessment off (no qc.busco_lineage)"
    if( doCompleasm || doBusco ) {
        def cache = cacheDir(params.qc.busco_download_path, 'busco')
        def lineages = cache.resolve('lineages')
        def odb = lineage.tokenize('_')[1]
        // what compleasm and BUSCO read; compleasm checks the .done markers
        def ready = ["${lineage}/dataset.cfg", "${lineage}.done", "eukaryota_${odb}.done",
                     'placement_files.done', 'file_versions.tsv.done'].every { f -> lineages.resolve(f).exists() }
        def busco_dir = ready ? channel.value(cache) :
            DOWNLOAD_BUSCO_LINEAGE(lineage, cache).done.map { _d -> cache }.first()
        def longest = LONGEST_ISOFORM(stem, final_gff3, genome)
        def busco_proteins = FILTER_BUSCO_PROTEINS(longest.proteins).proteins.first()
        def summaries = [busco_genome: empty, busco_proteins: empty, compleasm_genome: empty, compleasm_proteins: empty]
        if( doCompleasm ) {
            def cg = COMPLEASM_GENOME(genome, lineage, busco_dir)
            summaries.compleasm_genome   = cg.summary
            summaries.compleasm_proteins = COMPLEASM_PROTEINS(busco_proteins, lineage, busco_dir).summary
            versions = versions.mix(cg.versions)
        }
        if( doBusco ) {
            def bg = BUSCO_GENOME(genome, lineage, busco_dir)
            summaries.busco_genome   = bg.summary
            summaries.busco_proteins = BUSCO_PROTEINS(busco_proteins, lineage, busco_dir).summary
            versions = versions.mix(bg.versions)
        }
        def completeness = COMPLETENESS_SUMMARY(summaries.busco_genome, summaries.busco_proteins,
                                                summaries.compleasm_genome, summaries.compleasm_proteins)
        qc_files = qc_files.mix(completeness.tsv)
    }

    // ---- hint support and statistics
    def support = channel.value([])
    if( truthy(params.qc.gene_support) && mode != 'abinitio' ) {
        support = GENE_SUPPORT(final_gff3, hints).tsv
        qc_files = qc_files.mix(support)
    }
    if( truthy(params.qc.statistics) ) {
        def stats = GENE_SET_STATISTICS(final_gff3, support)
        qc_files = qc_files.mix(stats.text, stats.plots.flatten())
    }

    // ---- OMArk and gffcompare
    if( truthy(params.qc.omark) ) {
        if( !params.qc.omamer_db ) error "qc.omark = true needs the OMAmer database: qc.omamer_db (e.g. LUCA.h5 of https://omabrowser.org/All/)."
        def db = file(params.qc.omamer_db.toString())
        if( !db.exists() ) error "qc.omamer_db: no such file: ${db}"
        // the NCBI taxonomy as ete3's taxa.sqlite in the cache directory, built
        // once on the submitting host from taxdump.tar.gz (downloaded unless present)
        def taxaDir = cacheDir(params.qc.ete_taxa_path, 'ncbi_taxonomy')
        def taxa = taxaDir.resolve('taxa.sqlite').exists() ? channel.value(taxaDir) :
            DOWNLOAD_NCBI_TAXONOMY(taxaDir).done.map { _d -> taxaDir }.first()
        def omark = OMARK(proteins, final_gff3, db, taxa)
        qc_files = qc_files.mix(omark.summary)
        versions = versions.mix(omark.versions)
    }
    if( params.qc.reference_annotation ) {
        def reference = file(params.qc.reference_annotation.toString())
        if( !reference.exists() ) error "qc.reference_annotation: no such file: ${reference}"
        def cmp = GFFCOMPARE(final_gtf, reference)
        qc_files = qc_files.mix(cmp.stats)
        versions = versions.mix(cmp.versions)
    }

    // ---- ncRNA
    def ncrna_files = channel.empty()
    def ncrna_gff3 = null
    if( truthy(params.ncrna.run) ) {
        def rfamDir = cacheDir(params.ncrna.rfam_dir, 'rfam')
        def rfamReady = ['Rfam.cm', 'Rfam.clanin'].every { f -> rfamDir.resolve(f).exists() }
        if( params.ncrna.rfam_dir && !rfamReady )
            error "ncrna.rfam_dir: no Rfam.cm and Rfam.clanin in ${rfamDir}; leave ncrna.rfam_dir unset to download Rfam 15.1."
        def rfam = rfamReady ? channel.value(rfamDir) : DOWNLOAD_RFAM(rfamDir).done.map { _d -> rfamDir }.first()
        def rrna = BARRNAP(stem, genome)
        def trna = TRNASCAN(stem, genome)
        def chunks = SPLIT_GENOME_NCRNA(genome, 20000000, 20).chunks.flatten()
        def infernal = INFERNAL_TO_GFF3(stem, CMSCAN(chunks, rfam).tblout.collect())
        def lnc = channel.value([])
        if( transcripts && truthy(params.ncrna.lncrna) ) {
            def assemblies = asm_short.mix(asm_long).collect()
            def feelnc = FEELNC(stem, STRINGTIE_MERGE_ALL(assemblies).gtf, final_gtf, genome)
            lnc = FEELNC_TO_GFF3(stem, feelnc.gtf).gff
            ncrna_files = ncrna_files.mix(lnc)
        }
        ncrna_gff3 = MERGE_NCRNA(stem, final_gff3, rrna.gff, trna.gff, infernal.gff, lnc).gff3.first()
        ncrna_files = ncrna_files.mix(rrna.gff, trna.gff, infernal.gff)
        versions = versions.mix(rrna.versions, trna.versions)
    }

    // ---- GO terms
    def fantasia_files = channel.empty()
    def go_files = channel.empty()
    if( truthy(params.fantasia.run) ) {
        ['hf_cache_dir', 'lookup_dir'].each { key ->
            if( !params.fantasia[key] ) error "fantasia.run = true needs fantasia.${key} (see docs/postprocessing.md)."
            if( !file(params.fantasia[key].toString()).isDirectory() ) error "fantasia.${key}: not a directory: ${params.fantasia[key]}"
        }
        def results = FANTASIA_ANNOTATE(proteins, file(params.fantasia.hf_cache_dir.toString()),
                                        file(params.fantasia.lookup_dir.toString())).results.first()
        def summary = FANTASIA_SUMMARY(results)
        fantasia_files = fantasia_files.mix(summary.summary, summary.plot)
        go_files = go_files.mix(FANTASIA_DECORATE_CODING(stem, final_gff3, results).gff3)
        if( ncrna_gff3 != null )
            go_files = go_files.mix(FANTASIA_DECORATE_NCRNA("${stem}_with_ncRNA", ncrna_gff3, results).gff3)
    }

    // ---- software versions and report
    def pipeline_lines = ["Paludamentum\t${run_info.version}\t-"]
    def sv = SOFTWARE_VERSIONS(pipeline_lines, versions.collect().ifEmpty([]))
    qc_files = qc_files.mix(sv.tsv)
    if( truthy(params.qc.report) ) {
        def top = fin.gff3.mix(fin.gtf, fin.proteins, fin.cds, go_files)
        if( ncrna_gff3 != null ) top = top.mix(ncrna_gff3)
        def info = run_info + [stem: stem, busco_lineage: lineage]
        def params_yaml = file("${params.outdir}/params.yaml")   // written by the launcher before the run
        REPORT(top.collect(), qc_files.collect(), ncrna_files.collect().ifEmpty([]),
               fantasia_files.collect().ifEmpty([]), file("${params.outdir}/citations.md"),
               hints.collect().ifEmpty([]), params_yaml.exists() ? params_yaml : [], info)
    }

    emit:
    gff3 = fin.gff3
}
