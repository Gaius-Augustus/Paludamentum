// References of the tools that the pipeline runs, and the citation file of a run.
//
// The DOIs were checked against the metadata registered at doi.org. The same
// list is the References section of the README; tests/test_launcher.py checks
// that the two agree.

// key: [tool, use (what the pipeline does with it), ref, doi or url]
def references() {
    return [
        paludamentum: [tool: 'Paludamentum', use: 'the pipeline',
            ref: 'Gabriel L, Hoff KJ. Paludamentum: evidence integration pipeline for the Gaius-Augustus gene finders. 2026.',
            url: 'https://github.com/Gaius-Augustus/Paludamentum'],
        nextflow: [tool: 'Nextflow', use: 'workflow engine',
            ref: 'Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nature Biotechnology. 2017;35(4):316-319.',
            doi: '10.1038/nbt.3820'],
        tiberius: [tool: 'Tiberius', use: 'gene prediction',
            ref: 'Gabriel L, Becker F, Hoff KJ, Stanke M. Tiberius: end-to-end deep learning with an HMM for gene prediction. Bioinformatics. 2024;40(12):btae685.',
            doi: '10.1093/bioinformatics/btae685'],
        tiberius_clades: [tool: 'Tiberius', use: 'models of the non-mammalian clades',
            ref: 'Gabriel L, Brůna T, Kaur A, Krishnan A, Ortmann F, Salamov A, Talbot S, Becker F, Krieg R, Wheat CW, Grigoriev IV, Stanke M, Hoff KJ. Accurate *ab initio* gene prediction in eukaryotes with Tiberius in multiple clades. bioRxiv. 2026.',
            doi: '10.64898/2026.04.24.720536'],
        vipsania: [tool: 'Vipsania', use: 'gene prediction',
            ref: 'Krieg R, Becker F, Saenko S, Diehl J, Stanke M. Vipsania: Unsupervised Deep Gene Finding. bioRxiv. 2026.',
            doi: '10.64898/2026.08.26.747235'],
        orthodb: [tool: 'OrthoDB v12', use: 'protein database',
            ref: 'Tegenfeldt F, Kuznetsov D, Manni M, Berkeley M, Zdobnov EM, Kriventseva EV. OrthoDB and BUSCO update: annotation of orthologs with wider sampling of genomes. Nucleic Acids Research. 2025;53(D1):D516-D522.',
            doi: '10.1093/nar/gkae987'],
        miniprot: [tool: 'miniprot', use: 'protein to genome alignment',
            ref: 'Li H. Protein-to-genome alignment with miniprot. Bioinformatics. 2023;39(1):btad014.',
            doi: '10.1093/bioinformatics/btad014'],
        galba: [tool: 'miniprot-boundary-scorer and miniprothint', use: 'scoring of the protein alignments, protein hints',
            ref: 'Brůna T, Li H, Guhlin J, Honsel D, Herbold S, Stanke M, Nenasheva N, Ebel M, Gabriel L, Hoff KJ. Galba: genome annotation with miniprot and AUGUSTUS. BMC Bioinformatics. 2023;24(1):327.',
            doi: '10.1186/s12859-023-05449-z'],
        diamond: [tool: 'DIAMOND', use: 'protein homology of the ORFs, selection of source species in large protein databases',
            ref: 'Buchfink B, Reuter K, Drost HG. Sensitive protein alignments at tree-of-life scale using DIAMOND. Nature Methods. 2021;18(4):366-368.',
            doi: '10.1038/s41592-021-01101-x'],
        hisat2: [tool: 'HISAT2', use: 'alignment of short reads',
            ref: 'Kim D, Paggi JM, Park C, Bennett C, Salzberg SL. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nature Biotechnology. 2019;37(8):907-915.',
            doi: '10.1038/s41587-019-0201-4'],
        minimap2: [tool: 'minimap2', use: 'alignment of Iso-Seq reads',
            ref: 'Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018;34(18):3094-3100.',
            doi: '10.1093/bioinformatics/bty191'],
        samtools: [tool: 'SAMtools', use: 'processing of alignments',
            ref: 'Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. GigaScience. 2021;10(2):giab008.',
            doi: '10.1093/gigascience/giab008'],
        augustus: [tool: 'AUGUSTUS (bam2hints)', use: 'intron hints from read alignments',
            ref: 'Stanke M, Diekhans M, Baertsch R, Haussler D. Using native and syntenically mapped cDNA alignments to improve *de novo* gene finding. Bioinformatics. 2008;24(5):637-644.',
            doi: '10.1093/bioinformatics/btn013'],
        varus: [tool: 'VARUS', use: 'sampled and aligned reads of the pyVARUS input',
            ref: 'Stanke M, Bruhn W, Becker F, Hoff KJ. VARUS: sampling complementary RNA reads from the sequence read archive. BMC Bioinformatics. 2019;20(1):558.',
            doi: '10.1186/s12859-019-3182-x'],
        stringtie: [tool: 'StringTie', use: 'transcript assembly',
            ref: 'Kovaka S, Zimin AV, Pertea GM, Razaghi R, Salzberg SL, Pertea M. Transcriptome assembly from long-read RNA-seq alignments with StringTie2. Genome Biology. 2019;20(1):278.',
            doi: '10.1186/s13059-019-1910-1'],
        gffread: [tool: 'GffRead and GffCompare', use: 'format conversion, protein sequences, comparison with a reference annotation',
            ref: 'Pertea G, Pertea M. GFF Utilities: GffRead and GffCompare. F1000Research. 2020;9:304.',
            doi: '10.12688/f1000research.23297.2'],
        transdecoder: [tool: 'TransDecoder', use: 'ORFs of the assembled transcripts',
            ref: 'Haas BJ, Papanicolaou A, Yassour M, et al. De novo transcript sequence reconstruction from RNA-seq using the Trinity platform for reference generation and analysis. Nature Protocols. 2013;8(8):1494-1512.',
            doi: '10.1038/nprot.2013.084'],
        td2: [tool: 'TD2', use: 'ORFs of the assembled transcripts',
            ref: 'Mao A, Ji HJ, Haas BJ, Salzberg SL, Sommer MJ. TD2: finding protein coding regions in transcripts. bioRxiv. 2025.',
            doi: '10.1101/2025.04.13.648579'],
        psauron: [tool: 'PSAURON', use: 'ORF scores of TD2',
            ref: 'Sommer MJ, Zimin AV, Salzberg SL. PSAURON: a tool for assessing protein annotation across a broad range of species. NAR Genomics and Bioinformatics. 2025;7(1):lqae189.',
            doi: '10.1093/nargab/lqae189'],
        drusilla: [tool: 'Drusilla', use: 'ORFs of the assembled transcripts as high-confidence genes, filter of the gene predictions',
            ref: 'Gabriel L, Hoff KJ. Annotating Eukaryotic Genomes by Combining Deep Learning with extrinsic evidence. Poster, German Conference on Bioinformatics (GCB). 2026.',
            doi: '10.13140/RG.2.2.24444.91521'],
        lightgbm: [tool: 'LightGBM', use: 'filter of the gene predictions',
            ref: 'Ke G, Meng Q, Finley T, Wang T, Chen W, Ma W, Ye Q, Liu TY. LightGBM: A Highly Efficient Gradient Boosting Decision Tree. Advances in Neural Information Processing Systems 30 (NIPS 2017). 2017.',
            url: 'https://proceedings.neurips.cc/paper_files/paper/2017/hash/6449f44a102fde848669bdd9eb6b76fa-Abstract.html'],
        genometools: [tool: 'GenomeTools', use: 'validation of the published GFF3 files',
            ref: 'Gremme G, Steinbiss S, Kurtz S. GenomeTools: a comprehensive software library for efficient processing of structured genome annotations. IEEE/ACM Transactions on Computational Biology and Bioinformatics. 2013;10(3):645-656.',
            doi: '10.1109/TCBB.2013.68'],
        busco: [tool: 'BUSCO', use: 'completeness of genome and proteome',
            ref: 'Manni M, Berkeley MR, Seppey M, Simão FA, Zdobnov EM. BUSCO Update: novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes. Molecular Biology and Evolution. 2021;38(10):4647-4654.',
            doi: '10.1093/molbev/msab199'],
        compleasm: [tool: 'compleasm', use: 'completeness of genome and proteome',
            ref: 'Huang N, Li H. compleasm: a faster and more accurate reimplementation of BUSCO. Bioinformatics. 2023;39(10):btad595.',
            doi: '10.1093/bioinformatics/btad595'],
        omark: [tool: 'OMArk', use: 'consistency and completeness of the proteome',
            ref: 'Nevers Y, Warwick Vesztrocy A, Rossier V, Train CM, Altenhoff A, Dessimoz C, Glover NM. Quality assessment of gene repertoire annotations with OMArk. Nature Biotechnology. 2025;43(1):124-133.',
            doi: '10.1038/s41587-024-02147-w'],
        omamer: [tool: 'OMAmer', use: 'protein families of the proteome for OMArk',
            ref: 'Rossier V, Warwick Vesztrocy A, Robinson-Rechavi M, Dessimoz C. OMAmer: tree-driven and alignment-free protein assignment to subfamilies outperforms closest sequence approaches. Bioinformatics. 2021;37(18):2866-2873.',
            doi: '10.1093/bioinformatics/btab219'],
        // no papers: the citations that the authors ask for (CITATION.cff of
        // pybarrnap; barrnap --citation, release 0.9 that pybarrnap re-implements)
        pybarrnap: [tool: 'pybarrnap', use: 'rRNA genes',
            ref: 'Shimoyama Y. pybarrnap: Python implementation of barrnap. 2024.',
            url: 'https://github.com/moshi4/pybarrnap'],
        barrnap: [tool: 'barrnap', use: 'method and rRNA models that pybarrnap re-implements',
            ref: 'Seemann T. barrnap 0.9: BAsic Rapid Ribosomal RNA Predictor. 2018.',
            url: 'https://github.com/tseemann/barrnap'],
        trnascan: [tool: 'tRNAscan-SE', use: 'tRNA genes',
            ref: 'Chan PP, Lin BY, Mak AJ, Lowe TM. tRNAscan-SE 2.0: improved detection and functional classification of transfer RNA genes. Nucleic Acids Research. 2021;49(16):9077-9096.',
            doi: '10.1093/nar/gkab688'],
        infernal: [tool: 'Infernal', use: 'ncRNA genes of the Rfam families',
            ref: 'Nawrocki EP, Eddy SR. Infernal 1.1: 100-fold faster RNA homology searches. Bioinformatics. 2013;29(22):2933-2935.',
            doi: '10.1093/bioinformatics/btt509'],
        rfam: [tool: 'Rfam', use: 'RNA families database (release 15.1)',
            ref: 'Ontiveros-Palacios N, Cooke E, Nawrocki EP, Triebel S, Marz M, Rivas E, Griffiths-Jones S, Petrov AI, Bateman A, Sweeney B. Rfam 15: RNA families database in 2025. Nucleic Acids Research. 2025;53(D1):D258-D267.',
            doi: '10.1093/nar/gkae1023'],
        feelnc: [tool: 'FEELnc', use: 'lncRNA genes',
            ref: 'Wucher V, Legeai F, Hédan B, et al. FEELnc: a tool for long non-coding RNA annotation and its application to the dog transcriptome. Nucleic Acids Research. 2017;45(8):e57.',
            doi: '10.1093/nar/gkw1306'],
        fantasia: [tool: 'FANTASIA', use: 'GO terms from protein language model embeddings',
            ref: 'Martínez-Redondo GI, Perez-Canales FM, Carbonetto B, Fernández JM, Barrios-Núñez I, Vázquez-Valls M, Cases I, Rojas AM, Fernández R. FANTASIA leverages language models to decode the functional dark proteome across the animal tree of life. Communications Biology. 2025;8(1):1227.',
            doi: '10.1038/s42003-025-08651-2'],
        fantasia_suite: [tool: 'FANTASIA suite (FANTASIA-Lite)', use: 'GO term annotation software',
            ref: 'Pérez-Canales FM, Domínguez-Rodríguez À, Carbonetto B, Fernández R, Cases I, Rojas AM. FANTASIA suite: a reproducible and configurable framework for embedding-based functional annotation of proteins. NAR Genomics and Bioinformatics. 2026;8(3):lqag106.',
            doi: '10.1093/nargab/lqag106'],
        prott5: [tool: 'ProtT5', use: 'protein language model of FANTASIA',
            ref: 'Elnaggar A, Heinzinger M, Dallago C, Rehawi G, Wang Y, Jones L, Gibbs T, Feher T, Angerer C, Steinegger M, Bhowmik D, Rost B. ProtTrans: Toward Understanding the Language of Life Through Self-Supervised Learning. IEEE Transactions on Pattern Analysis and Machine Intelligence. 2022;44(10):7112-7127.',
            doi: '10.1109/TPAMI.2021.3095381'],
    ]
}

// Keys of references() for a run. run: [mode, genefinder (name or null),
// tiberiusModel (params.tiberius.model_cfg or null), hc ('drusilla' |
// 'transdecoder'), orfFinder ('td1' | 'td2'), rescue, odb12, shortFastq,
// shortBam, shortVarus, isoFastq, isoVarus, busco, compleasm, omark,
// gffcompare, ncrna, lncrna, fantasia] (all but the first five boolean).
def citationKeys(Map run) {
    def keys = evidenceCitationKeys(run)
    // post-processing of the final annotation (runs with a gene finder)
    if( run.genefinder ) {
        if( !keys.contains('gffread') ) keys << 'gffread'
        keys << 'genometools'
        if( run.busco ) keys << 'busco'
        if( run.compleasm ) keys << 'compleasm'
        if( run.omark ) keys << 'omark' << 'omamer'
        if( run.ncrna ) {
            keys << 'pybarrnap' << 'barrnap' << 'trnascan' << 'infernal' << 'rfam'
            if( run.lncrna && run.mode in ['rnaseq', 'isoseq', 'mixed'] ) keys << 'feelnc'
        }
        if( run.fantasia ) keys << 'fantasia' << 'fantasia_suite' << 'prott5'
    }
    return keys
}

// Keys of the gene finder and the evidence steps of a run (see citationKeys)
def evidenceCitationKeys(Map run) {
    def keys = ['paludamentum', 'nextflow']
    def transcripts = run.mode in ['rnaseq', 'isoseq', 'mixed']
    def drusilla = transcripts && run.hc == 'drusilla'
    def shortReads = run.mode in ['rnaseq', 'mixed']
    def longReads = run.mode in ['isoseq', 'mixed']

    if( run.genefinder == 'tiberius' ) {
        keys << 'tiberius'
        def model = run.tiberiusModel ? run.tiberiusModel.toString().tokenize('/').last() : ''
        if( !model.startsWith('mammalia') ) keys << 'tiberius_clades'
    }
    if( run.genefinder == 'vipsania' ) keys << 'vipsania'
    if( run.mode == 'abinitio' ) return keys

    if( run.odb12 ) keys << 'orthodb'
    keys << 'miniprot' << 'galba'
    if( (transcripts && !drusilla) || (run.genefinder && run.odb12) ) keys << 'diamond'

    if( shortReads && run.shortFastq ) keys << 'hisat2'
    if( longReads && run.isoFastq ) keys << 'minimap2'
    def aligned = (shortReads && (run.shortFastq || run.shortBam)) || (longReads && run.isoFastq)
    if( aligned || (drusilla && run.rescue) ) keys << 'samtools'
    if( aligned ) keys << 'augustus'
    if( (shortReads && run.shortVarus) || (longReads && run.isoVarus) ) keys << 'varus'
    if( transcripts ) keys << 'stringtie'
    keys << 'gffread'

    if( drusilla ) {
        keys << 'drusilla' << 'lightgbm'
        // the hint rescue runs Tiberius, also when Vipsania is the gene finder
        if( run.rescue && !keys.contains('tiberius') ) keys << 'tiberius'
    } else if( transcripts ) {
        if( run.orfFinder == 'td2' ) keys << 'td2' << 'psauron'
        else keys << 'transdecoder'
    }
    return keys
}

// One reference as a Markdown list item; a reference whose short citation
// has a letter after the year (shortCite) says so, for the methods text.
def citationLine(String key) {
    def r = references()[key]
    def link = r.doi ? "https://doi.org/${r.doi}" : r.url
    def cite = shortCite(key)
    def asCited = cite ==~ /.*\d{4}[a-z]/ ? " (cited as ${cite})" : ''
    return "- **${r.tool}** (${r.use}). ${r.ref} ${link}${asCited}"
}

// ---- short citations (methods.md)

// Surname of an author of a reference string ('Di Tommaso P' -> 'Di Tommaso').
def citeSurname(String author) {
    def parts = author.trim().tokenize(' ')
    return parts.size() > 1 ? parts.subList(0, parts.size() - 1).join(' ') : author.trim()
}

// Year of a reference string (the year after the venue: '. 2017;' or
// '. 2026.'), or 'n.d.' for a reference without one.
def citeYear(String ref) {
    def years = ref.findAll(/\. (?:19|20)\d\d[;.]/)
    return years ? years.last().substring(2, 6) : 'n.d.'
}

// Authors and year of a reference string: 'Li, 2023', 'Nawrocki & Eddy, 2013',
// 'Gabriel et al., 2024'.
def citeAuthorYear(String ref) {
    def end = ref.indexOf('. ')
    def authors = (end > 0 ? ref.substring(0, end) : ref).tokenize(',').collect { a -> a.trim() }.findAll { a -> a }
    def etAl = authors && authors.last() == 'et al'
    if( etAl ) authors = authors.subList(0, authors.size() - 1)
    def first = citeSurname(authors[0])
    def names = etAl || authors.size() > 2 ? "${first} et al." :
        authors.size() == 2 ? "${first} & ${citeSurname(authors[1])}" : first
    return "${names}, ${citeYear(ref)}"
}

// Short citation of a reference of references() for the methods text
// (lib_nf/methods.nf): citeAuthorYear, with a, b, ... after the year in the
// order of references() when references share authors and year
// ('Gabriel & Hoff, 2026a' Paludamentum, 'Gabriel & Hoff, 2026b' Drusilla).
def shortCite(String key) {
    def refs = references()
    if( !refs[key] ) error "shortCite: no reference '${key}' in lib_nf/citations.nf."
    def own = citeAuthorYear(refs[key].ref.toString())
    if( own.endsWith('n.d.') ) return own
    def same = refs.findAll { k, r -> citeAuthorYear(r.ref.toString()) == own }.keySet().toList()
    return same.size() > 1 ? own + ('abcdefghijklmnopqrstuvwxyz'[same.indexOf(key)]) : own
}

// Content of <outdir>/citations.md.
def citationsText(Map run, String version) {
    def lines = [
        '# References of this Paludamentum run',
        '',
        "Paludamentum ${version}, mode `${run.mode}`, gene finder: ${run.genefinder ?: 'none'}.",
        '',
        'This file lists the software and data that this run used, selected from the',
        'inputs and the parameters of the run. Please cite them in a publication',
        'that uses its results.',
        '',
    ]
    citationKeys(run).each { key -> lines << citationLine(key) }
    return lines.join('\n') + '\n'
}
