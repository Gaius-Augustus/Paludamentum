// References of the tools that the pipeline runs, and the citation file of a run.
//
// The DOIs were checked against the metadata registered at doi.org. The same
// list is the References section of the README; tests/test_launcher.py checks
// that the two agree.

// key: [tool, use (what the pipeline does with it), ref, doi or url]
def references() {
    return [
        paludamentum: [tool: 'Paludamentum', use: 'the pipeline',
            ref: 'Gabriel L, Hoff KJ. Paludamentum: evidence integration pipeline for the Gaius-Augustus gene finders.',
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
        gffread: [tool: 'GffRead', use: 'format conversion, protein sequences',
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
    ]
}

// Keys of references() for a run. run: [mode, genefinder (name or null),
// tiberiusModel (params.tiberius.model_cfg or null), hc ('drusilla' |
// 'transdecoder'), orfFinder ('td1' | 'td2'), rescue, odb12, shortFastq,
// shortBam, shortVarus, isoFastq, isoVarus] (all but the first five boolean).
def citationKeys(Map run) {
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

// One reference as a Markdown list item.
def citationLine(Map r) {
    def link = r.doi ? "https://doi.org/${r.doi}" : r.url
    return "- **${r.tool}** (${r.use}). ${r.ref} ${link}"
}

// Content of <outdir>/citations.md.
def citationsText(Map run, String version) {
    def refs = references()
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
    citationKeys(run).each { key -> lines << citationLine(refs[key]) }
    return lines.join('\n') + '\n'
}
