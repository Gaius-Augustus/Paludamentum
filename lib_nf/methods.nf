// Methods text of a run: <outdir>/methods.md, a paper-style description of
// what this run did, in execution order, with inline citations.
//
// The text is built from the same run map as citations.md (main.nf), and its
// citations are the short forms of references() (shortCite in lib_nf/citations.nf), so
// that the two files never disagree. Every key that the text cites must be one
// of citationKeys(run); methodsCitationProblems(run) lists the ones that are
// not, and tests/test_methods.py checks it for many runs.
//
// Format (read by bin/paludamentum_report.py): the first line
// '# Methods of this Paludamentum run', a blank line, then paragraphs
// separated by blank lines, one line each. Inline Markdown only (*italic*,
// **bold**, `code`); no lists, no further headings, no links.

include { references; citationKeys; shortCite; withProteins } from './citations.nf'

// ---- citations

// '(A, 2020; B et al., 2021)' for the keys, which are noted in ctx.used.
def cite(Map ctx, List keys) {
    keys.each { k -> if( !ctx.used.contains(k) ) ctx.used << k }
    return '(' + keys.collect { k -> shortCite(k) }.join('; ') + ')'
}

// ---- text helpers

// 'a', 'a and b', 'a, b and c'
def joinAnd(List items) {
    return joinWith(items, 'and')
}

// 'a', 'a or b', 'a, b or c'
def joinOr(List items) {
    return joinWith(items, 'or')
}

def joinWith(List items, String word) {
    def xs = items.findAll { x -> x }.collect { x -> x.toString() }
    if( xs.size() < 2 ) return xs ? xs[0] : ''
    return xs.subList(0, xs.size() - 1).join(', ') + " ${word} " + xs.last()
}

// A length in bp: '25,000 bp', '20 Mb'
def bpText(value) {
    def n = value as Long
    if( n >= 1000000 && n % 1000000 == 0 ) return "${n.intdiv(1000000)} Mb"
    return "${String.format('%,d', n)} bp"
}

// A number as written in the params (0.5, 3, 1e-5)
def numText(value) {
    return value == null ? '' : value.toString()
}

def capitalized(String s) {
    return s ? s.substring(0, 1).toUpperCase() + s.substring(1) : s
}

// ---- paragraphs

// Pipeline, gene finder and the ab initio predictions.
def methodsGeneFinder(Map ctx, Map run, String version) {
    def s = []
    def evidence = withProteins(run) ? [
        abinitio: 'without extrinsic evidence',
        proteins: 'with protein evidence',
        rnaseq:   'with protein and short-read RNA-Seq evidence',
        isoseq:   'with protein and Iso-Seq long-read evidence',
        mixed:    'with protein, short-read RNA-Seq and Iso-Seq long-read evidence',
    ][run.mode] : [
        abinitio: 'without extrinsic evidence',
        rnaseq:   'with short-read RNA-Seq evidence and without protein evidence',
        isoseq:   'with Iso-Seq long-read evidence and without protein evidence',
        mixed:    'with short-read RNA-Seq and Iso-Seq long-read evidence and without protein evidence',
    ][run.mode]
    if( run.mode in ['rnaseq', 'isoseq', 'mixed'] && run.stringtieFiles && !transcriptReads(run) )
        evidence = withProteins(run) ? 'with protein evidence and transcript assemblies' :
                                       'with transcript assemblies and without protein evidence'
    s << "The genome was annotated with Paludamentum ${version} ${cite(ctx, ['paludamentum'])}, " +
         "a pipeline implemented in Nextflow ${cite(ctx, ['nextflow'])}, in mode `${run.mode}` (${evidence})."

    if( !run.genefinder ) {
        s << 'No gene finder was run: the run prepared the extrinsic evidence and the high-confidence genes ' +
             'described below, but no final annotation.'
        return s.join(' ')
    }

    def clade = run.clade ? ", a model of the clade ${run.clade}" : ''
    if( run.genefinder == 'both' ) {
        def tibModel = run.tiberiusModel ? "`${file(run.tiberiusModel.toString()).baseName}`" : null
        def vipModel = run.model ? "`${run.model}`" : null
        if( run.tiberiusResult )
            s << "An existing Tiberius ${cite(ctx, ['tiberius'])} prediction (`${run.tiberiusResult}`) was reused " +
                 'instead of running the gene finder.'
        else {
            def clades = tibModel && run.tiberiusModel && !run.tiberiusModel.toString().startsWith('mammalia') ?
                " ${cite(ctx, ['tiberius_clades'])}" : ''
            s << "Genes were predicted *ab initio* with Tiberius ${cite(ctx, ['tiberius'])}" +
                 (tibModel ? " using the model ${tibModel}${clades}" : '') + '.'
        }
        if( run.vipsaniaResult )
            s << "An existing Vipsania ${cite(ctx, ['vipsania'])} prediction (`${run.vipsaniaResult}`) was reused."
        else {
            s << "A second set of *ab initio* predictions was made with Vipsania ${cite(ctx, ['vipsania'])}" +
                 (vipModel ? " using the pretrained model ${vipModel}" : '') + '.'
            if( run.finetune )
                s << 'Vipsania was first finetuned on the target genome, which it then annotated as a whole.'
        }
        return s.join(' ')
    }

    def tool = run.genefinder == 'tiberius' ? 'Tiberius' : 'Vipsania'
    def model = run.model ? "`${run.model}`" : null
    def cladeTxt = run.clade && run.model && run.clade.toString().toLowerCase() != run.model.toString().toLowerCase() ?
        ", a model of the clade ${run.clade}" : ''
    def resultFile = run.genefinder == 'tiberius' ? run.tiberiusResult : run.vipsaniaResult
    if( resultFile ) {
        s << "An existing prediction of ${tool} ${cite(ctx, [run.genefinder])}, `${resultFile}`" +
             (model ? ", made with the model ${model}${cladeTxt}," : '') +
             ' was used as the *ab initio* gene set instead of running the gene finder; ' +
             'its genes and transcripts were renumbered.'
        return s.join(' ')
    }
    if( run.genefinder == 'tiberius' ) {
        def clades = model && !run.model.toString().startsWith('mammalia') ? " ${cite(ctx, ['tiberius_clades'])}" : ''
        s << "Genes were predicted *ab initio* with Tiberius ${cite(ctx, ['tiberius'])}" +
             (model ? " using the model ${model}${clades}${cladeTxt}" : '') + '.'
    } else {
        s << "Genes were predicted *ab initio* with Vipsania ${cite(ctx, ['vipsania'])}" +
             (model ? " using the pretrained model ${model}${cladeTxt}" : '') + '.'
        if( run.finetune )
            s << 'The model was first finetuned on the target genome, which Vipsania then annotated as a whole.'
    }
    return s.join(' ')
}

// True if the run aligned or brought reads (not only StringTie assemblies).
def transcriptReads(Map run) {
    def shortReads = run.mode in ['rnaseq', 'mixed'] && (run.shortFastq || run.shortBam || run.shortVarus)
    def longReads = run.mode in ['isoseq', 'mixed'] && (run.isoFastq || run.isoVarus)
    return shortReads || longReads || run.mixVarus
}

// Protein database, miniprot, miniprothint.
def methodsProteins(Map ctx, Map run) {
    def s = []
    def sources = []
    if( run.proteinFiles == 1 ) sources << 'the provided protein sequences'
    else if( run.proteinFiles > 1 ) sources << "${run.proteinFiles} provided protein files"
    if( run.odb12 ) {
        def parts = run.odb12Partitions ?: []
        sources << "the OrthoDB v12 partition${parts.size() == 1 ? '' : 's'} " +
                   joinAnd(parts.collect { p -> "*${p}*" }) + " ${cite(ctx, ['orthodb'])}"
    }
    s << "The protein database consisted of ${joinAnd(sources)}" +
         (sources.size() > 1 ? ', concatenated into one file.' : '.')
    if( run.genefinder && run.odb12 )
        s << 'If the database held more than one million proteins, all with OrthoDB identifiers, it was first ' +
             'reduced to the proteins of 13 species: the *ab initio* proteins were searched against it with DIAMOND ' +
             "${cite(ctx, ['diamond'])} in `blastp --very-sensitive` mode with an e-value of at most 1e-5, and the 13 species with " +
             'the most best hits (identity at least 30 %, query coverage at least 50 %) were kept.'
    def matrix = run.scoringMatrix ? " with the scoring matrix `${run.scoringMatrix}`" : ''
    s << "The proteins were aligned to the genome with miniprot ${cite(ctx, ['miniprot'])}, and the alignments " +
         "were scored with miniprot-boundary-scorer${matrix}. miniprothint ${cite(ctx, ['galba'])} filtered the " +
         'scored alignments (`--topNperSeed 10 --minScoreFraction 0.5`) and derived high-confidence hints from them, ' +
         'and intron and CDS-part hints were made from its alignments with `aln2hints.pl`.'
    return s.join(' ')
}

// Read alignment, intron hints, StringTie assemblies.
def methodsTranscripts(Map ctx, Map run) {
    def s = []
    def shortReads = run.mode in ['rnaseq', 'mixed']
    def longReads = run.mode in ['isoseq', 'mixed']
    // Evidence modes always use ONE StringTie assembly of all reads (the
    // combine_genesets flow). Keeps the same assembly rule the Drusilla flow
    // used before the unification.
    def oneAssembly = run.mode in ['rnaseq', 'isoseq', 'mixed']
    def aligned = []

    if( shortReads && run.shortFastq ) {
        def kinds = joinAnd([run.shortPaired ? 'paired-end' : null, run.shortSingle ? 'single-end' : null])
        if( run.shortSra )
            s << 'RNA-Seq runs were downloaded from the NCBI Sequence Read Archive with `prefetch` and `fasterq-dump`.'
        s << "Short ${kinds} RNA-Seq reads were aligned to the genome with HISAT2 ${cite(ctx, ['hisat2'])}, " +
             "and the alignments were sorted with SAMtools ${cite(ctx, ['samtools'])}."
        aligned << 'short-read'
    }
    if( shortReads && run.shortBam ) {
        s << 'Provided short-read alignments (BAM) were used as provided, sorted by coordinate with SAMtools' +
             (ctx.used.contains('samtools') ? '' : " ${cite(ctx, ['samtools'])}") +
             ' where their header did not declare that order.'
    }
    if( longReads && run.isoFastq ) {
        if( run.isoSra )
            s << 'Iso-Seq runs were downloaded from the NCBI Sequence Read Archive with `prefetch` and `fasterq-dump`.'
        s << "Iso-Seq reads were aligned to the genome with minimap2 ${cite(ctx, ['minimap2'])} " +
             'with `-ax splice:hq -uf`, and the alignments were sorted with SAMtools' +
             (ctx.used.contains('samtools') ? '.' : " ${cite(ctx, ['samtools'])}.")
        aligned << 'Iso-Seq'
    }
    if( aligned )
        s << "${capitalized(joinAnd(aligned))} libraries with an alignment rate below ${numText(run.minAlignmentRate)} % " +
             '(`samtools flagstat`) were discarded.'
    if( shortReads && (run.shortFastq || run.shortBam) && run.strandFix )
        s << 'The library strandedness of every short-read alignment file was inferred from its spliced reads ' +
             '(read orientation against the XS strand tag of the aligner), and the XS tag was set on all reads of ' +
             'stranded libraries (`drusilla strand-fix-bam`), so that StringTie could assign a strand to ' +
             'single-exon transcripts and separate overlapping antisense transcripts.'
    def bams = (shortReads && (run.shortFastq || run.shortBam)) || (longReads && run.isoFastq)
    if( bams )
        s << 'The alignments of all libraries of a read type were merged with SAMtools, and intron hints were ' +
             "extracted from them with bam2hints of AUGUSTUS ${cite(ctx, ['augustus'])} and " +
             '`filterIntronsFindStrand.pl`, which assigns the strand of an intron from its splice site dinucleotides.'
    def varus = []
    if( shortReads && run.shortVarus ) varus << 'short reads'
    if( longReads && run.isoVarus ) varus << 'Iso-Seq reads'
    if( varus )
        s << "For the ${joinAnd(varus)}, pyVARUS output was used: reads that VARUS ${cite(ctx, ['varus'])} had " +
             'sampled from the Sequence Read Archive and aligned to the genome, with the StringTie assembly and ' +
             'the intron hints that pyVARUS made from them; the genome MD5 checksum of each pyVARUS run was ' +
             'checked against the genome of this run, and intron hints of several sources were summed per intron.'

    // StringTie: one merged assembly of all reads (combine_genesets flow)
    def st = cite(ctx, ['stringtie'])
    if( oneAssembly ) {
        if( run.stringtieFiles )
            s << "The provided StringTie ${st} assembly of all reads was used as the transcript assembly" +
                 (transcriptReads(run) ? '; the read alignments only gave the coverage track and the intron hints.' :
                                         '; no reads were aligned.')
        else if( run.mode == 'mixed' && run.mixVarus )
            s << "The transcript assembly was the StringTie ${st} assembly of the short-read and the Iso-Seq " +
                 'alignments of pyVARUS together (`varus assemble --short --long`).'
        else if( run.mode == 'mixed' )
            s << "The merged short-read and Iso-Seq alignments were assembled together into one transcript set with " +
                 "StringTie ${st}, using `stringtie --mix`."
        else if( (run.mode == 'rnaseq' && run.shortVarus) || (run.mode == 'isoseq' && run.isoVarus) )
            s << "The StringTie ${st} assembly of the pyVARUS run was the transcript assembly."
        else if( run.mode == 'rnaseq' )
            s << "The merged short-read alignments were assembled into one transcript set with StringTie ${st}."
        else
            s << "The merged Iso-Seq alignments were assembled into one transcript set with StringTie ${st}, using `stringtie -L`."
    }
    if( !bams && !varus )
        s << (withProteins(run) ? 'Without reads, the hints of the run were protein hints only.' :
                                  'Without reads and proteins, the run had no hints.')
    return s.join(' ')
}

// High-confidence genes and their merge with the ab initio predictions.
def methodsHighConfidence(Map ctx, Map run) {
    def s = []
    def transcripts = run.mode in ['rnaseq', 'isoseq', 'mixed']
    def stopCheck = 'The CDS of the HC genes were extended to include the stop codon, and genes without a valid ' +
                    'terminal stop codon or with an in-frame internal stop codon were removed.'
    if( !transcripts ) {
        s << 'Without transcript evidence, the high-confidence (HC) genes were the gene structures of the ' +
             "high-confidence protein alignments selected by miniprothint ${cite(ctx, ['galba'])}."
        s << stopCheck
        return s.join(' ')
    }
    // Unified HC flow (combine_genesets): both ORF callers + the four
    // DIAMOND searches + bedtools genomecov coverage + bam2hints introns.
    def dmd = cite(ctx, ['diamond'])
    def dr  = cite(ctx, ['drusilla'])
    def td2 = cite(ctx, ['td2'])
    def psa = cite(ctx, ['psauron'])
    s << "High-confidence (HC) genes were derived from two ORF callers on the StringTie transcript assembly: " +
         "TD2 ${td2} (which scores candidate ORFs with PSAURON ${psa}" +
         (run.td2PredictArgs ? ", `TD2.Predict ${run.td2PredictArgs}`" : '') + ") and Drusilla ${dr}."
    s << "For TD2, 5' partial and internal ORFs were also trimmed to their first methionine, both versions were " +
         "searched against the protein database with DIAMOND ${dmd} blastp (five targets per query), and " +
         'the comparison of their hits decided whether an ORF was revised to the trimmed, complete start; the ' +
         'revised ORFs were converted back to genome coordinates with `cdna_alignment_orf_to_genome_orf.pl` of ' +
         "TransDecoder ${cite(ctx, ['transdecoder'])}."
    def dru = []
    if( run.drusillaBothStrands )
        dru << 'scanning both strands of transcripts without an assembled strand'
    if( run.drusillaFlank )
        dru << "completing ORFs that reached a transcript end with up to ${numText(run.drusillaFlank)} bp of genomic sequence"
    if( run.drusillaOrfsPerTx != null ) {
        def nOrfs = run.drusillaOrfsPerTx as Integer
        if( nOrfs == 0 )     dru << 'reporting every ORF of a transcript'
        else if( nOrfs > 1 ) dru << "reporting up to ${numText(nOrfs)} ORFs per transcript"
    }
    s << 'For Drusilla, assembled transcripts were pre-filtered by length, StringTie coverage and TPM; Drusilla ' +
         'predicted their ORFs' + (dru ? " (${joinAnd(dru)})" : '') +
         ', and isoforms whose CDS was a subsequence of another were collapsed.'
    def insect = run.clade?.toString()?.toLowerCase() in ['insecta', 'lepidoptera']
    def dropFlag = insect ? ' (plus `--drop-lone-single-exon tiberius`)' : ''
    // the coverage track needs a BAM; pyVARUS directories and a given assembly without reads have none
    def bams = (run.mode in ['rnaseq', 'mixed'] && (run.shortFastq || run.shortBam)) ||
               (run.mode in ['isoseq', 'mixed'] && run.isoFastq)
    def varusHints = (run.mode in ['rnaseq', 'mixed'] && run.shortVarus) || (run.mode in ['isoseq', 'mixed'] && run.isoVarus)
    def rna = bams ? ', bedtools genomecov produced a per-base read-depth bedGraph of the merged read alignments, and ' +
                     'bam2hints produced intron hints' :
              varusHints ? ', and the intron hints of the pyVARUS runs were used (without alignments there was no coverage track)' :
                           ' (without reads, no coverage track and no intron hints entered the integration)'
    s << "The four gene sets (Tiberius, Vipsania, revised TD2, Drusilla) were integrated by `combine_genesets.py`: " +
         "DIAMOND blastp was run on each gene set's proteins (14-column outfmt, including query and subject lengths)" +
         rna + '. combine_genesets.py merged transcripts with identical CDS, clustered them into loci ' +
         'on the same strand, tiered each unique structure from its supporting tool set, adjusted the tier with ' +
         ((bams || varusHints) ? 'RNA, protein' : 'protein') + ' and DIAMOND features, and kept each locus whose best ' +
         "model reached a tier of ${numText(run.minLocus)}${dropFlag}; " +
         'single-exon Tiberius-only models were demoted on insect runs. The result is the final annotation (no ' +
         'separate merge with the *ab initio* predictions is needed, because they already enter the integration).'
    s << stopCheck
    return s.join(' ')
}

// Sanity filter, UTRs, published files.
def methodsPostprocessing(Map ctx, Map run) {
    def s = []
    def transcripts = run.mode in ['rnaseq', 'isoseq', 'mixed']
    def checks = 'transcripts without CDS, with a CDS length (less the phase of a 5\'-partial start) that is not ' +
                 'a multiple of three, with overlapping CDS segments or CDS segments outside their exons, or with ' +
                 'an internal stop codon'
    if( run.sanityFilter )
        s << "A sanity filter removed ${checks}, and genes left without transcripts; a CDS that ended right before " +
             'a stop codon was extended to include it.'
    else
        s << "A sanity filter flagged ${checks}, without changing the annotation."
    if( transcripts && run.utr ) {
        // One merged StringTie assembly feeds UTRs (no separate long-read assembly in the combined flow)
        def longread = false
        s << 'Coding transcripts without UTRs received UTRs from the StringTie assemblies of the run: a multi-exon ' +
             'transcript matched an assembled transcript that contained all its CDS introns, a single-exon ' +
             'transcript an overlapping assembled transcript without an intron in its CDS, and of several matches ' +
             'the longest was used' + (longread ? ', a match in an Iso-Seq assembly taking precedence over one in a short-read assembly' : '') +
             '. The parts of the matched exons outside the CDS became the UTRs, which reached at most ' +
             "${bpText(run.maxUtrExtension)} beyond the CDS and stopped at the nearest gene on the same strand " +
             'that had its own match or UTRs; no CDS was changed.'
    }
    s << "The annotation was written as GFF3 and validated with `gt gff3validator` of GenomeTools " +
         "${cite(ctx, ['genometools'])}, and GffRead ${cite(ctx, ['gffread'])} converted it to GTF and extracted the " +
         'protein and coding sequences of all transcripts.'
    return s.join(' ')
}

// Completeness, hint support, OMArk, gffcompare.
def methodsQuality(Map ctx, Map run) {
    def s = []
    def tools = []
    if( run.compleasm ) tools << "compleasm ${cite(ctx, ['compleasm'])}"
    if( run.busco ) tools << "BUSCO ${cite(ctx, ['busco'])}"
    if( tools )
        s << "The completeness of the genome and of the proteome (the longest coding isoform of each gene) was " +
             "assessed with ${joinAnd(tools)} against the lineage `${run.buscoLineage}`."
    def hintKinds = joinAnd([transcriptReads(run) ? 'RNA-Seq' : null, withProteins(run) ? 'protein' : null])
    if( run.geneSupport && run.mode != 'abinitio' && hintKinds )
        s << "For each coding transcript, the CDS introns and CDS segments supported by ${hintKinds} hints were counted."
    if( run.statistics )
        s << 'Gene structure statistics of the annotation were computed.'
    if( run.omark )
        s << "The consistency and completeness of the proteome were assessed with OMArk ${cite(ctx, ['omark'])} " +
             "after OMAmer ${cite(ctx, ['omamer'])} had placed the proteins of all isoforms into gene families."
    if( run.gffcompare )
        s << 'The coding sequences of the annotation were compared with those of the reference annotation using ' +
             "GffCompare ${cite(ctx, ['gffread'])} with `--strict-match -e 3 -T` on the CDS features."
    return s.join(' ')
}

// rRNA, tRNA, Rfam and lncRNA genes.
def methodsNcrna(Map ctx, Map run) {
    def s = []
    def lncrna = run.lncrna && run.mode in ['rnaseq', 'isoseq', 'mixed']
    s << "Non-coding RNA genes were annotated as follows: rRNA genes with pybarrnap ${cite(ctx, ['pybarrnap'])}, " +
         "the Python implementation of barrnap 0.9 ${cite(ctx, ['barrnap'])}, with `--kingdom euk`, tRNA genes with tRNAscan-SE " +
         "${cite(ctx, ['trnascan'])} in eukaryotic mode" +
         (run.trnascanHighConfidence ? ' with the EukHighConfidenceFilter' : '') +
         ", and other ncRNA genes with `cmscan` of Infernal ${cite(ctx, ['infernal'])} against the covariance " +
         "models of Rfam 15.1 ${cite(ctx, ['rfam'])} with their gathering thresholds, keeping the best-scoring " +
         'hit where hits of one Rfam clan overlapped.'
    if( lncrna )
        s << "Long non-coding RNAs were identified with FEELnc ${cite(ctx, ['feelnc'])}: merged StringTie " +
             'transcripts of at least 200 bp with more than one exon that did not overlap the annotation were ' +
             'the candidates, and those without coding potential according to a random forest trained on the ' +
             'annotated mRNAs and shuffled copies of them (`--mode=shuffle`) were kept; the classifier needs at ' +
             'least 100 candidates and 100 annotated transcripts to be trained.'
    def types = lncrna ? 'rRNA, tRNA, Rfam and lncRNA genes' : 'rRNA, tRNA and Rfam genes'
    s << "The ${types} were added to the annotation in this order of priority; an ncRNA gene was dropped if more " +
         'than half of its exon length overlapped coding exons on the same strand, or more than half of its span ' +
         'an ncRNA gene of higher priority.'
    return s.join(' ')
}

// GO terms with FANTASIA-Lite.
def methodsGoTerms(Map ctx, Map run) {
    return "Gene Ontology terms were assigned with FANTASIA ${cite(ctx, ['fantasia'])}, as implemented in " +
           "FANTASIA-Lite of the FANTASIA suite ${cite(ctx, ['fantasia_suite'])}: ProtT5 ${cite(ctx, ['prott5'])} " +
           'embeddings of all proteins were compared with the FANTASIA lookup table, and GO terms with a score ' +
           "of at least ${numText(run.fantasiaMinScore)} were added to the mRNAs and genes of the GFF3 as " +
           '`Ontology_term` attributes.'
}

// [paragraphs: the paragraphs of the methods, keys: the cited keys of references()]
def methodsParagraphs(Map run, String version) {
    def ctx = [used: []]
    def ps = [methodsGeneFinder(ctx, run, version)]
    if( run.mode != 'abinitio' ) {
        if( withProteins(run) ) ps << methodsProteins(ctx, run)
        if( run.mode in ['rnaseq', 'isoseq', 'mixed'] ) ps << methodsTranscripts(ctx, run)
        ps << methodsHighConfidence(ctx, run)
    }
    if( run.genefinder ) {
        ps << methodsPostprocessing(ctx, run)
        def qc = methodsQuality(ctx, run)
        if( qc ) ps << qc
        if( run.ncrna ) ps << methodsNcrna(ctx, run)
        if( run.fantasia ) ps << methodsGoTerms(ctx, run)
    }
    return [paragraphs: ps, keys: ctx.used]
}

// Cited keys of the methods of a run that citationKeys(run) does not select
// (empty for a correct text).
def methodsCitationProblems(Map run) {
    def selected = citationKeys(run)
    return methodsParagraphs(run, 'test').keys.findAll { k -> !selected.contains(k) }
}

// Content of <outdir>/methods.md.
def methodsText(Map run, String version) {
    def lines = ['# Methods of this Paludamentum run', '']
    lines << methodsParagraphs(run, version).paragraphs.join('\n\n')
    return lines.join('\n') + '\n'
}
