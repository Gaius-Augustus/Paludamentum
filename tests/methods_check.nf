// Helper of tests/test_methods.py: the methods text, its cited keys and the
// citation keys of many run maps at once, without running a pipeline.
//
//   nextflow run tests/methods_check.nf --cases cases.json --out result.json
//
// cases.json: a list of run maps as main.nf builds them. result.json: per case
// [text, cited, selected, problems], and the short citations of all references.

include { references; citationKeys } from '../lib_nf/citations.nf'
include { methodsText; methodsParagraphs; methodsCitationProblems; shortCite } from '../lib_nf/methods.nf'

workflow {
    main:
    def cases = new groovy.json.JsonSlurper().parse(file(params.cases.toString()).toFile())
    def results = cases.collect { run ->
        [
            text: methodsText(run, '9.9.9'),
            cited: methodsParagraphs(run, '9.9.9').keys,
            selected: citationKeys(run),
            problems: methodsCitationProblems(run),
        ]
    }
    def refs = references().collectEntries { key, r -> [(key): [ref: r.ref, short: shortCite(key)]] }
    file(params.out.toString()).text = groovy.json.JsonOutput.toJson([cases: results, references: refs])
}
