// Using DSL-2
nextflow.enable.dsl=2

// Import modules
include { validate } from './modules/validate'
include { test } from './modules/test'
include { collect } from './modules/collect'

// Main workflow
workflow {

    log.info"""
Differential Expression

    manifest:           ${params.manifest}
    counts:             ${params.counts}
    algorithm:          ${params.algorithm}
    comp_col:           ${params.comp_col}
    comp_ref:           ${params.comp_ref}
    group_cols:         ${params.group_cols}
    filter:             ${params.filter}
    output_folder:      ${params.output_folder}
    web_folder:         ${params.web_folder}
    min_count:          ${params.min_count}
    min_total_count:    ${params.min_total_count}
    large_n:            ${params.large_n}
    min_prop:           ${params.min_prop}
    fdr_method:         ${params.fdr_method}
    container__pandas:  ${params.container__pandas}
    container__deseq2:  ${params.container__deseq2}
    container__edgeR:   ${params.container__edgeR}
    """

    // Validate the contents of --counts and align the
    // column order with rows in --manifest
    validate()

    // Run the indicated test library on the counts table
    test(
        validate.out
    )

    collect(
        test.out.results,
        test.out.filtered
    )
}