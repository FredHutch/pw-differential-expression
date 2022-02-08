// Using DSL-2
nextflow.enable.dsl=2

// Import modules
include { validate } from './modules/validate'
include { test } from './modules/test'
include { collect } from './modules/collect'

// Main workflow
workflow {

    // Validate the contents of --counts and align the
    // column order with rows in --manifest
    validate()

    // Run the indicated test library on the counts table
    test(
        validate.out
    )

    collect(
        test.out
    )
}