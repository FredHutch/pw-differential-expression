// Using DSL-2
nextflow.enable.dsl=2

// Import modules
include { validate } from './modules/validate'
include { test } from './modules/test'
include { collect } from './modules/collect'

// Main workflow
workflow {
    validate | test | collect
}