// Collect all of the results into a single table
process all {
    container "${params.container__pandas}"
    label "io_limited"
    publishDir "${params.output_folder}", mode: "copy", overwrite: true
    
    input:
    path "*"

    output:
    path "DE_results.csv"

    script:
    // Run the script in templates/collect_all.py
    template "collect_all.py"

}

process anndata {
    container "${params.container__pandas}"
    label "io_limited"
    publishDir "${params.output_folder}", mode: "copy", overwrite: true, pattern: "*.h5ad"
    publishDir "${params.output_folder}", mode: "copy", overwrite: true, pattern: "*.csv"
    publishDir "${params.web_folder}", mode: "copy", overwrite: true, pattern: "*.zarr", enabled: "${params.web_folder}" != "false"
    publishDir "${params.web_folder}", mode: "copy", overwrite: true, pattern: "*.vt.json", enabled: "${params.web_folder}" != "false"

    input:
    path "DE_results.csv"
    tuple path("manifest.csv"), path("counts.csv")

    output:
    path "*.h5ad", emit: h5ad
    path "*.zarr", hidden: true, emit: zarr
    path "*.vt.json", emit: vt_json
    path "*.csv", emit: csv

    """#!/bin/bash
set -e
make_anndata.py
    """

}

process manifest {
    publishDir "${params.web_folder}", mode: "copy", overwrite: true, enabled: "${params.web_folder}" != "false"
    container "${params.container__pandas}"
    label "io_limited"

    input:
    path "*"

    output:
    path "chart.manifest.json"

    """#!/bin/bash
set -e
chart_manifest.py
    """

}

workflow collect {
    take:
    // A collection of CSVs with results from a differential expression test
    results_csv_ch
    // The filtered counts and manifest used to run the tests
    filtered_ch

    main:

    // Collect all of the results
    all(results_csv_ch.toSortedList())

    // Format as AnnData
    anndata(all.out.toSortedList(), filtered_ch)

    // Format the chart.manifest.json
    manifest(anndata.out.vt_json.toSortedList())

}