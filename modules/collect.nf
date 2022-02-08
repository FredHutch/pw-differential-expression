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

workflow collect {
    take:
    // A collection of CSVs with results from a differential expression test
    csv_ch

    main:

    // Collect all of the results
    all(csv_ch.toSortedList())

}