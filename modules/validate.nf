// Validate the gene count tables
process validate_tables {
    container "${params.container__pandas}"
    label "io_limited"
    
    input:
    // Input file will be placed in the working directory with this name
    tuple path("manifest.csv"), path("${gene_count_filename}.raw"), val(gene_count_filename)

    output:
    // If validation was successful, the output will be written with this path
    path "${gene_count_filename}"

    script:
    // Run the script in bin/validate_tables.py
    """
    validate_tables.py
    """

}


workflow validate {

    main:

        // Set up a channel which will be used to place input files
        Channel
            .empty()
            .set { gene_counts }

        // If the parameter `salmon_merged_counts` was provided
        if ( params.salmon_merged_counts ){

            // Add that file to the channel for validation
            gene_counts.bind([
                file(params.manifest),
                file(params.salmon_merged_counts),
                "salmon.merged.gene_counts.tsv"
            ])

        }

        // Validate any of the files found
        validate_tables(
            Channel
                .fromPath(params.salmon_merged_counts)
        )

    emit:
    validate_tables.out        

}