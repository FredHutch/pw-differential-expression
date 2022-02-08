// Validate the gene count tables
process validate_counts {
    container "${params.container__pandas}"
    label "io_limited"
    
    input:
    // Input file will be placed in the working directory with this name
    path "manifest.csv"
    path "counts.raw"

    output:
    // If validation was successful, the output will be written with this path
    path "counts.csv"

    script:
    // Run the script in bin/validate_counts.py
    """
    validate_counts.py
    """

}


workflow validate {

    main:

        log.info"""${params.manifest}"""

        // Validate the counts file
        validate_counts(
            Channel.fromPath(params.manifest),
            Channel.fromPath(params.counts)
        )

    emit:
    validate_counts.out        

}