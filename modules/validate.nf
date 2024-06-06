// Validate the metadata table, and reformat it as appropriate
// to drive downstream comparisons
process manifest {

    container "${params.container__pandas}"
    label "io_limited"
    publishDir "${params.output_folder}/manifest/", mode: "copy", overwrite: true
    
    input:
    // Input file will be placed in the working directory with this name
    path "input_manifest.csv"

    output:
    // The output file(s) will contain the comparison column name in the file name
    path "*.manifest.csv", emit: for_de
    path "manifest.csv", emit: full

    script:
    // Run the script in templates/validate_manifest.py
    template "validate_manifest.py"

}

// Validate the gene count tables
process counts {
    container "${params.container__pandas}"
    label "io_limited"
    
    input:
    // Input file will be placed in the working directory with this name
    tuple path(counts_table), path(manifest_table)

    output:
    // If validation was successful, the output will be written with this path
    tuple path("validated.${manifest_table.name}"), path("counts.csv")

    script:
    // Run the script in templates/validate_counts.py
    template "validate_counts.py"

}


workflow validate {

    main:

        // Make sure that an output folder was defined
        if ( params.output_folder == false ) {
            throw new Exception("""Must specify parameter: output_folder""")
        }

        // Make sure that a comparison column was defined
        if ( params.comp_col == "" ) {
            throw new Exception("""Must specify parameter: comp_col""")
        }

        // Validate the contents of the manifest
        
        // If a categorical comparison was defined, split up
        // the manifest to yield each of the appropriate pairwise
        // comparisons.
        
        // If a filtering expression was specified, apply that filtering
        // before performing any additional transformations.
        manifest(
            Channel.fromPath("${params.manifest}", checkIfExists: true)
        )

        // Validate the counts file
        counts(
            Channel
                .fromPath("${params.counts}")
                .combine(
                    manifest.out.for_de.flatten()
                )
        )

    emit:
    counts.out        

}