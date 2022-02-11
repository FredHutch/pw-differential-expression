// Filter genes with the filterbyExpr package
process filter {
    container "${params.container__edgeR}"
    label "io_limited"
    
    input:
    tuple path(manifest), path("raw.counts.csv")

    output:
    tuple path("${manifest.name}"), path("counts.csv")

    script:
    template "filterbyExpr.R"

}

// Run the DESeq2 algorithm
process deseq2 {
    container "${params.container__deseq2}"
    label "mem_medium"
    publishDir "${params.output_folder}", mode: "copy", overwrite: true
    
    input:
    // Input file will be placed in the working directory with this name
    tuple path(manifest), path(counts)

    output:
    // If validation was successful, the output will be written with this path
    path "*.DEseq2.csv"

    script:
    // Run the script in templates/run_deseq2.R
    template "run_deseq2.R"

}

// Run the edgeR algorithm
process edgeR {
    container "${params.container__edgeR}"
    label "mem_medium"
    publishDir "${params.output_folder}", mode: "copy", overwrite: true
    
    input:
    // Input file will be placed in the working directory with this name
    tuple path(manifest), path(counts)

    output:
    // If validation was successful, the output will be written with this path
    path "*.edgeR.csv"

    script:
    // Run the script in templates/run_edgeR.R
    template "run_edgeR.R"

}

workflow test {
    take:
    // Table of gene counts paired with the manifest, 
    // validated to conform to the same order of specimens in each
    counts_ch

    main:

    // Filter the counts table with filterbyExpr
    filter(counts_ch)

    // The statistical test applied to the data will be determined
    // by the parameters selected by the user
    if ( params.algorithm == "deseq2" ){
        
        deseq2(filter.out)
        csv = deseq2.out

    } else if ( params.algorithm == "edgeR" ){
        
        edgeR(filter.out)
        csv = edgeR.out 

    } else if ( params.algorithm == "limma_voom" ){
        
        limma_voom(filter.out)
        csv = limma_voom.out

    } else {
        throw new Exception("""
    ERROR:
    Algorithm not recognized: ${params.algorithm}
    Supported options: deseq2, edgeR, limma_voom
        """)
    }

    emit:
    csv
}