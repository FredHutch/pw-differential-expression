// Filter genes with the filterbyExpr package
process filter {
    container "${params.container__edger}"
    label "io_limited"
    
    input:
    tuple path(manifest), path("raw.counts.csv")

    output:
    tuple path("${manifest.name}"), path("counts.csv")

    script:
    template "filterbyExpr.R"

}

    
workflow test {
    take:
    // Table of gene counts paired with the manifest, 
    // validated to conform to the same order of specimens in each
    counts_ch

    main:

    // Filter the counts table with filterbyExpr
    filter(counts_ch)

    }
}