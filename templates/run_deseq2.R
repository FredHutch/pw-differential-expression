#!/usr/bin/env Rscript

library("DESeq2")
library("BiocParallel")
register(MulticoreParam(${task.cpus}))


# Get the names of the files to process
manifest_fp = "${manifest}"
counts_fp = "${counts}"

# Read in the manifest and counts table
manifest = read.table(manifest_fp, header=TRUE, sep=",", row.names=1)
counts = read.table(counts_fp, header=TRUE, sep=",", row.names=1)
cnames = names(counts)

# Make sure that all counts are integers
counts = data.frame(
    lapply(counts,as.integer),
    row.names = rownames(counts)
)
names(counts) = cnames

# Split up the manifest filename, which has the format
# "{comp_column}.[continuous|categorical].manifest.csv"
manifest_fields = strsplit(manifest_fp, split = "[.]")[[1]]
stopifnot(length(manifest_fields) == 4)

# The first field in the manifest filename is the column indicating the
# metadata field to use in the formula
test_col = manifest_fields[1]

# Any additional grouping columns will be provided with the Nextflow parameter `group_cols`
group_cols = strsplit("${params.group_cols}", split = ",")[[1]]

# If >=1 grouping columns were provided
if ( group_cols[1] != "" ){

    # Make formula which uses the test column in the last position
    design = formula(paste("~", paste(group_cols, collapse = " + "), "+", test_col))

} else {
    # Otherwise, the formula will just have the test column
    design = formula(paste("~", test_col))
}

# Make the DEseq2 dataset
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = manifest,
    design = design
)

# Run the analysis
dds <- DESeq(dds)

# Get the results
res <- results(dds)

# Format as a DataFrame
res_df <- as.data.frame(res)

# Add the FDR-adjusted q-value
res_df\$qvalue <- p.adjust(res_df\$pvalue, method="${params.fdr_method}")

# Write out the results
write.csv(
    res_df, 
    file=paste(test_col, "DEseq2.csv", sep="."),
    quote=FALSE
)
