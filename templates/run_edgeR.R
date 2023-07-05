#!/usr/bin/env Rscript

# Ref: https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

library(limma)
library(edgeR)

# Get the names of the files to process
manifest_fp = "${manifest}"
counts_fp = "${counts}"

# Read in the manifest and counts table
manifest = read.table(manifest_fp, header=TRUE, sep=",", row.names=1, comment.char="")
counts = read.table(counts_fp, header=TRUE, sep=",", row.names=1, comment.char="")
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
if ( length(group_cols) > 0 ){

    # Make formula which uses the test column in the last position
    model_formula = formula(paste("~", paste(group_cols, collapse = " + "), "+", test_col))

} else {
    # Otherwise, the formula will just have the test column
    model_formula = formula(paste("~", test_col))
}

# Make a design object
design = model.matrix(model_formula, data=manifest)

# Combine the counts and the metadata
dat = DGEList(counts=counts)

# Estimate the dispersions
disp = estimateDisp(dat, design)

# Perform quasi-likelihood F-tests
fit = glmQLFit(disp, design)
qlf = glmQLFTest(fit)

# Get the output table
res_df = qlf\$table

# Add the FDR-adjusted q-value
res_df\$QValue <- p.adjust(res_df\$PValue, method="${params.fdr_method}")

# Write out the results
write.csv(
    res_df, 
    file=paste(test_col, "edgeR.csv", sep="."),
    quote=FALSE
)
