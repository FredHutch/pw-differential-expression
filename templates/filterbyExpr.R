#!/usr/bin/env Rscript

library(limma)
library(edgeR)

# Get the name of the manifest from Nextflow
manifest_fp = "${manifest}"

# Read in the manifest and counts table
manifest = read.table(manifest_fp, header=TRUE, sep=",", row.names=1, comment.char="")
counts = read.table("raw.counts.csv", header=TRUE, sep=",", row.names=1, comment.char="")

# Split up the manifest filename, which has the format
# "validated.{comp_column}.[continuous|categorical].manifest.csv"
manifest_fields = strsplit(manifest_fp, split = "[.]")[[1]]
stopifnot(length(manifest_fields) == 5)

starting_counts = nrow(counts)

# If the comparison is continuous
if (manifest_fields[3] == "continuous"){

    # Treat the group of samples as belonging to a single group
    group = rep(c('dummy_group'), times=ncol(counts))
    names(group) = names(counts)

    # Get a vector showing which genes to keep
    keep = filterByExpr(
        counts,
        group=group,
        min.count=${params.min_count},
        min.total.count=${params.min_total_count},
        large.n=${params.large_n},
        min.prop=${params.min_prop},
    )

# Otherwise, the comparison is categorical
}else{

    # Provide the vector of group membership
    # Get a vector showing which genes to keep
    keep = filterByExpr(
        counts,
        group=manifest[[manifest_fields[2]]],
        min.count=${params.min_count},
        min.total.count=${params.min_total_count},
        large.n=${params.large_n},
        min.prop=${params.min_prop},
    )

}

# Subset the table to just those genes which survived the filter
counts = counts[keep,]

print(paste(nrow(counts), "/", starting_counts, "genes pass the filter"))

# Write to a file
write.csv(counts, "counts.csv", row.names = TRUE)