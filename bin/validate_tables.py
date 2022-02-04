#!/usr/bin/python3

import os
import pandas as pd
import logging


def read_and_validate(
    fp=None,
    sep=",",
    fail_if_missing=False,
    validate_f=None,
    manifest=None
):
    """Check if a file exists. If so, read it in and validate its contents."""

    assert fp is not None, "Must provide value for `fp`"
    assert validate_f is not None, "Must provide function for `validate_f`"

    # The input file will have the '.raw' suffix
    raw_fp = fp + ".raw"

    # If the file does not exist
    if not os.path.exists(raw_fp):

        # If the missing file should trigger a failure
        if fail_if_missing:
            assert os.path.exists(raw_fp), f"ERROR: File not found {raw_fp}"

        # Otherwise, just stop the process
        return

    # Read in the file
    logging.info(f"Reading in {raw_fp}")
    df = pd.read_csv(raw_fp, sep=sep)
    logging.info(f"Read in {df.shape[0]:,} lines")

    # Validate the contents of the file
    df = validate_f(df)

    # No matter what, all of the column and index labels must be strings
    df = df.rename(
        columns=str,
        index=str
    )

    # If a manifest was provided
    if manifest is not None:
        # Align it with the order of samples in the manifest
        # and save to a file with the extension (.validated)
        align_and_save(
            df,
            manifest=manifest,
            # Output file does not have that '.raw' suffix
            fp=fp
        )

    return df


def align_and_save(df, manifest=None, fp=None):
    """Make sure that the order of samples in `df` matches the order in the manifest."""

    logging.info(f"Analyzing a table of {df.shape[0]:,} rows and {df.shape[1]:,} columns")
    logging.info(f"Ultimately saving to {fp}")

    logging.info(f"Number of rows in manifest: {manifest.shape[0]:,}")

    # Find the set of specimen names which are shared in the manifest and `df`
    shared = list(set(df.columns.values) & set(manifest.index.values))
    logging.info(f"Number of shared specimen labels: {len(shared):,}")

    assert len(shared) > 0, "No specimen labels are shared with the manifest"

    msg = "All of the values in the manifest must have values in `df`"
    assert len(shared) == manifest.shape[0], msg

    # Adjust the order of columns in `df` to match the rows in `manifest`
    df = df.reindex(
        columns=manifest.index.values
    )

    logging.info(f"Saving to {fp}")
    df.to_csv(fp)


def validate_manifest(df):
    """Validate the contents of the manifest table."""

    assert df.shape[0] > 0, "Manifest must have >1 row"
    assert df.shape[1] > 1, "Manifest must have >1 column"

    # Get the name of the first column
    sample_cname = df.columns.values[0]

    logging.info(f"The first column of manifest.csv is labeled {sample_cname}")

    # The first column must only have unique values
    msg = "The first column of manifest.csv must only have unique values."
    assert df[sample_cname].shape[0] == df[sample_cname].unique().shape[0], msg

    # Set the index on the first column
    df.set_index(sample_cname, inplace=True)

    return df


def validate_salmon_merged_gene_counts(df):
    """Validate the file named "salmon.merged.gene_counts.tsv", output by nf-core/rnaseq."""

    # The first column is named "gene_id"
    assert df.columns.values[0] == "gene_id"

    # The second column is named "gene_name"
    assert df.columns.values[1] == "gene_name"

    # Drop the gene name
    df = df.drop(columns=["gene_name"])

    # The gene_id must be unique
    assert df["gene_id"].unique().shape[0] == df.shape[0], "Every gene_id must be unique"

    # Set the index on the first column
    df.set_index("gene_id", inplace=True)

    return df


if __name__ == "__main__":

    # Check for the manifest file
    manifest = read_and_validate(
        fp="manifest.csv",
        fail_if_missing=True,
        validate_f=validate_manifest
    )

    # For each of the supported filetypes, validate it if the file exists

    # Files output by nf-core/rnaseq as "salmon.merged.gene_counts.tsv"
    salmon_merged_gene_counts = read_and_validate(
        # Read the input with the .raw suffix, but write out to the file without that suffix
        fp="salmon.merged.gene_counts.tsv",
        sep="\t",
        validate_f=validate_salmon_merged_gene_counts,
        manifest=manifest
    )    