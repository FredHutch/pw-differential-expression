#!/usr/bin/python3
"""Validate the contents of a counts file using an associated metadata table."""

import os
import pandas as pd
import logging

def validate_counts(
    manifest_csv="manifest.csv",
    manifest_sep=",",
    counts_input="counts.raw",
    counts_sep="\t",
    counts_output="counts.csv"
):

    # Set up logging
    logFormatter = logging.Formatter(
        '%(asctime)s %(levelname)-8s [validate_counts] %(message)s'
    )
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Write to STDOUT
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    # Make sure that all of the expected files are present
    for fp in [manifest_csv, counts_input]:
        assert os.path.exists(fp), f"File not found: {fp}"

    # Read in the manifest, using the first column as the index
    manifest = pd.read_csv(manifest_csv, index_col=0, sep=manifest_sep)

    # Log the specimens defined in the manifest
    logger.info("Specimens defined in the manifest:")
    for n in manifest.index.values:
        logger.info(n)

    # Read in the counts, using the first column as the index
    counts = pd.read_csv(counts_input, index_col=0, sep=counts_sep)

    # Log the columns in the counts table
    logger.info("Columns in counts table:")
    for n in counts.columns.values:
        logger.info(n)

    # Make sure that every row in the manifest has a corresponding
    # column in the counts file
    
    # Get the sets of index and column values from each
    manifest_rows = set(manifest.index.values)
    counts_cols = set(counts.columns.values)

    # See if there are any rows in the manfiest which are missing
    # in the columns from the counts
    missing_specimens = manifest_rows - counts_cols

    # Raise an error if there are any missing
    msg = f"Missing specimens from counts columns: {', '.join(list(missing_specimens))}"
    assert len(missing_specimens) == 0, msg

    # Reorder the columns of the counts to match the rows of the manifest
    counts = counts.reindex(
        columns=manifest.index.values
    )

    # Write out the counts file to a CSV
    counts.to_csv(counts_output)


if __name__ == "__main__":

    # Validate that the counts file has the expected format, and
    # write out a copy which has the column order aligned with the
    # order of specimens in the manifest
    validate_counts()
