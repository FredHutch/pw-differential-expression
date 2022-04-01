#!/usr/bin/python3
"""Validate the contents of a counts file using an associated metadata table."""

import os
import pandas as pd
import logging


def get_sep(fp):
    """Return the separator value which should be used, based on the file extension."""
    
    # Remove the '.gz', if any
    if fp.endswith('.gz'):
        fp = fp[:-3]

    # If the extension is .csv
    if fp.endswith('.csv'):

        # The separator is ','
        return ','

    # If the extension is .tsv
    elif fp.endswith('.tsv'):

        # The separator is '\t'
        return '\t'

    else:

        msg = f"Did not recognize file extension: {fp.split('.')[-1]}"
        raise Exception(msg)


def validate_counts(
    # The path to the manifest and counts table will be filled in by Nextflow prior to execution
    manifest_csv="${manifest_table}",
    counts_input="${counts_table}",
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
    logger.info(f"Reading in {manifest_csv}")
    manifest = pd.read_csv(manifest_csv, index_col=0, sep=get_sep(manifest_csv))

    # Log the specimens defined in the manifest
    logger.info("Specimens defined in the manifest:")
    for n in manifest.index.values:
        logger.info(n)

    # Read in the counts, using the first column as the index
    logger.info(f"Reading in {counts_input}")
    counts = pd.read_csv(counts_input, index_col=0, sep=get_sep(counts_input))

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

    # Tell the user if there are any missing
    if len(missing_specimens) > 0:
        print(f"Missing specimens from counts columns: {', '.join(list(missing_specimens))}")

        # Drop those specimens from the manifest
        manifest = manifest.drop(index=list(missing_specimens))

    # Reorder the columns of the counts to match the rows of the manifest
    counts = counts.reindex(
        columns=manifest.index.values
    )

    # Write out the counts file to a CSV
    logger.info(f"Writing out {counts_output}")
    counts.to_csv(counts_output)

    # Write out the manifest file to a CSV
    manifest_output = "${manifest_table.name}"
    logger.info(f"Writing out {manifest_output}")
    manifest.to_csv(manifest_output)


if __name__ == "__main__":

    # Validate that the counts file has the expected format, and
    # write out a copy which has the column order aligned with the
    # order of specimens in the manifest
    validate_counts()
