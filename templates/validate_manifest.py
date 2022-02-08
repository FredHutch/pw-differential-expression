#!/usr/bin/env python3

import logging
import pandas as pd
import os


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


def get_params():
    """Get the values defined in the Nextflow params."""

    # Required: Column used for comparisons
    comp_col = "${params.comp_col}"
    assert comp_col != 'false', "Must specify parameter: comp_col"
    assert ' ' not in comp_col, "Comparison column name cannot contain spaces"
    
    # Reference value used for categorical comparisons
    comp_ref = "${params.comp_ref}"

    # If no value was provided, use a null value
    if comp_ref == "false":
        comp_ref = None

    # List of columns to use for batch correction
    group_cols = "${params.group_cols}".split(",")

    # If no value was specified, replace with an empty list
    if len(group_cols) == 1 and group_cols[0] == "false":
        group_cols = []

    # Make sure that grouping columns do not contain spaces
    for cname in group_cols:
        assert ' ' not in cname, f"Grouping column name cannot contain spaces ({cname})"

    # Optional filtering expression to be applied
    filter = "${params.filter}"

    # If no value was provided, use a null value
    if filter == "false":
        filter = None

    return comp_col, comp_ref, group_cols, filter

def validate_manifest(manifest="manifest.csv"):

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

    # Make sure that the input file exists
    assert os.path.exists(manifest), f"File not found: {manifest}"

    # Read in the manifest
    logger.info(f"Reading in {manifest}")
    df = pd.read_csv(manifest, sep=get_sep(manifest), index_col=0)

    # Make sure that there are rows in the manifest
    assert df.shape[0] > 1, "Manifest does not contain enough rows"

    # Get the values defined by the user in the `params` scope of the workflow
    logger.info("Parsing parameters from nextflow")
    comp_col, comp_ref, group_cols, filter = get_params()

    # FILTERING
    # If a `filter` parameter was defined
    if filter is not None:

        # Filter the manifest with that boolean expression
        logger.info(f"Applying filter: {filter}")
        df = df.query(filter)

        # Make sure that multiple rows pass the filter
        assert df.shape[0] > 1, f"Filter excludes too many rows: {filter}"

    # GROUPING
    # If grouping columns were specified
    if len(group_cols) > 0:

        logger.info(f"Validating grouping columns: {', '.join(group_cols)}")

        # Make sure that those columns are in the table
        for cname in group_cols:
            msg = f"Manifest does not contain grouping column: '{cname}'"
            assert cname in df.columns.values, msg

        # Drop any specimens which lack values for those columns
        df = df.reindex(
            index=df.reindex(
                columns=group_cols
            ).dropna(
            ).index.values
        )

        # Make sure that multiple rows have grouping values
        msg = f"Not enough specimens have valid grouping information"
        assert df.shape[0] > 1, msg

    # COMPARISON COLUMN
    # Make sure that the table contains the grouping column
    msg = f"Manifest does not contain column: {comp_col}"
    assert comp_col in df.columns.values, msg

    # If the column is all numeric
    logger.info(f"Checking for all numeric values in {comp_col}")
    try:

        df = df.assign(
            **{
                comp_col: df[comp_col].apply(float)
            }
        )

        logger.info("Values appear to all be numeric")

        # There should not be a `comp_ref` value
        msg = f"Column ({comp_col} is numeric - `comp_ref` not allowed"
        assert comp_ref is None, msg

        # Write out a table which indicates the comparison column in the file name
        fp = f"{comp_col}.continuous.manifest.csv"
        print(f"Writing out file to {fp}")
        df.to_csv(fp)

    except:

        # If the column is not all numeric
        logger.info("Values are not all numeric")

        # There must be a `comp_ref` defined
        msg = f"Column ({comp_col}) is not numeric, `comp_ref` must be defined"
        assert comp_ref is not None, msg

        # The value of `comp_ref` must be present in the `comp_col` column
        comp_ref_count = (df[comp_col] == comp_ref).sum()
        msg = f"Found value ({comp_ref}) in column ({comp_col}) {comp_ref_count} times"
        assert comp_ref_count > 0, msg

        # Iterate over each of the unique values in the `comp_col` column
        for comp_val, comp_val_count in df[comp_col].value_counts().items():

            # Skip the comparison reference value
            if comp_val == comp_ref:
                continue

            logger.info(f"Formatting a table to compare {comp_val} vs. {comp_ref}")

            # Make a DataFrame which only contains those rows where the
            # value in `comp_col` is either this value, or the `comp_ref` value

            comp_df = df.loc[
                df[comp_col].isin(
                    [comp_val, comp_ref]
                )
            ]

            # Remove the `comp_col` column, and replace it with
            # a column named for `comp_val`, containing either 0 or 1
            comp_df = comp_df.drop(
                columns=[comp_col]
            ).assign(
                **{
                    comp_val: comp_df[comp_col].apply(
                        dict(
                            comp_ref=0,
                            comp_val=1
                        ).get
                    )
                }
            )

            # Write out this table as a CSV
            fp = f"{comp_val}.categorical.manifest.csv"
            logger.info(f"Writing out file to {fp}")
            comp_df.to_csv(fp)


if __name__ == "__main__":
    validate_manifest()
