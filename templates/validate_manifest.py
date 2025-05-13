#!/usr/bin/env python3

import logging
import pandas as pd
import os


def get_sep(fp):
    """
    Return the separator value which should be used,
    based on the file extension.
    """

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
    assert comp_col != '', "Must specify parameter: comp_col"
    assert ' ' not in comp_col, "Comparison column name cannot contain spaces"

    # Reference value used for categorical comparisons
    comp_ref = "${params.comp_ref}"

    # If no value was provided, use a null value
    if comp_ref == "":
        comp_ref = None

    # List of columns to use for batch correction
    group_cols = "${params.group_cols}".split(",")

    # If no value was specified, replace with an empty list
    if len(group_cols) == 1 and group_cols[0] == "":
        group_cols = []

    # Make sure that grouping columns do not contain spaces
    for cname in group_cols:
        msg = f"Grouping column name cannot contain spaces ({cname})"
        assert ' ' not in cname, msg

    # Optional filtering expression to be applied
    filter = "${params.filter}"

    # If no value was provided, use a null value
    if filter == "":
        filter = None

    return comp_col, comp_ref, group_cols, filter


def sample_mask_initial_numeral(s):
    """
    Any sample names which start with numerals
    will have an X prepended in the counts.
    Update the sample name to match.
    """
    if s.startswith(('1', '2', '3', '4', '5', '6', '7', '8', '9', '0')):
        print(f"Sample '{s}' will be modified to X{s}")
        return f"X{s}"
    else:
        return s


def validate_manifest(manifest="input_manifest.csv"):

    # Set up logging
    logFormatter = logging.Formatter(
        '%(asctime)s %(levelname)-8s [validate_manifest] %(message)s'
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

    # SAMPLE IDS
    # Prepend 'X' to any samples which start with numerals
    df = df.rename(
        index=sample_mask_initial_numeral
    )

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
        msg = "Not enough specimens have valid grouping information"
        assert df.shape[0] > 1, msg

    # COMPARISON COLUMN
    # Make sure that the table contains the grouping column
    msg = f"Manifest does not contain column: {comp_col}"
    assert comp_col in df.columns.values, msg

    # Write out the full set of values
    df.to_csv("manifest.csv")

    # If the column is all numeric
    logger.info(f"Checking for all numeric values in {comp_col}")
    try:

        df = df.assign(
            **{
                comp_col: df[comp_col].apply(float)
            }
        )

        logger.info("Values appear to all be numeric")

        # Using a value with spaces or periods will introduce
        # errors later on when R tries to read it in
        new_comp_col = comp_col.replace(" ", "_").replace(".", "_")
        df = df.rename(columns=dict({comp_col: new_comp_col}))
        comp_col = new_comp_col

        # There should not be a `comp_ref` value
        msg = f"Column ({comp_col} is numeric - `comp_ref` not allowed"
        assert comp_ref is None, msg

        # Write out a table which indicates the
        # comparison column in the file name
        fp = f"{comp_col}.continuous.manifest.csv"
        logger.info(f"Writing out file to {fp}")
        df.to_csv(fp)

    except Exception as e:

        # Log the exception
        logger.info(str(e))

        # If the column is not all numeric
        logger.info("Values are not all numeric")

        # There must be a `comp_ref` defined
        msg = f"Column ({comp_col}) is not numeric, `comp_ref` must be defined"
        assert comp_ref is not None, msg

        # The value of `comp_ref` must be present in the `comp_col` column
        ref_n = (df[comp_col] == comp_ref).sum()
        msg = f"Found value ({comp_ref}) in column ({comp_col}) {ref_n} times"
        assert ref_n > 0, msg

        # Make sure there are no extra spaces in the comp_col column
        df[comp_col] = df[comp_col].str.strip()

        # Iterate over each of the unique values in the `comp_col` column
        for comp_val in df[comp_col].unique():

            # Skip the comparison reference value
            if comp_val == comp_ref:
                continue

            msg = f"Formatting a table to compare {comp_val} vs. {comp_ref}"
            logger.info(msg)

            # Make a DataFrame which only contains those rows where the
            # value in `comp_col` is either this value, or the `comp_ref` value

            comp_df = df.loc[
                df[comp_col].isin(
                    [comp_val, comp_ref]
                )
            ]
            logger.info(f"Using {comp_df.shape[0]:,} / {df.shape[0]:,} samples for this comparison")

            # Using a value with spaces or periods will introduce errors later on when R tries to read it in
            comp_val_sanitized = comp_val.replace(" ", "_").replace(".", "_").replace("-", "_")

            # Remove the `comp_col` column, and replace it with
            # a column named for `comp_val_sanitized`, containing either 0 or 1
            comp_df = comp_df.drop(
                columns=[comp_col]
            ).assign(
                **{
                    comp_val_sanitized: comp_df[comp_col].apply(
                        {
                            comp_ref: 0,
                            comp_val: 1
                        }.get
                    )
                }
            )

            # Write out this table as a CSV
            fp = f"{comp_val_sanitized}.categorical.manifest.csv"
            logger.info(f"Writing out file to {fp}")
            comp_df.to_csv(fp)


if __name__ == "__main__":
    validate_manifest()
