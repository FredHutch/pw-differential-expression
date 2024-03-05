#!/usr/bin/env python3

import numpy as np
import os
import pandas as pd
import logging

# Set up logging
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [collect_all.py] %(message)s'
)
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Write to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)

# Keep a list of all of the data that we've found
dat = []


def neg_log10_pvalue(df: pd.DataFrame) -> pd.Series:
    """
    Calculate the negative log10 of the pvalue
    """
    # Find the smallest non-negative p-value
    min_p = df['pvalue'][df['pvalue'] > 0].min()
    return -df['pvalue'].clip(lower=min_p).apply(np.log10)


for fp in os.listdir("."):

    if fp.endswith(".csv"):

        logging.info(f"Parsing {fp}")

        # Parse the details of the analysis from the file name
        variable, method = fp[:-len(".csv")].split(".", 1)

        logging.info(f"Reading in {method} results for {variable} from {fp}")

        dat.append(
            pd.read_csv(
                fp
            ).rename(
                columns={
                    "PValue": "pvalue",
                    "P.Value": "pvalue",
                    "log2FoldChange": "logFC",
                    "QValue": "qvalue",
                    "adj.P.Val": "qvalue"
                }
            ).assign(
                method=method,
                variable=variable,
                neg_log10_pvalue=neg_log10_pvalue
            )
        )

# Concatenate all results
df = pd.concat(dat)

# Remove any rows with missing values
df = df.dropna()

# Fix the empty header for the row names
# also make sure that all outputs have the same names for
# `pvalue` and `logFC`
df = df.rename(
    columns={
        "Unnamed: 0": "gene_id"
    }
)

# Write out to CSV
df.to_csv("DE_results.csv", index=None)
