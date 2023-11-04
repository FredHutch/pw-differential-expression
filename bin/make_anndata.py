#!/usr/bin/env python3

from anndata import AnnData
import json
import logging
import numpy as np
import pandas as pd
import scanpy as sc
from vitessce import (
    VitessceConfig,
    CoordinationType as ct,
    Component as cm,
    AnnDataWrapper
)
from vitessce.data_utils.anndata import optimize_adata

# Set up logging
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [make_anndata.py] %(message)s'
)
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Write to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)


def make_anndata(
    category: str,
    res: pd.DataFrame,
    manifest: pd.DataFrame,
    counts: pd.DataFrame,
    scale_factor=1e6
) -> AnnData:

    # Make an AnnData object (scaling to CPM)
    # Transposing so that samples are on the .obs axis
    logger.info("Scaling input data")
    adata = AnnData(
        (scale_factor * counts / counts.sum()).T,
        dtype=np.float32
    )
    sc.pp.log1p(adata)

    # Run UMAP ordination
    logger.info("Running UMAP & PCA")
    sc.pp.neighbors(adata, use_rep='X', metric='euclidean')
    sc.tl.umap(adata)
    # Run PCA
    sc.tl.pca(adata)

    logger.info("Annotating AnnData object")

    # Add the sample annotation
    adata.obs[category] = manifest[category]

    # Set a floor on the pvalues to remove zero values
    res = res.assign(
        pvalue=clip_zeros(res['pvalue']),
        qvalue=clip_zeros(res['qvalue']),
    )

    # Format the results as an indexed table
    res = (
        res
        .set_index("gene_id")
        .assign(
            neg_log10_qvalue=lambda d: -d['qvalue'].apply(np.log10),
            top_significant=lambda d: top_significant(d),
            mean_abund=adata.to_df().mean()
        )
    )

    # Add the differential abundance analysis results
    for kw, val in res.items():
        adata.var[kw] = val

    # Format the volcano plot and MA plot with plotting coordinates
    adata.varm["results"] = (
        res
        .reindex(columns=["mean_abund", "logFC", "neg_log10_padj"])
        .apply(scale_values)
        .values
    )

    return adata


def clip_zeros(r: pd.Series) -> pd.Series:
    # Find the lowest non-zero pvalue
    threshold = r.loc[r > 0].dropna().min()
    return r.clip(lower=threshold).fillna(1)


def top_significant(df: pd.DataFrame, n=100) -> pd.Series:
    """Identify the most highly significant genes."""
    # Calculate a score for each gene using the log10(qvalue) and log(fold_change)
    score = df["neg_log10_qvalue"].abs() * df["logFC"].abs()
    threshold = score.sort_values().tail(n).min()
    return (score >= threshold).apply(int)


def scale_values(r: pd.Series):
    """Scale values to have a range of 1"""
    if r.max() > r.min():
        return r / (r.max() - r.min())
    else:
        return r


def write_vitessce(
    category,
    desc="Identification of genes differentially expressed between groups",
    schema_version="1.0.16"
):

    logger.info("Setting up Vitessce config")
    vc = VitessceConfig(
        schema_version=schema_version,
        name=f"Differential Expression by {category}",
        description=desc
    )

    # Configure the dataset of samples
    samples_dataset = (
        vc
        .add_dataset(name=f"Samples: {category}")
        .add_object(
            AnnDataWrapper(
                adata_url=f"{category}.samples.zarr",
                obs_embedding_paths=["obsm/X_umap", "obsm/X_pca"],
                obs_embedding_names=["UMAP", "PCA"],
                obs_set_paths=[f"obs/{category}"],
                obs_set_names=[category],
                obs_feature_matrix_path="X",
                feature_filter_path="var/top_significant",
                coordination_values=dict(
                    obsType="Sample",
                    featureType="Gene"
                )
            )
        )
    )
    genes_dataset = (
        vc
        .add_dataset(name=f"Genes: {category}")
        .add_object(
            AnnDataWrapper(
                adata_url=f"{category}.genes.zarr",
                obs_feature_matrix_path="X",
                obs_embedding_paths=["obsm/results", "obsm/results"],
                obs_embedding_names=["Volcano", "MA Plot"],
                obs_embedding_dims=[[1, 2], [0, 1]],
                obs_set_paths=["obs/top_significant"],
                obs_set_names=["Differentially Expressed"],
                coordination_values=dict(
                    obsType="Gene",
                    featureType="Sample"
                )
            )
        )
    )

    # Set up the scatterplot with the volcano plot for genes
    scatter_volcano = vc.add_view(
        cm.SCATTERPLOT,
        dataset=genes_dataset,
        mapping="Volcano"
    )

    # Set up the scatterplot with the MA plot for genes
    scatter_ma = vc.add_view(
        cm.SCATTERPLOT,
        dataset=genes_dataset,
        mapping="MA Plot"
    )

    # Set up the scatterplot with PCA on the samples
    scatter_pca = vc.add_view(
        cm.SCATTERPLOT,
        dataset=samples_dataset,
        mapping="PCA"
    )

    # Show the distribution of a cell across samples
    show_gene_abund = vc.add_view(
        cm.OBS_SET_FEATURE_VALUE_DISTRIBUTION,
        dataset=samples_dataset
    )

    # Allow the user to select a gene
    select_gene = vc.add_view(cm.FEATURE_LIST, dataset=samples_dataset)

    # Show all values in a heatmap
    heatmap = vc.add_view(cm.HEATMAP, dataset=samples_dataset)

    for elems, obs_type, feature_type in [
        ([scatter_ma, scatter_volcano], "Gene", "Sample"),
        ([scatter_pca, show_gene_abund, select_gene, heatmap], "Sample", "Gene")
    ]:
        set_radius(vc, 5, elems)
        isolate_selections(vc, elems)
        set_obs_type(vc, elems, obs_type)
        set_feature_type(vc, elems, feature_type)

    # Lay out the entire display
    vc.layout(
        (scatter_volcano | select_gene | scatter_pca) /
        (scatter_ma | show_gene_abund | heatmap)
    )

    logger.info("Writing out Vitessce config")
    write_vt_config(vc, category)


def set_radius(vc, px, elems):

    # Customize the scatterplots, point size and zoom
    obs_radius, obs_radius_mode = vc.add_coordination(
        ct.EMBEDDING_OBS_RADIUS,
        ct.EMBEDDING_OBS_RADIUS_MODE
    )
    for elem in elems:
        elem.use_coordination(obs_radius)
        elem.use_coordination(obs_radius_mode)

    obs_radius.set_value(px)
    obs_radius_mode.set_value("manual")


def isolate_selections(vc: VitessceConfig, elems):

    for coord_scope_type, val in [
        (ct.OBS_SET_SELECTION, None),
        (ct.OBS_SET_HIGHLIGHT, None),
        (ct.OBS_SET_COLOR, None),
        (ct.OBS_COLOR_ENCODING, "cellSetSelection")
    ]:
        coord_scope = vc.add_coordination(coord_scope_type).pop()
        for elem in elems:
            elem.use_coordination(coord_scope)
        if val is not None:
            coord_scope.set_value(val)


def set_elem_attr(vc: VitessceConfig, ct_type, elems, val):

    coord_scope = vc.add_coordination(ct_type).pop()
    for elem in elems:
        elem.use_coordination(coord_scope)
    coord_scope.set_value(val)


def set_obs_type(vc: VitessceConfig, elems, obs_type):
    set_elem_attr(vc, ct.OBS_TYPE, elems, obs_type)


def set_feature_type(vc: VitessceConfig, elems, feature_type):
    set_elem_attr(vc, ct.FEATURE_TYPE, elems, feature_type)


def write_vt_config(vc: VitessceConfig, category: str):
    """Write out the configuration to JSON."""

    # Write out as JSON
    vt_fp = f"{category}.vt.json"
    with open(vt_fp, "w") as handle:
        # Get the configuration
        config = vc.to_dict(base_url=".")
        json.dump(config, handle, indent=4)


def save_anndata(adata: AnnData, category: str):
    """
    Save in two orientations:
    {category}.samples.h5ad - samples are obs, genes are var
    {category}.genes.h5ad - genes are obs, samples are var
    """

    logger.info("Optimizing AnnData object for serialization")
    gene_cols = [
        "pvalue",
        "qvalue",
        "logFC",
        "neg_log10_qvalue",
        "mean_abund",
        "top_significant"
    ]

    samples_adata = optimize_adata(
        adata,
        obs_cols=[category],
        obsm_keys=["X_umap", "X_pca"],
        var_cols=gene_cols
    )
    logger.info("Writing samples to h5ad and zarr")
    samples_adata.write_h5ad(f"{category}.samples.h5ad", compression="gzip")
    samples_adata.write_zarr(f"{category}.samples.zarr")

    genes_adata = optimize_adata(
        adata.T,
        obs_cols=gene_cols,
        var_cols=[category],
        obsm_keys=["results"]
    )
    logger.info("Writing genes to h5ad and zarr")
    genes_adata.write_h5ad(f"{category}.genes.h5ad", compression="gzip")
    genes_adata.write_zarr(f"{category}.genes.zarr")


if __name__ == "__main__":

    # Read in the data
    logger.info("Reading input data")
    DE_results = pd.read_csv("DE_results.csv")
    manifest = pd.read_csv("manifest.csv", index_col=0)
    counts = pd.read_csv("counts.csv", index_col=0)

    # Process each of the DA analyses
    for category, res in DE_results.groupby("variable"):
        logger.info(f"Processing results for '{category}'")

        # Make the AnnData object
        adata = make_anndata(category, res, manifest, counts)

        # Save to H5AD
        # Note: This will save the data in two both orientations
        # {category}.samples.h5ad - samples are obs, genes are var
        # {category}.genes.h5ad - genes are obs, samples are var
        save_anndata(adata, category)

        # Write out a vitessce configuration
        write_vitessce(category)
    logger.info("Done")
