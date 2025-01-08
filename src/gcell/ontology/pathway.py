"""
Module for pathway analysis using gprofiler.

Classes
-------
Pathways: Class to represent pathways loaded from a gmt file.

Functions
---------
read_gmt_file: Read gmt file and return a dataframe with pathway id as index and gene set as value.
get_tf_pathway: Get pathway for a pair of transcription factors.
plot_geneset: Plot gene set expression.
fisher_exact_test: Perform Fisher's exact test.
hypergeometric_test: Perform hypergeometric test.
plot_fold_enrichment_with_significance: Plot fold enrichment with significance.
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from gprofiler import GProfiler
from matplotlib import pyplot as plt
from scipy.stats import fisher_exact, hypergeom

from ..rna.gencode import Gencode
from ..rna.gene import GeneSets

gp = GProfiler(return_dataframe=True)


def read_gmt_file(gmt_file) -> pd.DataFrame:
    """
    Read gmt file and return a dataframe with pathway id as index and gene set as value.

    Parameters
    ----------
    gmt_file: Path to the gmt file

    Returns
    -------
    DataFrame: A dataframe with pathway id as index and gene set as value
    """
    pathway_dict = {}
    with Path(gmt_file).open() as f:
        for line in f:
            line = line.strip().split("\t")
            pathway_dict[line[0]] = (line[0].split(":")[0], line[1], line[2:])

    pathway_df = pd.DataFrame.from_dict(
        pathway_dict, orient="index", columns=["source", "description", "genes"]
    )
    return pathway_df


class Pathways:
    """A class to represent pathways loaded from a gmt file."""

    def __init__(
        self,
        gmt_file=None,
        assembly=None,
        gencode_version=None,
        annotation_dir=None,
        config=None,
    ):
        """
        Initialize the Pathways class.

        Parameters
        ----------
        gmt_file : str or Path
            Path to the gmt file.
        assembly : str, optional
            The genome assembly version.
        gencode_version : str, optional
            The Gencode version.
        annotation_dir : str or Path, optional
            Directory containing the annotation files.
        config : Config, optional
            A configuration object. If provided, it overrides the other parameters.
        """
        if config:
            assembly = config.get("assembly")
            gencode_version = config.get("gencode_version")
            annotation_dir = Path(config.get("annotation_dir"))
            gmt_file = config.get("pathway_file")
            gmt_file = annotation_dir / gmt_file
            if not gmt_file.exists():
                raise FileNotFoundError(f"Pathway file not found at {gmt_file}")

        self.gencode = Gencode(
            assembly=assembly, version=gencode_version, gtf_dir=annotation_dir
        )
        self.pathways = read_gmt_file(gmt_file)

    def get_gene_sets(self, genes) -> GeneSets:
        """Get GeneSets object from the list of genes."""
        return self.gencode.get_genes(genes)

    def __repr__(self) -> str:
        """Return source of the pathways and number of pathways for each source."""
        return self.pathways.groupby("source").size().to_string()

    def query_pathway(self, query_str="cell cycle") -> pd.DataFrame:
        """Query a pathway using a string."""
        return self.pathways.query("description.str.contains(@query_str)")


def get_tf_pathway(tf1, tf2, cell, filter_str="term_size<1000 & term_size>100"):
    df = cell.gene_annot.query("pred>0")
    tf1_genes = cell.gene_annot.iloc[
        cell.gene_by_motif.data[tf1].sort_values().tail(10000).index.values
    ]
    tf2_genes = cell.gene_annot.iloc[
        cell.gene_by_motif.data[tf2].sort_values().tail(10000).index.values
    ]
    intersect_genes = tf1_genes.merge(tf2_genes, on="gene_name").gene_name.unique()
    tf1_genes = tf1_genes.gene_name.unique()
    tf2_genes = tf2_genes.gene_name.unique()
    background = df.query("pred>0").gene_name.unique()
    # keep only specific genes by remove intersect genes
    tf1_genes = np.setdiff1d(tf1_genes, intersect_genes)
    tf2_genes = np.setdiff1d(tf2_genes, intersect_genes)
    go_tf1 = gp.profile(
        organism="hsapiens",
        query=list(tf1_genes),
        user_threshold=0.05,
        no_evidences=False,
        background=list(background),
    )
    go_tf2 = gp.profile(
        organism="hsapiens",
        query=list(tf2_genes),
        user_threshold=0.05,
        no_evidences=False,
        background=list(background),
    )
    go_intersect = gp.profile(
        organism="hsapiens",
        query=list(intersect_genes),
        user_threshold=0.05,
        no_evidences=False,
        background=list(background),
    )
    go_tf1_filtered = go_tf1.query(filter_str)
    go_tf2_filtered = go_tf2.query(filter_str)
    go_intersect_filtered = go_intersect.query(filter_str)

    return (
        tf1_genes,
        tf2_genes,
        intersect_genes,
        go_tf1_filtered,
        go_tf2_filtered,
        go_intersect_filtered,
    )


def plot_geneset(
    intersect_genes, ng_ball_exp, geneset, geneset_name, sample_category_name="G183S"
) -> plt.Figure:
    """
    Plot gene set expression.

    Parameters
    ----------
    intersect_genes: List of intersect genes
    ng_ball_exp: DataFrame of gene expression
    geneset: List of genes in the gene set
    geneset_name: Name of the gene set
    sample_category_name: Name of the sample category

    Returns
    -------
    Figure: A plot of gene set expression
    """
    data = ng_ball_exp.loc[
        [g for g in intersect_genes if g in ng_ball_exp.index and g in geneset]
    ]
    data = data.loc[data.std(1) > 0.1]
    sample_name = data.columns
    # replace multiple '_' to one
    sample_name = [re.sub("_+", "_", s) for s in sample_name]
    data.columns = sample_name
    sample_category = ["FAMIALL" in s.split("_")[0] for s in sample_name]
    # sample_category_color use Tab20
    sample_category_color = sns.color_palette("tab20", len(set(sample_category)))
    sample_category_color = dict(zip(set(sample_category), sample_category_color))
    sample_category_color = [sample_category_color[s] for s in sample_category]
    sns.clustermap(
        data,
        cmap="RdBu_r",
        figsize=(7, 10),
        vmax=5,
        vmin=-5,
        row_cluster=True,
        col_cluster=True,
        z_score=0,
        col_colors=sample_category_color,
        cbar_kws={"label": "Z-score"},
        xticklabels=True,
        yticklabels=True,
        method="ward",
    )
    from matplotlib.patches import Patch

    legend_elements = [
        Patch(
            facecolor=sns.color_palette("tab20", len(set(sample_category)))[i],
            edgecolor="k",
            label=s,
        )
        for i, s in enumerate(set(sample_category))
    ]
    plt.legend(
        handles=legend_elements,
        title=sample_category_name,
        bbox_to_anchor=(-0.5, 1),
        loc="upper right",
    )
    # set title
    plt.title(f"{geneset_name}")
    return plt


def fisher_exact_test(set1, set2, background) -> tuple[float, float, float]:
    """
    Perform Fisher's exact test.

    Parameters
    ----------
    set1: List of genes in the first set
    set2: List of genes in the second set
    background: List of genes in the background

    Returns
    -------
    tuple: A tuple containing the p-value, fold enrichment, and odds ratio
    """
    contingency_table = [
        [len(set1.intersection(set2)), len(set1.difference(set2))],
        [len(set2.difference(set1)), len(background.difference(set1.union(set2)))],
    ]

    # Perform Fisher's exact test
    odds_ratio, p_value = fisher_exact(contingency_table, alternative="greater")
    # return p_value and common genes
    return p_value, len(set1.intersection(set2)) / len(set2), odds_ratio


def hypergeometric_test(set1, set2, background):
    """
    Perform hypergeometric test.

    Parameters
    ----------
    set1: List of genes in the first set
    set2: List of genes in the second set
    background: List of genes in the background

    Returns
    -------
    tuple: A tuple containing the p-value and fold enrichment
    """
    set1 = set(set1)
    set2 = set(set2)
    background = set(background)
    p_value = hypergeom.sf(
        len(set1.intersection(set2)), len(background), len(set1), len(set2)
    )
    fold_enrichment = (
        len(set1.intersection(set2)) / len(set1) / len(set2) * len(background)
    )
    return p_value, fold_enrichment


def plot_fold_enrichment_with_significance(
    gene_lists, gene_list_names, df_genes, gene_annot, ax=None
) -> plt.Axes:
    """
    Plot fold enrichment with significance.

    Parameters
    ----------
    gene_lists: List of gene lists
    gene_list_names: List of gene list names
    df_genes: DataFrame of genes
    gene_annot: DataFrame of gene annotations
    ax: Axes to plot on

    Returns
    -------
    Axes: Axes with the plot
    """
    p_values = {}
    ratios = {}
    for gene_list, gene_list_name in zip(gene_lists, gene_list_names):
        p_values[gene_list_name], ratios[gene_list_name] = hypergeometric_test(
            gene_list, df_genes, gene_annot.query("pred>0").gene_name.unique()
        )
    if ax is None:
        fig, ax = plt.subplots(figsize=(1.5, 2))
    # use set_colors
    sns.barplot(
        x=gene_list_names,
        y=[ratios[gene_list_name] for gene_list_name in gene_list_names],
        ax=ax,
    )
    ax.set_ylabel("Fold enrichment")
    ax.set_ylim(0, 1.5)
    # rotate x ticks
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    # add star to significant group
    for gene_list_name in gene_list_names:
        if p_values[gene_list_name] < 0.05 / len(gene_list_names):
            ax.text(
                gene_list_names.index(gene_list_name),
                ratios[gene_list_name],
                "*",
                ha="center",
                va="center",
            )
    return ax
