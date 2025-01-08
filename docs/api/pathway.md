```{eval-rst}
.. module:: gcell.ontology.pathway
```

# Pathway Module API Documentation

```{eval-rst}
.. currentmodule:: gcell.ontology.pathway
```

This document provides an overview of the `pathway.py` module, which includes classes and functions for pathway analysis using gprofiler.

## Classes

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: ../generated/

   Pathways
```

### Pathways

- **Description**: Class to represent pathways loaded from a gmt file.

## Functions

### read_gmt_file

- **Description**: Read gmt file and return a dataframe with pathway id as index and gene set as value.
- **Parameters**:
  - `gmt_file`: Path to the gmt file.
- **Returns**: DataFrame with pathway id as index and gene set as value.

### get_tf_pathway

- **Description**: Get pathway for a pair of transcription factors.
- **Parameters**:
  - `tf1`: First transcription factor.
  - `tf2`: Second transcription factor.
  - `cell`: Cell object containing gene annotations.
  - `filter_str`: Filter string for pathway analysis.
- **Returns**: Tuple containing gene sets and filtered pathway data.

### plot_geneset

- **Description**: Plot gene set expression.
- **Parameters**:
  - `intersect_genes`: List of intersect genes.
  - `ng_ball_exp`: DataFrame of gene expression.
  - `geneset`: List of genes in the gene set.
  - `geneset_name`: Name of the gene set.
  - `sample_category_name`: Name of the sample category.
- **Returns**: A plot of gene set expression.

### fisher_exact_test

- **Description**: Perform Fisher's exact test.
- **Parameters**:
  - `set1`: List of genes in the first set.
  - `set2`: List of genes in the second set.
  - `background`: List of genes in the background.
- **Returns**: Tuple containing the p-value, fold enrichment, and odds ratio.

### hypergeometric_test

- **Description**: Perform hypergeometric test.
- **Parameters**:
  - `set1`: List of genes in the first set.
  - `set2`: List of genes in the second set.
  - `background`: List of genes in the background.
- **Returns**: Tuple containing the p-value and fold enrichment.

### plot_fold_enrichment_with_significance

- **Description**: Plot fold enrichment with significance.
- **Parameters**:
  - `gene_lists`: List of gene lists.
  - `gene_list_names`: List of gene list names.
  - `df_genes`: DataFrame of genes.
  - `gene_annot`: DataFrame of gene annotations.
  - `ax`: Axes to plot on.
- **Returns**: Axes with the plot.

## Usage

The `pathway.py` module is designed to be used in genomic data analysis pipelines where pathway analysis is required. It provides a comprehensive set of tools for analyzing pathways and gene sets.

For more detailed usage examples, refer to the tutorials and examples provided in the documentation.
