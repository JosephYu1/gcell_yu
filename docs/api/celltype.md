```{eval-rst}
.. module:: gcell.cell.celltype
```

# Celltype Module API Documentation

```{eval-rst}
.. currentmodule:: gcell.cell.celltype
```

This document provides an overview of the `celltype.py` module, which is used for cell type analysis in genomic data. The module includes several classes and methods for analyzing cell type-specific genomic data, including gene expression, motif analysis, and causal relationships between genes.

## Classes

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: ../generated/

   gcell.cell.celltype.Celltype
   gcell.cell.celltype.GETCellType
   gcell.cell.celltype.GETHydraCellType
```

### Celltype

- **Description**: Base class for cell type analysis, providing core functionality for analyzing cell type-specific genomic data.
- **Parameters**:
  - `features`: Array of feature names/identifiers.
  - `num_region_per_sample`: Number of regions per sample.
  - `celltype`: Name/identifier of the cell type.
  - `data_dir`: Directory containing data files.
  - `interpret_dir`: Directory for interpretation results.
  - `assets_dir`: Directory for assets/resources.
  - `input`: Whether to load input data.
  - `jacob`: Whether to load Jacobian data.
  - `embed`: Whether to load embedding data.
  - `num_cls`: Number of classes.
  - `s3_file_sys`: S3 filesystem object for remote storage.

### GETCellType

- **Description**: Extended cell type class with additional functionality, used for backwards compatibility for demo website.
- **Parameters**:
  - `celltype`: Cell type name.
  - `config`: Configuration object.
  - `s3_file_sys`: S3 file system object.

### GETHydraCellType

- **Description**: Cell type class optimized for hydra model analysis, including zarr-based data storage and efficient processing of large datasets.
- **Parameters**:
  - `celltype`: Name/identifier of the cell type.
  - `zarr_path`: Path to zarr data.
  - `prediction_target`: Prediction target.

## Methods

### Celltype Methods

- `load_gene_annot()`: Load gene annotations from feather file.
- `get_gene_idx(gene_name)`: Get the index of a gene in the gene list.
- `get_tss_idx(gene_name)`: Get the TSS index in the peak annotation.
- `get_gene_jacobian(gene_name, multiply_input)`: Get the jacobian of a gene.
- `get_input_data(peak_id, focus, start, end)`: Get input data from self.input_all using a slice.
- `get_tss_jacobian(jacob, tss, multiply_input)`: Get the jacobian of a TSS.
- `gene_jacobian_summary(gene, axis, multiply_input, stats)`: Summarizes the Jacobian for a given gene.

### GETCellType Methods

- Inherits all methods from `Celltype`.

### GETHydraCellType Methods

- Inherits all methods from `Celltype`.
- `get_gene_by_motif(overwrite)`: Calculate gene by motif data for hydra-based cell-type.

## Usage

This module is designed to be used in genomic data analysis pipelines where cell type-specific analysis is required. It provides a comprehensive set of tools for analyzing gene expression, motif relationships, and causal interactions in genomic data.

For more detailed usage examples, refer to the tutorials and examples provided in the documentation.
