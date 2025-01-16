# gcell

[![Documentation Status](https://readthedocs.org/projects/gcell/badge/?version=latest)](https://gcell.readthedocs.io/en/latest/)

`gcell` is a comprehensive toolkit for genomic data analysis, focusing on cell type-specific regulatory analysis, DNA sequence manipulation, protein structure prediction, and pathway analysis. It integrates various modules to facilitate the study of different aspects of gene expression regulation.

# News
- I dropped `graphviz` and `pygraphviz` dependency due to complexity of maintain installing across different platform. It only affects the network layout for drawing causal graph. `nx.spring_layout(G)` is used instead of `nx.nx_agraph.graphviz_layout(G)` now, which is uglier unforturnately.
- Feature: Now you can load pre-infered cell types on getdemo website easily:
```python
from gcell.cell.celltype import GETDemoLoader
g = GETDemoLoader()
print(g.available_celltypes) # this gives you a list of cell type names
g.load_celltype('Plasma Cell')
```
- Fix: `zarr` has been limited to `<3.0.0` to avoid s3 problem


## Goal
The long term goal of this package is to create a open-source, community-involved python-centric playground/tool-set for future AI Agent to discover new biology through predictive model.


## Key Modules

- **Celltype**: Tools for cell type analysis, including gene expression and motif analysis. Basis for `get_model` interpretation analysis.
- **DNA**: Functions for DNA sequence manipulation, motif scanning, and track visualization.
- **RNA**: Classes for handling GENCODE gene annotations and GTF files.
- **Protein**: Functionality for protein domain analysis (Uniprot, InterPro) and AlphaFold2 predictions parsing as well as retrieve protein-protein interaction networks from the STRING database..
- **Pathway**: Tools for pathway (GO, KEGG, Reactome, etc.) analysis using gprofiler.


## Installation

```
pip install git+https://github.com/GET-Foundation/gcell.git@main
```

## License

`gcell` is open-source software licensed under the MIT License. See the LICENSE file for more details.
