# Protein Module

```{eval-rst}
.. module:: gcell.protein
```

The Protein module provides functionality for working with protein sequences, structures, and AlphaFold2 predictions.

```{eval-rst}
.. currentmodule:: gcell.protein
```

## Core Functions

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: ../generated/

   get_seq_from_gene_name
   get_seq_from_uniprot_id
   get_lddt_from_gene_name
   get_lddt_from_uniprot_id
   get_uniprot_from_gene_name
```

### Protein Class

Base class for working with protein sequences and AlphaFold2 predictions.

Key methods:
- `plot_plddt()`: Plots pLDDT scores with optional domain annotations
- `plotly_plddt()`: Interactive Plotly version of pLDDT plot
- `get_domain_from_uniprot()`: Gets domain information from UniProt

```{eval-rst}
.. currentmodule:: gcell.protein.uniprot
```

## UniProt API

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: ../generated/

   UniProtAPI
```

### UniProtAPI

Class for interacting with the UniProt REST API.

Key methods:
- `get_uniprot_id()`: Gets UniProt ID from gene name
- `get_protein_sequence()`: Gets protein sequence from UniProt ID
- `get_domains()`: Gets domain information for a protein
- `get_protein_info()`: Gets detailed protein information

```{eval-rst}
.. currentmodule:: gcell.protein.af2
```

## AlphaFold2 Analysis

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: ../generated/

   AFResult
   AFMonomer
   AFHomodimer
   AFMultimer
   AFPairseg
```

### AFResult

Base class for analyzing AlphaFold2 prediction results.

Key attributes:
- `plddt`: Per-residue pLDDT confidence scores
- `pae`: Predicted Aligned Error matrix
- `iptm`: Interface TM-score
- `ptm`: Predicted TM-score

### AFPairseg

Class for analyzing AlphaFold2 predictions of protein segment pairs.

Key methods:
- `plot_plddt_gene1()`: Plots pLDDT scores for first protein
- `plot_plddt_gene2()`: Plots pLDDT scores for second protein
- `plot_score_heatmap()`: Plots heatmap of interaction scores

## Data Management

```{eval-rst}
.. currentmodule:: gcell.protein.data
```

### Organism Mappings

```python
# Map organism names to UniProt identifiers
organism_to_uniprot = {
    "human": "HUMAN_9606",
    "mouse": "MOUSE_10090",
    "rat": "RAT_10116",
    # ...
}

# Map UCSC genome builds to organisms
ucsc_to_organism = {
    "hg38": "human",
    "mm10": "mouse",
    "hg19": "human",
    # ...
}
```

### Data Loading

The module implements lazy loading of protein data:
- UniProt sequences
- Gene name to UniProt ID mappings
- AlphaFold pLDDT scores
- UniProt XML schema

Data is downloaded automatically when first needed.
