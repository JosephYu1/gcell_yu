```{eval-rst}
.. module:: gcell.protein
```

# Protein Module API Documentation

```{eval-rst}
.. currentmodule:: gcell.protein
```

The Protein module provides functionality for working with protein sequences, structures, and AlphaFold2 predictions.

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

- **Description**: Base class for working with protein sequences and AlphaFold2 predictions.
- **Key Methods**:
  - `plot_plddt()`: Plots pLDDT scores with optional domain annotations.
  - `plotly_plddt()`: Interactive Plotly version of pLDDT plot.
  - `get_domain_from_uniprot()`: Gets domain information from UniProt.

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

- **Description**: Class for interacting with the UniProt REST API.
- **Key Methods**:
  - `get_uniprot_id()`: Gets UniProt ID from gene name.
  - `get_protein_sequence()`: Gets protein sequence from UniProt ID.
  - `get_domains()`: Gets domain information for a protein.
  - `get_protein_info()`: Gets detailed protein information.

```{eval-rst}
.. currentmodule:: gcell.protein.stringdb
```

## STRING Database API

### get_string_network

- **Description**: Retrieve protein-protein interaction network data from the STRING database.
- **Parameters**:
  - `identifiers`: List of protein identifiers (gene symbols or STRING IDs).
  - `species`: Species name for the proteins (default: "human").
  - `network_type`: Type of interaction network ("physical" or "functional").
  - `required_score`: Minimum required interaction score (0-1000).
  - `save_image`: Path to save the network visualization image.
- **Returns**: DataFrame containing the interaction network data.

```{eval-rst}
.. currentmodule:: gcell.protein.interpro
```

## InterPro API

### InterProAPI

- **Description**: A class to interact with the InterPro API for protein domain and feature data.
- **Key Methods**:
  - `get_protein_entries(uniprot_id)`: Get all InterPro entries for a protein.
  - `get_entry_info(entry_id)`: Get detailed information about an InterPro entry.
  - `search_by_name(query, entry_type)`: Search InterPro entries by name/description.
  - `get_domain_architecture(uniprot_id)`: Get the domain architecture of a protein.
  - `get_protein_families(uniprot_id)`: Get protein family classifications.
  - `get_protein_go_terms(uniprot_id)`: Get GO terms associated with protein's InterPro entries.

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

- **Description**: Base class for analyzing AlphaFold2 prediction results.
- **Key Attributes**:
  - `plddt`: Per-residue pLDDT confidence scores.
  - `pae`: Predicted Aligned Error matrix.
  - `iptm`: Interface TM-score.
  - `ptm`: Predicted TM-score.

### AFPairseg

- **Description**: Class for analyzing AlphaFold2 predictions of protein segment pairs.
- **Key Methods**:
  - `plot_plddt_gene1()`: Plots pLDDT scores for first protein.
  - `plot_plddt_gene2()`: Plots pLDDT scores for second protein.
  - `plot_score_heatmap()`: Plots heatmap of interaction scores.

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
