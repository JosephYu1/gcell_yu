```{eval-rst}
.. module:: gcell.dna
```

# DNA Module API Documentation

```{eval-rst}
.. currentmodule:: gcell.dna
```

The DNA module provides functionality for working with DNA sequences, motifs, mutations, and genomic regions.

## Classes

```{eval-rst}
.. currentmodule:: gcell.dna.sequence
```

### DNASequence

- **Description**: Base class for representing and manipulating DNA sequences. Inherits from Seq in Biopython.
- **Key Methods**:
  - `get_reverse_complement()`: Returns the reverse complement sequence.
  - `padding(left=0, right=0, target_length=0)`: Pads sequence with N's.
  - `mutate(pos, alt)`: Creates new sequence with mutation at specified position.
  - `one_hot`: Property that returns one-hot encoded sequence.

### DNASequenceCollection

- **Description**: Collection of DNASequence objects with methods for batch operations.
- **Key Methods**:
  - `from_fasta(filename)`: Creates collection from FASTA file.
  - `mutate(pos_list, alt_list)`: Applies mutations to sequences.
  - `scan_motif(motifs)`: Scans sequences for motif matches.
  - `save_zarr()`: Saves sequences in Zarr format.

```{eval-rst}
.. currentmodule:: gcell.dna.motif
```

## Motif Analysis

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: ../generated/

   Motif
   MotifCluster
```

### Motif

- **Description**: Base class for transcription factor binding site motifs.
- **Key Methods**:
  - `plot_logo()`: Plots sequence logo for motif.
  - `__repr__()`: String representation of motif.

### MotifCluster

- **Description**: Collection of related motifs.
- **Key Methods**:
  - `get_gene_name_list()`: Gets associated gene names.

```{eval-rst}
.. currentmodule:: gcell.dna.track
```

## Track Visualization

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: ../generated/

   Track
```

### Track

- **Description**: Class for visualizing genomic tracks.
- **Key Methods**:
  - `plot_tracks()`: Plots tracks with optional gene annotations.
  - `plot_tracks_with_motif_density()`: Plots tracks with motif density.
  - `generate_bedgraph()`: Exports track as bedGraph.
  - `generate_bigwig()`: Exports track as BigWig.
