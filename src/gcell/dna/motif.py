"""
Motif module for handling transcription factor binding site (TFBS) motifs.

This module provides base classes and functions for working with TFBS motifs, including
parsing PFM files, scanning sequences for motifs, and plotting sequence logos.

Classes
-------
Motif: Represents a single TFBS motif
MotifCluster: Represents a cluster of TFBS motifs
MotifCollection: A collection of TFBS motifs

Functions
---------
pfm_conversion: Convert a PFM file to log-odds matrices
print_results: Print the results of a motif scan
prepare_scanner: Prepare a scanner for scanning motifs
"""

import numpy as np
import pandas as pd
import seqlogo
from MOODS.parsers import pfm_to_log_odds
from MOODS.scan import Scanner
from MOODS.tools import threshold_from_p_with_precision


def pfm_conversion(
    filename, lo_bg=[2.977e-01, 2.023e-01, 2.023e-01, 2.977e-01], ps=0.01
):
    """
    Convert a PFM file to log-odds matrices.

    Parameters
    ----------
    filename : str
        The path to the PFM file
    lo_bg : list, optional
        The background frequencies for the log-odds matrices, defaults to [0.2977, 0.2023, 0.2023, 0.2977]
    ps : float, optional
        The precision for the threshold, defaults to 0.01

    Returns
    -------
    bool, list
        Whether the conversion was successful and the log-odds matrices
    """
    mat = pfm_to_log_odds(str(filename), lo_bg, ps)
    if len(mat) != 4:
        return False, mat
    else:
        return True, mat


def print_results(header, seq, matrices, matrix_names, results):
    """
    Print the results of a motif scan.

    Parameters
    ----------
    header : str
        The header of the sequence
    seq : str
        The sequence to scan
    matrices : list
        The log-odds matrices to use for the scan
    matrix_names : list
        The names of the matrices
    results : list
        The results of the scan
    """
    # split results into forward and reverse strands
    fr = results[: len(matrix_names)]
    rr = results[len(matrix_names) :]

    # mix the forward and reverse strand results
    mixed_results = [
        [(r.pos, r.score, "+", ()) for r in fr[i]]
        + [(r.pos, r.score, "-", ()) for r in rr[i]]
        for i in range(len(matrix_names))
    ]

    output = []
    for matrix, matrix_name, result in zip(matrices, matrix_names, mixed_results):
        # determine the length of the matrix
        matrix_length = len(matrix[0]) if len(matrix) == 4 else len(matrix[0]) + 1

        # sort the results by position
        sorted_results = sorted(result, key=lambda r: r[0])

        # create the output for each result
        for r in sorted_results:
            strand = r[2]
            pos = r[0]
            hitseq = seq[pos : pos + matrix_length]
            output.append([header, matrix_name, pos, strand, r[1], hitseq])
    return output


def prepare_scanner(matrices_all, bg=[2.977e-01, 2.023e-01, 2.023e-01, 2.977e-01]):
    """
    Prepare scanner for scanning motif.

    Parameters
    ----------
    matrices_all : list
        The log-odds matrices to use for the scan
    bg : list, optional
        The background frequencies for the log-odds matrices, defaults to [0.2977, 0.2023, 0.2023, 0.2977]
    """
    scanner = Scanner(7)
    scanner.set_motifs(
        matrices_all,
        bg,
        [threshold_from_p_with_precision(m, bg, 0.0001, 200, 4) for m in matrices_all],
    )
    return scanner


class Motif:
    """Base class for TFBS motifs.

    Parameters
    ----------
    id
        The ID of the motif
    gene_name
        The name of the gene
    dbd
        The DNA binding domain of the motif
    database
        The database the motif is from
    cluster_id
        The ID of the cluster the motif is in
    cluster_name
        The name of the cluster the motif is in
    pfm
        The path to the PFM file
    """

    def __init__(
        self,
        id: str,
        gene_name: str,
        dbd: str,
        database: str,
        cluster_id: str,
        cluster_name: str,
        pfm: str,
    ):
        self.id = id
        self.gene_name = gene_name
        self.dbd = dbd
        self.database = database
        self.cluster_id = cluster_id
        self.cluster_name = cluster_name
        pfm = pd.read_csv(pfm, sep="\t", header=None).T
        pfm.columns = ["A", "C", "G", "T"]
        self.pfm = seqlogo.CompletePm(pfm=seqlogo.Pfm(pfm))

    def __repr__(self) -> str:
        return f"Motif(id={self.id}, gene_name={self.gene_name}, dbd={self.dbd}, database={self.database}, cluster_id={self.cluster_id}, cluster_name={self.cluster_name})"

    def plot_logo(
        self,
        filename=None,
        format="png",
        size="large",
        ic_scale=True,
        ic_ref=0.2,
        ylabel="",
        show_xaxis=False,
        show_yaxis=False,
        show_ends=False,
        rotate_numbers=False,
        color_scheme="classic",
        logo_title="",
        fineprint="",
    ):
        """plot seqlogo of motif using pfm file

        Parameters
        ----------
        filename : str, optional
            The path to the file to save the logo, defaults to None
        format : str, optional
            The format to save the logo, defaults to "png"
        size : str, optional
            The size of the logo, defaults to "large"
        ic_scale : bool, optional
            Whether to scale the information content, defaults to True
        ic_ref : float, optional
            The reference information content, defaults to 0.2
        ylabel : str, optional
            The label for the y-axis, defaults to ""
        show_xaxis : bool, optional
            Whether to show the x-axis, defaults to False
        show_yaxis : bool, optional
            Whether to show the y-axis, defaults to False
        show_ends : bool, optional
            Whether to show the ends, defaults to False
        rotate_numbers : bool, optional
            Whether to rotate the numbers, defaults to False
        color_scheme : str, optional
            The color scheme to use, defaults to "classic"
        logo_title : str, optional
            The title of the logo, defaults to ""
        fineprint : str, optional
            The fineprint to add to the logo, defaults to ""

        Notes
        -----
        If logo_title is "id", the cluster ID will be used as the title.

        For more detailed documentation, see http://weblogo.threeplusone.com/manual.html
        """
        if logo_title == "id":
            logo_title = self.cluster_id
        return seqlogo.seqlogo(
            self.pfm,
            filename=filename,
            format=format,
            size=size,
            ic_scale=ic_scale,
            ic_ref=ic_ref,
            ylabel=ylabel,
            show_xaxis=show_xaxis,
            show_yaxis=show_yaxis,
            show_ends=show_ends,
            rotate_numbers=rotate_numbers,
            color_scheme=color_scheme,
            logo_title=logo_title,
            fineprint=fineprint,
        )


class MotifCluster:
    """Base class for TFBS motif clusters.

    Parameters
    ----------
    id
        The ID of the motif cluster
    name
        The name of the motif cluster
    seed_motif
        The ID of the seed motif
    motifs
        The motifs in the cluster
    annotations
        The annotations for the motif cluster
    """

    def __init__(self):
        self.id = None
        self.name = None
        self.seed_motif = None
        self.motifs = MotifCollection()
        self.annotations = None

    def __repr__(self) -> str:
        return f"MotifCluster(id={self.id}, name={self.name}, seed_motif={self.seed_motif})"

    def get_gene_name_list(self):
        """Get list of gene names for TFs in the motif cluster"""
        return np.concatenate([motif.gene_name for motif in self.motifs.values()])


class MotifClusterCollection:
    """List of TFBS motif clusters.

    Properties
    ----------
    annotations : pd.DataFrame
        The annotations for the motif clusters
    """

    def __init__(self):
        super().__init__()
        self.annotations = None


class MotifCollection(dict):
    """Dictionary of TFBS motifs."""

    def __init__(self):
        super().__init__()

    def __repr__(self) -> str:
        return "\n".join(self.keys())

    def get_motif_list(self):
        """Get list of motifs."""
        return list(self.keys())

    def get_motif(self, motif_id):
        """Get motif by ID."""
        return self[motif_id]
