"""
Module for handling genomic tracks.

Classes
-------
Track: A class to represent and manipulate genomic track data.

Functions
---------
get_iqr_threshold: Get the IQR threshold for a specific track for peak calling.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
from matplotlib import gridspec

from .genome import ChromSize


def get_iqr_threshold(
    track, min_length=100, merge_length=100, return_indices=True, iqr_factor=1.2
):
    """
    Get the IQR threshold for a specific track for peak calling.

    Parameters
    ----------
    track : np.ndarray
        The track data.
    min_length : int, optional
        The minimum length of a peak. Defaults to 100.
    merge_length : int, optional
        The maximum distance to merge consecutive peaks. Defaults to 100.
    return_indices : bool, optional
        Whether to return the indices of the peaks. Defaults to True.
    """
    Q1, Q3 = np.percentile(track, [25, 75])
    IQR = Q3 - Q1
    threshold = Q3 + iqr_factor * IQR
    if not return_indices:
        return threshold
    idx = np.where(track > threshold)[0]
    starts = []
    ends = []
    from itertools import groupby
    from operator import itemgetter

    for k, g in groupby(enumerate(idx), lambda ix: ix[0] - ix[1]):
        g = list(map(itemgetter(1), g))
        starts.append(g[0])
        ends.append(g[-1])
    # merge consecutive indices if they are within 100bp

    for i in range(len(starts) - 1, -1, -1):
        if starts[i] - ends[i - 1] < merge_length:
            ends[i - 1] = ends[i]
            del starts[i]
            del ends[i]
    for i in range(len(starts) - 1, -1, -1):
        if ends[i] - starts[i] < min_length:
            del starts[i]
            del ends[i]
    return np.array(list(zip(starts, ends)))


class Track:
    """
    A class to represent and manipulate genomic track data.

    This class handles genomic track data, including normalization, convolution,
    and visualization. It supports multiple tracks and provides methods for
    plotting, generating bedGraph and BigWig files.
    """

    def __init__(
        self, chrom, start, end, assembly, tracks, normalize_factor=None, conv_size=50
    ):
        """
        Initialize a Track object.

        Parameters
        ----------
        chrom : str
            The chromosome name.
        start : int
            The start position of the track.
        end : int
            The end position of the track.
        assembly : str
            The genome assembly name.
        tracks : dict
            A dictionary of track data.
        normalize_factor : dict, optional
            A dictionary of normalization factors for each track. Defaults to None.
            conv_size (int, optional): The size of the convolution window. Defaults to 50.
        """
        self.chrom = chrom
        self.start = start
        self.end = end
        self.assembly = assembly
        self.tracks = tracks
        self.ids = list(self.tracks.keys())
        self.tracks_agg = np.sum(list(self.tracks.values()), axis=0)
        self.normalize_factor = normalize_factor
        if isinstance(self.normalize_factor, dict):
            self.normalize_factor_agg = np.sum(
                list(self.normalize_factor.values()), axis=0
            )
        elif isinstance(self.normalize_factor, pd.Series):
            self.normalize_factor_agg = self.normalize_factor.sum()
        if self.normalize_factor is not None:
            for label, y in self.tracks.items():
                self.tracks[label] = y / self.normalize_factor[label] * 1e8
            self.tracks_agg = self.tracks_agg / self.normalize_factor_agg * 1e8
        # set conv_size
        self._conv_size = conv_size
        self.convoluted_tracks = {
            label: np.convolve(y, np.ones(self.conv_size) / self.conv_size, mode="same")
            for label, y in self.tracks.items()
        }
        self.convoluted_tracks_agg = np.convolve(
            self.tracks_agg, np.ones(self.conv_size) / self.conv_size, mode="same"
        )
        self._sanity_check()

    def __repr__(self):
        """
        Return a string representation of the Track object.

        Returns
        -------
        str
            A string representation of the Track object.
        """
        return (
            f"Track(chrom='{self.chrom}', start={self.start}, end={self.end}, "
            f"labels={list(self.tracks.keys())}, conv_size={self.conv_size})"
        )

    def _sanity_check(self):
        """
        Perform a sanity check on the track data.

        Raises
        ------
        ValueError: If the length of any track doesn't match the expected length.
        """
        expected_length = self.end - self.start
        for label, y in self.tracks.items():
            if len(y) != expected_length:
                raise ValueError(
                    f"Length of y for {label} ({len(y)}) does not match end-start ({expected_length})"
                )

    @property
    def conv_size(self):
        """
        Get the convolution window size.

        Returns
        -------
        int
            The convolution window size.
        """
        return self._conv_size

    @conv_size.setter
    def conv_size(self, size):
        """
        Set a new convolution size and recalculate convolved tracks.

        Parameters
        ----------
        size : int
            The new convolution window size.
        """
        self._conv_size = size
        self.convoluted_tracks = {
            label: np.convolve(y, np.ones(self.conv_size) / self.conv_size, mode="same")
            for label, y in self.tracks.items()
        }

    def plot_tracks(
        self,
        color_dict=None,
        gene_annot=None,
        genes_to_highlight=None,
        beds=None,
        out_file=None,
        center=False,
    ):
        """
        Plot the convolved tracks and optionally add gene annotations.

        Parameters
        ----------
        color_dict : dict, optional
            A dictionary mapping track labels to colors. Defaults to None.
        gene_annot : Gencode, optional
            A GTF object containing gene annotations. Defaults to None.
        genes_to_highlight : list, optional
            A list of gene names to highlight. Defaults to None.
        beds : list, optional
            A list of BED objects containing regions to plot. Defaults to None.
        out_file : str, optional
            The path to save the output plot. Defaults to None.
        """
        num_tracks = len(self.convoluted_tracks)
        additional_tracks = 0
        num_subplots = num_tracks
        if isinstance(beds, np.ndarray):
            beds = [beds]
        if beds is not None:
            num_subplots += len(beds)
            additional_tracks += len(beds)

        if gene_annot is not None and (
            self.end - self.start < 1e6 or genes_to_highlight is not None
        ):
            num_subplots += 1
            additional_tracks += 1

        fig = plt.figure(figsize=(8, 1 * num_tracks + 0.2 * additional_tracks))
        gs = gridspec.GridSpec(
            num_subplots, 1, height_ratios=[1] * num_tracks + [0.2] * additional_tracks
        )

        # Adjust the space between categories
        gs.update(hspace=0.5)

        if color_dict is None:
            colors = plt.cm.tab10.colors
            color_dict = {
                label: colors[i % 10]
                for i, label in enumerate(self.convoluted_tracks.keys())
            }

        # Plot convoluted tracks
        max_y = max([max(y) for y in self.convoluted_tracks.values()])
        if center:
            mean_y = np.stack([y for y in self.convoluted_tracks.values()]).mean(axis=0)
            if len(mean_y) > 1e5:
                mean_y = mean_y[::100]
        for i, (label, y_conv) in enumerate(self.convoluted_tracks.items()):
            axes = []
            ax = fig.add_subplot(
                gs[i],
                sharex=axes[0] if axes else None,
                sharey=axes[0] if axes else None,
            )
            x = np.arange(self.start, self.end)
            # get mean per 100bp for faster plotting if the length is too long
            if len(x) > 1e5:
                x = x[::100]
                y_conv = y_conv[::100]

            if center:
                y_conv = y_conv - mean_y

            ax.fill_between(x, y_conv, color=color_dict[label])
            ax.set_ylabel(
                f"{label}",
                rotation=0,
                color=color_dict[label],
                ha="right",
                fontsize=12,
                va="center",
            )
            ax.set_xlabel("")
            ax.set_xticks([])
            if center:
                ax.set_ylim(bottom=-max_y, top=max_y)
            else:
                ax.set_ylim(bottom=0, top=max_y)
            ax.set_xlim(self.start, self.end)
            # ax.set_yticks([])  # no y ticks
            axes.append(ax)

        y_max_lim = fig.get_axes()[-1].get_ylim()[1]

        # Plot BED regions
        gs.update(hspace=0.2)
        if beds is not None:
            for bed in beds:
                if isinstance(bed, pd.DataFrame):
                    bed = bed.query(
                        "Chromosome==@self.chrom & Start>@self.start & End<@self.end"
                    )[["Start", "End"]].values
                filtered_bed = bed[(bed[:, 0] > self.start) & (bed[:, 1] < self.end)]
                ax_bed = fig.add_subplot(gs[-additional_tracks], sharex=axes[0])
                for start, end in filtered_bed:
                    ax_bed.fill_between(
                        [start, end], 0, y_max_lim, color="lightgray", alpha=1
                    )

                ax_bed.set_yticks([])
                ax_bed.set_ylabel("")
                ax_bed.set_xlabel("")
                ax_bed.set_xticks([])
                # ax_bed.spines['top'].set_visible(False)
                # ax_bed.spines['right'].set_visible(False)
                # ax_bed.spines['bottom'].set_visible(False)
                # ax_bed.spines['left'].set_visible(False)
                additional_tracks -= 1

        if gene_annot is not None and (
            self.end - self.start < 1e6 or genes_to_highlight is not None
        ):
            gene_color_dict = {
                gene: colors[(i + 1) % 10]
                for i, gene in enumerate(gene_annot.gtf.gene_name.unique())
            }
            gene_positions = {}
            ax_gene_pos = fig.add_subplot(gs[-1], sharex=axes[0])
            for gene, df in gene_annot.query_region(
                self.chrom, self.start, self.end
            ).groupby("gene_name"):
                if genes_to_highlight is None or gene in genes_to_highlight:
                    for i, row in df.iterrows():
                        ax_gene_pos.plot(
                            row["Start"],
                            1,
                            marker="^",
                            markersize=5,
                            color=gene_color_dict[gene],
                        )
                        gene_positions[gene] = df["Start"].mean()

            # Sort genes by their positions
            sorted_genes = sorted(gene_positions.items(), key=lambda x: x[1])

            # Calculate jittered positions
            jitter_range = (self.end - self.start) * 0.01  # 1% of the total range
            for i, (gene, position) in enumerate(sorted_genes):
                jitter = (
                    (i / (len(sorted_genes) - 1) - 0.5) * jitter_range
                    if len(sorted_genes) > 1
                    else 0
                )
                jittered_position = position + jitter
                ax_gene_pos.text(
                    jittered_position,
                    0.8,
                    gene,
                    verticalalignment="top",
                    horizontalalignment="center",
                    color=gene_color_dict[gene],
                    rotation=-45,
                    ha="left",
                    fontsize=12,
                )

            ax_gene_pos.set_yticks([])
            ax_gene_pos.set_ylabel("")
            ax_gene_pos.set_xlabel("")
            ax_gene_pos.set_xticks([])
            ax_gene_pos.set_ylim(0.5, 1)
            ax_gene_pos.spines["top"].set_visible(False)
            ax_gene_pos.spines["right"].set_visible(False)
            ax_gene_pos.spines["bottom"].set_visible(False)
            ax_gene_pos.spines["left"].set_visible(False)

        fig.text(
            0.5,
            1,
            f"[{self.assembly}] {self.chrom}:{self.start}-{self.end}",
            ha="center",
        )

        plt.tight_layout()
        if out_file is not None:
            plt.savefig(out_file)
        else:
            plt.show()

    def plot_tracks_with_motif_density(
        self,
        motif_zarr_path,
        cutoff_quantile=0.99,
        conv_size=100,
        motif_index=None,
        **kwargs,
    ):
        """
        Plot tracks with motif density.

        Parameters
        ----------
        motif_zarr_path : str
            The path to the motif Zarr file.
        cutoff_quantile : float, optional
            The quantile to cutoff the motif density. Defaults to 0.99.
        conv_size : int, optional
            The size of the convolution window. Defaults to 100.
        motif_index : int, optional
            The index of the motif to plot. Defaults to None.
        cutoff_quantile : float, optional
            The quantile to cutoff the motif density. Defaults to 0.99.
        conv_size : int, optional
            The size of the convolution window. Defaults to 100.
        """
        # Get motif density
        try:
            from caesar.io.zarr_io import DenseZarrIO
        except ImportError:
            "caesar is not installed, this function is not available"
            return
        m = DenseZarrIO(motif_zarr_path)
        motif_density = m.get_track(self.chrom, self.start, self.end)
        for i, quantile in enumerate(
            np.quantile(motif_density, cutoff_quantile, axis=0)
        ):
            motif_density[:, i][motif_density[:, i] < quantile] = quantile
            motif_density[:, i] = motif_density[:, i] - quantile
        if motif_index is not None:
            motif_density = motif_density[:, motif_index]
        else:
            motif_density = motif_density.sum(axis=1)

        # Convolve motif density
        motif_density = np.convolve(
            motif_density, np.ones((conv_size,)) / conv_size, mode="same"
        )

        # Add motif density to tracks
        self.convoluted_tracks["motif_density"] = (
            motif_density / motif_density.max() * self.convoluted_tracks_agg.max()
        )

        # Plot tracks with motif density
        self.plot_tracks(**kwargs)

    def get_local_iqr_peak(
        self, label, min_length=100, merge_length=100, iqr_factor=1.2
    ):
        """
        Get local IQR peaks for a specific track.

        Parameters
        ----------
        label : str
            The label of the track to get the local IQR peaks for, or 'agg' for the aggregated track.
        min_length : int, optional
            The minimum length of a peak. Defaults to 100.
        merge_length : int, optional
            The maximum distance to merge consecutive peaks. Defaults to 100.
        return_indices : bool, optional
            Whether to return the indices of the peaks. Defaults to True.

        Returns
        -------
        np.ndarray: An array of local IQR peaks.
        """
        if label == "agg":
            y = self.convoluted_tracks_agg
        else:
            y = self.convoluted_tracks[label]
        return (
            get_iqr_threshold(
                y,
                min_length=min_length,
                merge_length=merge_length,
                return_indices=True,
                iqr_factor=iqr_factor,
            )
            + self.start
        )

    def generate_bedgraph(self, label, out_file):
        """
        Generate a bedGraph file for a specific track.

        Parameters
        ----------
        label : str
            The label of the track to generate the bedGraph for.
        out_file : str
            The output file path.

        Raises
        ------
        ValueError: If the specified label is not found in the tracks.
        """
        if label not in self.tracks:
            raise ValueError(f"Label {label} not found in tracks")

        with Path(out_file).open("w") as f:
            for i, v in enumerate(self.convoluted_tracks[label]):
                f.write(f"{self.chrom}\t{i+self.start}\t{i+self.start+1}\t{v}\n")

    def generate_bigwig(self, label, out_file):
        """
        Generate a BigWig file for a specific track.

        Parameters
        ----------
        label : str
            The label of the track to generate the BigWig for.
        out_file : str
            The output file path.

        Raises
        ------
        ValueError: If the specified label is not found in the tracks.
        """
        if label not in self.tracks:
            raise ValueError(f"Label {label} not found in tracks")

        # Create a new BigWig file
        bw = pyBigWig.open(out_file, "w")

        # load chrom_size
        chrom_sizes = ChromSize(self.assembly, "../../data").dict

        # Add chromosome sizes (you may need to adjust this based on your assembly)
        bw.addHeader(list(chrom_sizes.items()))

        # Prepare the data
        starts = list(range(self.start, self.end))
        values = self.convoluted_tracks[label]
        starts = np.array(starts)[values > 0]
        values = np.array(values)[values > 0]
        chroms = np.array([self.chrom] * len(starts))
        # Add entries to the BigWig file
        bw.addEntries(chroms, starts, ends=starts + 1, values=values)

        # Close the file
        bw.close()
