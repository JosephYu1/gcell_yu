"""
Non-redundant DNA Motif Collection (Version 1.0)
==============================================

This module provides functionality to work with the non-redundant transcription factor binding site (TFBS)
motif clusters defined by Vierstra et al. The motifs are sourced from
https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/.

The module includes:
    - Loading and managing motif clusters
    - Converting between different motif formats
    - Scanning sequences for motif matches
    - Mapping between different gene naming conventions

Classes
-------
NrMotifV1
    Main class for working with the non-redundant motif collection.

Dependencies
-----------
- MOODS: Required for reverse complement functionality
- pandas: For data manipulation
- pickle: For serialization
- pathlib: For file path handling

Example
-------
>>> from gcell.dna.nr_motif_v1 import NrMotifV1
>>> motifs = NrMotifV1("path/to/motif/dir")
>>> # Get a specific motif cluster
>>> cluster = motifs.get_motif_cluster_by_name("CTCF")
"""

import os
import pickle
import tarfile
from pathlib import Path

import pandas as pd
import requests
from pandas import DataFrame

from .._settings import get_setting

try:
    from MOODS.tools import reverse_complement
except ImportError:
    print(
        "MOODS not installed. Please install MOODS to use the reverse_complement function."
    )

from .motif import (
    Motif,
    MotifCluster,
    MotifClusterCollection,
    MotifCollection,
    pfm_conversion,
    prepare_scanner,
)

annotation_dir = Path(get_setting("annotation_dir"))

other_gene_mapping = {
    "ARI5B": "ARID5B",
    "HXA13": "HOXA13",
    "ATF6A": "ATF6",
    "HNF6": "ONECUT1",
    "DMRTB": "DMRTB1",
    "COE1": "EBF1",
    "EVI1": "MECOM",
    "EWSR1-FLI1": "EWSR1-FLI1",
    "ITF2": "TCF4",
    "ARNTL": "BMAL1",
    "BHE40": "BHLHE40",
    "BHLHB3": "BHLHE41",
    "BHLHB2": "BHLHE40",
    "NGN2": "NEUROG2",
    "TWST1": "TWIST1",
    "NDF2": "NEUROD2",
    "NDF1": "NEUROD1",
    "ZNF238": "ZBTB18",
    "BHA15": "BHLHA15",
    "TFE2": "TCF3",
    "ANDR": "AR",
    "HXA10": "HOXA10",
    "HXC9": "HOXC9",
    "HXA9": "HOXA9",
    "HXA1": "HOXA1",
    "HXB7": "HOXB7",
    "HXB8": "HOXB8",
    "HXB13": "HOXB13",
    "HXB4": "HOXB4",
    "RAXL1": "RAX2",
    "MIX-A": "MIXL1",
    "CART1": "ALX1",
    "HEN1": "NHLH1",
    "KAISO": "ZBTB33",
    "MYBA": "MYBL1",
    "TF65": "RELA",
    "DUX": "DUX1",
    "STF1": "NR5A1",
    "ERR2": "ESRRB",
    "COT1": "NR2F1",
    "COT2": "NR2F2",
    "THA": "THRA",
    "THB": "THRB",
    "NR1A4": "RXRA",
    "RORG": "RORC",
    "PRGR": "PGR",
    "GCR": "NR3C1",
    "ERR3": "ESRRG",
    "ERR1": "ESRRA",
    "PO3F1": "POU3F1",
    "PO3F2": "POU3F2",
    "PO5F1": "POU5F1",
    "P53": "TP53",
    "P73": "TP73",
    "P63": "TP63",
    "POU5F1P1": "POU5F1B",
    "PO2F2": "POU2F2",
    "PO2F1": "POU2F1",
    "SUH": "RBPJ",
    "PEBB": "CBFB",
    "SRBP1": "SREBF1",
    "SRBP2": "SREBF2",
    "T": "TBXT",
    "BRAC": "TBXT",
    "TF2L1": "TFCP2L1",
    "TYY1": "YY1",
    "ZKSC1": "ZKSCAN1",
    "THA11": "THAP11",
    "OZF": "ZNF146",
    "ZNF306": "ZKSCAN3",
    "Z324A": "ZNF324",
    "AP2A": "TFAP2A",
    "AP2B": "TFAP2B",
    "AP2C": "TFAP2C",
    "TF7L1": "TCF7L1",
    "TF7L2": "TCF7L2",
    "STA5B": "STAT5B",
    "STA5A": "STAT5A",
    "BC11A": "BCL11A",
    "Z354A": "ZNF354A",
    "HINFP1": "HINFP",
    "PIT1": "POU1F1",
    "HTF4": "TCF12",
    "ZNF435": "ZSCAN16",
}


def fix_gene_name(x: str) -> str:
    if x.startswith("ZN") and not x.startswith("ZNF"):
        x = x.replace("ZN", "ZNF")
    if x.startswith("ZSC") and not x.startswith("ZSCAN"):
        x = x.replace("ZSC", "ZSCAN")
    if x.startswith("NF2L"):
        x = x.replace("NF2L", "NFE2L")
    if x.startswith("PKNX1"):
        x = "PKNOX1"
    if x.startswith("NKX") and "-" not in x:
        x = x[:-1] + "-" + x[-1]
    if x.startswith("PRD") and not x.startswith("PRDM"):
        x = x.replace("PRD", "PRDM")
    if x.startswith("NFAC"):
        x = x.replace("NFAC", "NFATC")
    if x.startswith("SMCA"):
        x = x.replace("SMCA", "SMARCA")
    if x.startswith("ZBT") and not x.startswith("ZBTB"):
        x = x.replace("ZBT", "ZBTB")
    # other mappings
    for k, v in other_gene_mapping.items():
        if x.startswith(k):
            x = v
    return x


class NrMotifV1(MotifClusterCollection):
    """
    A collection of non-redundant transcription factor binding site motif clusters.

    This class manages the motif clusters defined in the non-redundant motif collection v1.0.
    It provides functionality to load, access, and work with motif clusters and individual motifs.

    Parameters
    ----------
    motif_dir : str or Path
        Directory where motif data will be stored
    base_url : str, optional
        URL to fetch motif data from, by default points to Altius Institute resources

    Notes
    -----
    The class automatically downloads required motif data if not present in the specified directory.
    """

    def __init__(
        self,
        motif_dir: str | Path = annotation_dir / "nr_motif_v1",
        base_url: str = "https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/",
    ) -> None:
        super().__init__()
        self.motif_dir = Path(motif_dir)
        self.annotations: DataFrame = self.get_motif_data(self.motif_dir, base_url)
        matrices: list = []
        matrices_rc: list = []
        for motif in self.get_motif_list():
            filename = self.motif_dir / "pfm" / f"{motif}.pfm"
            valid = False
            if filename.exists():  # let's see if it's pfm
                valid, matrix = pfm_conversion(filename)
                matrices.append(matrix)
                matrices_rc.append(reverse_complement(matrix, 4))

        self.matrices: list = matrices
        self.matrices_all: list = self.matrices + matrices_rc
        self.matrix_names: list[str] = self.get_motif_list()
        self.cluster_names: list[str] = self.get_motifcluster_list()
        self.motif_to_cluster: dict[str, str] = (
            self.annotations[["Motif", "Name"]].set_index("Motif").to_dict()["Name"]
        )
        self.cluster_gene_list: dict[str, list[str]] = (
            self.get_motifcluster_list_genes()
        )

    # facility to export the instance as a pickle and load it back
    def __getstate__(self) -> dict:
        state = self.__dict__.copy()
        return state

    def __setstate__(self, state: dict) -> None:
        self.__dict__.update(state)

    def save_to_pickle(
        self, file_path: str | Path = annotation_dir / "nr_motif_v1.pkl"
    ) -> None:
        """
        Serialize the NrMotifV1 instance to a pickle file.

        Parameters
        ----------
        file_path : Path or str, optional
            Path where the pickle file will be saved
        """
        with Path(file_path).open("wb") as f:
            pickle.dump(self.__getstate__(), f)

    @classmethod
    def download_pickle(
        cls,
        file_path: str | Path = annotation_dir / "nr_motif_v1.pkl",
        url: str = "https://zenodo.org/record/14614892/files/nr_motif_v1.pkl?download=1",
    ):
        """
        Download the pickle file from Zenodo using pooch.

        Parameters
        ----------
        file_path : Path or str, optional
            Path where the pickle file will be saved
        url : str, optional
            URL to fetch the pickle file from
        """
        # Ensure the directory exists
        file_path.parent.mkdir(parents=True, exist_ok=True)
        from .._settings import POOCH

        # Use pooch to download the file
        POOCH.path = file_path.parent
        POOCH.registry[file_path.name] = None  # No hash available
        POOCH.urls[file_path.name] = url

        # Fetch the file
        POOCH.fetch(file_path.name)

    @classmethod
    def load_from_pickle(
        cls,
        file_path: str | Path = annotation_dir / "nr_motif_v1.pkl",
        motif_dir: str | Path | None = annotation_dir / "nr_motif_v1",
    ) -> "NrMotifV1":
        """
        Load a NrMotifV1 instance from a pickle file.

        Parameters
        ----------
        file_path : Path or str, optional
            Path to the pickle file
        motif_dir : Path or str, optional
            Directory containing motif data

        Returns
        -------
        NrMotifV1
            Loaded instance of NrMotifV1

        Notes
        -----
        Creates motif directory if it doesn't exist and downloads required data.
        """
        if not Path(file_path).exists():
            cls.download_pickle(file_path)
        with Path(file_path).open("rb") as f:
            state = pd.read_pickle(f)
        instance = cls.__new__(cls)
        instance.__setstate__(state)
        if motif_dir is not None:
            instance.motif_dir = motif_dir
            if (
                not instance.motif_dir.exists()
                or len(os.listdir(instance.motif_dir)) == 0
            ):
                # create the directory
                instance.motif_dir.mkdir(parents=True, exist_ok=True)
                instance.get_motif_data(motif_dir)
        return instance

    def get_motif_data(
        self,
        motif_dir: str | Path = annotation_dir / "nr_motif_v1",
        base_url: str = "https://zenodo.org/record/14635057/files/nr_motif_v1.tar?download=1",
    ) -> DataFrame:
        """
        Download and load motif data from the specified source.

        Parameters
        ----------
        motif_dir : Path
            Directory to store downloaded motif data
        base_url : str
            URL to fetch motif data from

        Returns
        -------
        pandas.DataFrame
            Motif annotations data

        Notes
        -----
        Downloads and extracts a tar file if not present locally. Combines motif annotations from multiple sources.
        """
        pfm_dir = motif_dir / "pfm"
        if pfm_dir.exists() and pfm_dir.is_dir() and len(os.listdir(pfm_dir)) > 0:
            pass
        else:
            tar_path = annotation_dir / "nr_motif_v1.tar"

            # Download the tar file
            print("Downloading motif data...")
            response = requests.get(base_url, stream=True)
            with Path(tar_path).open("wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)

            # Extract the tar file
            print("Extracting motif data...")
            with tarfile.open(tar_path, "r") as tar:
                tar.extractall(path=annotation_dir)

        annotations_file = motif_dir / "motif_annotations.csv"
        if annotations_file.exists():
            motif_annotations = pd.read_csv(annotations_file)
        else:
            a = pd.read_excel(f"{base_url}/motif_annotations.xlsx", sheet_name=1)
            b = pd.read_excel(f"{base_url}/motif_annotations.xlsx", sheet_name=0)
            motif_annotations = pd.merge(
                a, b, left_on="Cluster_ID", right_on="Cluster_ID"
            )
            motif_annotations.to_csv(annotations_file, index=False)
        return motif_annotations

    def get_motif_list(self) -> list[str]:
        """
        Get a sorted list of all motif IDs in the collection.

        Returns
        -------
        list
            Sorted list of motif IDs
        """
        return sorted(self.annotations.Motif.unique())

    def get_motif(self, motif_id: str) -> Motif:
        """
        Get a specific motif by its ID.

        Parameters
        ----------
        motif_id : str
            ID of the motif to retrieve

        Returns
        -------
        Motif
            Motif object containing ID, associated genes, DBD, database source, and other metadata
        """
        row = self.annotations[self.annotations.Motif == motif_id].iloc[0]
        return Motif(
            row.Motif,
            row.Motif.split("_")[0].split("+"),
            row.DBD,
            row.Database,
            row.Cluster_ID,
            row.Name,
            self.motif_dir / "pfm" / f"{row.Motif}.pfm",
        )

    def get_motifcluster_list(self) -> list[str]:
        """
        Get a sorted list of all motif cluster names.

        Returns
        -------
        list
            Sorted list of motif cluster names
        """
        return sorted(self.annotations.Name.unique())

    def get_motifcluster_list_genes(self) -> dict[str, list[str]]:
        cluster_gene_list: dict[str, list[str]] = {}
        for c in self.get_motifcluster_list():
            for g in self.get_motif_cluster_by_name(c).get_gene_name_list():
                if g.endswith("mouse"):
                    g = g.replace(".mouse", "").upper()
                else:
                    if c in cluster_gene_list:
                        cluster_gene_list[c].append(g.upper())
                    else:
                        cluster_gene_list[c] = [g.upper()]
            if c in cluster_gene_list:
                cluster_gene_list[c] = list(set(cluster_gene_list[c]))
        # fix the gene names in motif_gene_list
        for k in cluster_gene_list:
            cluster_gene_list[k] = [fix_gene_name(x) for x in cluster_gene_list[k]]
        return cluster_gene_list

    def get_motif_cluster_by_name(self, mc_name: str) -> MotifCluster:
        """
        Retrieve a motif cluster by its name.

        Parameters
        ----------
        mc_name : str
            Name of the motif cluster

        Returns
        -------
        MotifCluster
            Cluster object containing all associated motifs and metadata
        """
        mc = MotifCluster()
        mc.name = mc_name
        mc.annotations = self.annotations[self.annotations.Name == mc_name]
        mc.seed_motif = self.get_motif(mc.annotations.iloc[0].Seed_motif)
        mc.id = mc.annotations.iloc[0].Cluster_ID
        mc.motifs = MotifCollection()
        for motif_id in self.annotations[
            self.annotations.Name == mc_name
        ].Motif.unique():
            mc.motifs[motif_id] = self.get_motif(motif_id)
        return mc

    def get_motif_cluster_by_id(self, mc_id: str) -> MotifCluster:
        """
        Retrieve a motif cluster by its ID.

        Parameters
        ----------
        mc_id : str
            ID of the motif cluster

        Returns
        -------
        MotifCluster
            Cluster object containing all associated motifs and metadata
        """
        mc = MotifCluster()
        mc.name = mc_id
        mc.annotations = self.annotations[self.annotations.Cluster_ID == mc_id]
        mc.seed_motif = self.get_motif(self.annotations.iloc[0].Seed_motif)
        mc.name = self.annotations.iloc[0].Name

        mc.motifs = MotifCollection()
        for motif_id in self.annotations[
            self.annotations.Cluster_ID == mc_id
        ].Motif.unique():
            mc.motifs[motif_id] = self.get_motif(motif_id)
        return mc

    @property
    def scanner(self, bg: list[float] = [2.977e-01, 2.023e-01, 2.023e-01, 2.977e-01]):
        """
        Get a MOODS scanner configured with the current motif matrices.

        Parameters
        ----------
        bg : list, optional
            Background nucleotide frequencies [A,C,G,T]

        Returns
        -------
        MOODSScanner
            Scanner object configured for motif searching

        Notes
        -----
        Uses both forward and reverse complement matrices for comprehensive scanning.
        """
        return prepare_scanner(self.matrices_all, bg)

    def __repr__(self) -> str:
        """
        Return a string representation of the NrMotifV1 instance.

        Returns
        -------
        str
            A formatted string showing key information about the motif collection
        """
        n_motifs = len(self.matrix_names)
        n_clusters = len(self.cluster_names)
        n_matrices = len(self.matrices)

        return (
            f"NrMotifV1(\n"
            f"    motifs: {n_motifs:,}\n"
            f"    clusters: {n_clusters:,}\n"
            f"    matrices: {n_matrices:,}\n"
            f"    directory: {self.motif_dir}\n"
            f")"
        )
