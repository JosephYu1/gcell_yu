"""
UniProt API Client
=================

A module providing a comprehensive interface to the UniProt API for protein-related data retrieval
and XML downloads.

Classes
-------
UniProtAPI
    A class to interact with UniProt API for protein-related data retrieval and batch operations.

Example
-------
.. code-block:: python

    from gcell.protein.uniprot import UniProtAPI

    # Initialize the client
    client = UniProtAPI()

    # Get protein sequence
    sequence = client.get_protein_sequence("P12345")

    # Download XML files for multiple genes
    client.batch_download_xml(input_csv="genes.csv", output_dir="xml_files", max_workers=8)
"""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import requests
from Bio.Seq import Seq

from .._settings import get_setting


def get_default_xml_dir() -> Path:
    """Get the default directory for UniProt XML files from settings."""
    return Path(get_setting("cache_dir")) / "uniprot_xml"


class UniProtAPI:
    """A class to interact with UniProt API for protein-related data retrieval."""

    def __init__(self):
        """Initialize UniProt API client."""
        self.base_url = "https://rest.uniprot.org/uniprotkb"
        self.search_url = f"{self.base_url}/search"
        self.fetch_url = f"{self.base_url}/stream"

    def get_uniprot_id(
        self, gene_name: str, organism: str = "Homo sapiens", reviewed: bool = True
    ) -> list[str]:
        """
        Get UniProt ID(s) from gene name.

        Args:
            gene_name: Gene name to search for
            organism: Organism name (default: "Homo sapiens")
            reviewed: Whether to return only reviewed (SwissProt) entries (default: True)

        Returns:
            List of UniProt IDs
        """
        query = f'gene:{gene_name} AND organism_name:"{organism}"'
        if reviewed:
            query += " AND reviewed:true"

        params = {"query": query, "format": "json"}

        response = requests.get(self.search_url, params=params)
        response.raise_for_status()

        results = response.json()
        return [item["primaryAccession"] for item in results.get("results", [])]

    def get_protein_sequence(self, uniprot_id: str) -> Seq | None:
        """
        Get protein sequence from UniProt ID.

        Args:
            uniprot_id: UniProt ID

        Returns:
            Protein sequence as Seq object or None if not found
        """
        params = {"query": f"accession:{uniprot_id}", "format": "json"}

        response = requests.get(self.search_url, params=params)
        response.raise_for_status()

        results = response.json()
        if results.get("results"):
            sequence = results["results"][0]["sequence"]["value"]
            return Seq(sequence)
        return None

    def get_domains(
        self, uniprot_id: str, xml_dir: str | Path | None = None
    ) -> pd.DataFrame:
        """
        Get domain information from UniProt.

        Args:
            uniprot_id: UniProt ID
            xml_dir: Optional directory containing cached XML files

        Returns:
            DataFrame containing domain information
        """
        from .data import _get_schema  # Import here to avoid circular imports

        schema = _get_schema()

        if xml_dir is not None:
            xml_dir = Path(xml_dir)
            xml_file = xml_dir / f"{uniprot_id}.xml"
            if not xml_file.exists():
                url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
                response = requests.get(url)
                response.raise_for_status()
                xml_file.write_bytes(response.content)
            url = str(xml_file)
        else:
            url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"

        entry_dict = schema.to_dict(url)
        features = entry_dict["entry"][0]["feature"]
        df = []

        for feature in features:
            feature_type = feature["@type"]
            if feature_type == "chain":
                continue
            if "begin" in feature["location"]:
                if "@description" in feature:
                    feature_description = feature["@description"]
                else:
                    feature_description = feature["@type"]
                feature_begin = feature["location"]["begin"]["@position"]
                feature_end = feature["location"]["end"]["@position"]
            else:
                continue

            df.append(
                {
                    "feature_type": feature_type,
                    "feature_description": feature_description,
                    "feature_begin": int(feature_begin),
                    "feature_end": int(feature_end),
                }
            )

        return pd.DataFrame(df)

    def get_protein_info(self, uniprot_id: str) -> dict | None:
        """
        Get detailed protein information from UniProt ID.

        Args:
            uniprot_id: UniProt ID

        Returns:
            Dictionary containing protein information or None if not found
        """
        params = {"query": f"accession:{uniprot_id}", "format": "json"}

        response = requests.get(self.search_url, params=params)
        response.raise_for_status()

        results = response.json()
        if results.get("results"):
            return results["results"][0]
        return None

    def download_database(
        self,
        output_dir: str | Path,
        organism: str = "Homo sapiens",
        reviewed: bool = True,
    ) -> Path:
        """
        Download UniProt database for a specific organism.

        Args:
            output_dir: Directory to save the database
            organism: Organism name (default: "Homo sapiens")
            reviewed: Whether to download only reviewed (SwissProt) entries (default: True)

        Returns:
            Path to the downloaded file
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        query = f'organism_name:"{organism}"'
        if reviewed:
            query += " AND reviewed:true"

        params = {
            "query": query,
            "format": "tsv",
            "fields": "accession,id,gene_names,protein_name,organism_name,length,sequence",
        }

        output_file = output_dir / f"uniprot_{organism.replace(' ', '_')}.tsv"

        response = requests.get(self.fetch_url, params=params, stream=True)
        response.raise_for_status()

        output_file.write_bytes(response.content)

        return output_file

    def search_proteins(self, query: str, fields: list[str] = None) -> pd.DataFrame:
        """
        Search proteins using custom query and return results as DataFrame.

        Args:
            query: Search query in UniProt syntax
            fields: List of fields to retrieve (default: basic fields)

        Returns:
            DataFrame containing search results
        """
        if fields is None:
            fields = [
                "accession",
                "id",
                "gene_names",
                "protein_name",
                "organism_name",
                "length",
                "sequence",
            ]

        params = {"query": query, "format": "tsv", "fields": ",".join(fields)}

        response = requests.get(self.fetch_url, params=params)
        response.raise_for_status()

        # Create DataFrame from TSV response
        df = pd.read_csv(pd.StringIO(response.text), sep="\t")
        return df

    def download_xml(
        self, uniprot_id: str, output_dir: str | Path | None = None
    ) -> bool:
        """
        Download UniProt XML file for a given UniProt ID.

        Parameters
        ----------
        uniprot_id : str
            UniProt identifier
        output_dir : str or Path or None, optional
            Directory where XML files will be saved (default: cache_dir/uniprot_xml)

        Returns
        -------
        bool
            True if download was successful, False otherwise
        """
        output_dir = get_default_xml_dir() if output_dir is None else Path(output_dir)

        output_dir.mkdir(parents=True, exist_ok=True)

        try:
            response = requests.get(f"https://www.uniprot.org/uniprot/{uniprot_id}.xml")
            response.raise_for_status()
            output_file = output_dir / f"{uniprot_id}.xml"
            output_file.write_bytes(response.content)
            print(f"Successfully downloaded {uniprot_id}.xml")
            return True
        except Exception as e:
            print(f"Failed to download {uniprot_id}.xml: {str(e)}")
            return False

    def batch_download_xml(
        self,
        input_csv: str | Path,
        output_dir: str | Path | None = None,
        max_workers: int = 16,
    ) -> None:
        """
        Download UniProt XML files in batch for a list of gene names.

        Parameters
        ----------
        input_csv : str or Path
            Path to input CSV file containing gene names (one per line)
        output_dir : str or Path or None, optional
            Directory where XML files will be saved (default: cache_dir/uniprot_xml)
        max_workers : int, optional
            Number of parallel download workers (default: 16)

        Notes
        -----
        The input CSV file should contain one gene name per line without a header.
        The function will automatically map gene names to UniProt IDs and download
        the corresponding XML files.

        Example
        -------
        .. code-block:: python

            client = UniProtAPI()
            client.batch_download_xml(
                input_csv="genes.csv", output_dir="xml_files", max_workers=8
            )
        """
        output_dir = get_default_xml_dir() if output_dir is None else Path(output_dir)

        from .data import get_uniprot_from_gene_name

        # Load target list from CSV file
        target_list = pd.read_csv(input_csv, header=None).iloc[:, 0].tolist()

        # Get gene name to UniProt ID mapping
        genename_to_uniprot = get_uniprot_from_gene_name()

        def process_target(gene_name: str) -> None:
            """Process a single target gene and download its UniProt XML file."""
            if gene_name in genename_to_uniprot:
                uniprot_id = genename_to_uniprot[gene_name]
                self.download_xml(uniprot_id, output_dir)
            else:
                print(f"UniProt ID for {gene_name} not found")

        # Download XML files in parallel
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(process_target, gene_name) for gene_name in target_list
            ]
            for future in futures:
                future.result()
