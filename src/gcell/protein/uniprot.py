from pathlib import Path

import pandas as pd
import requests
from Bio.Seq import Seq


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
