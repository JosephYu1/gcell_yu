import logging

import pandas as pd
import requests

logger = logging.getLogger(__name__)


class InterProAPI:
    """A class to interact with InterPro API for protein domain and feature data."""

    def __init__(self):
        """Initialize InterPro API client."""
        self.base_url = "https://www.ebi.ac.uk/interpro/api"

    def _make_request(self, endpoint: str, params: dict | None = None) -> dict:
        """Make a request to the InterPro API.

        Parameters
        ----------
        endpoint: API endpoint to query
        params: Optional query parameters

        Returns
        -------
        JSON response from the API

        Raises
        ------
        requests.exceptions.RequestException: If the API request fails
        """
        url = f"{self.base_url}/{endpoint}"
        try:
            response = requests.get(url, params=params)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            logger.error(f"API request failed: {e}")
            raise

    def get_protein_entries(self, uniprot_id: str) -> pd.DataFrame:
        """Get all InterPro entries for a protein.

        Parameters
        ----------
        uniprot_id: UniProt accession ID

        Returns
        -------
            DataFrame containing entry information including:
            - Entry ID
            - Entry type
            - Entry name
            - Source database
            - Location on protein sequence
        """
        endpoint = f"entry/interpro/protein/uniprot/{uniprot_id}"

        try:
            data = self._make_request(endpoint)
            entries = []
            for result in data.get("results", []):
                entry_info = {
                    "entry_id": result.get("metadata", {}).get("accession"),
                    "entry_name": result.get("metadata", {}).get("name"),
                    "entry_type": result.get("metadata", {}).get("type"),
                    "source_database": result.get("metadata", {}).get(
                        "source_database"
                    ),
                    "description": result.get("metadata", {}).get("description"),
                }

                # Extract locations if available
                locations = result.get("proteins", [])[0].get(
                    "entry_protein_locations", []
                )
                if locations:
                    for loc in locations:
                        entry_loc = entry_info.copy()
                        entry_loc.update(
                            {
                                "start": loc.get("fragments", [{}])[0].get("start"),
                                "end": loc.get("fragments", [{}])[0].get("end"),
                                "score": loc.get("score"),
                            }
                        )
                        entries.append(entry_loc)
                else:
                    entries.append(entry_info)

            return pd.DataFrame(entries)

        except Exception as e:
            logger.error(f"Failed to get protein entries for {uniprot_id}: {e}")
            return pd.DataFrame()

    def get_entry_info(self, entry_id: str) -> dict:
        """Get detailed information about an InterPro entry.

        Parameters
        ----------
        entry_id: InterPro entry ID (e.g., IPR000001)

        Returns
        -------
        Dictionary containing entry information
        """
        endpoint = f"entry/interpro/{entry_id}"
        return self._make_request(endpoint)

    def search_by_name(self, query: str, entry_type: str | None = None) -> pd.DataFrame:
        """Search InterPro entries by name/description.

        Parameters
        ----------
        query: Search term
        entry_type: Optional filter by entry type (e.g., "Domain", "Family", "Repeat")

        Returns
        -------
        DataFrame containing matching entries
        """
        endpoint = "entry/interpro"
        params = {"search": query}
        if entry_type:
            params["type"] = entry_type

        try:
            data = self._make_request(endpoint, params)
            entries = []
            for result in data.get("results", []):
                entry_info = {
                    "entry_id": result.get("metadata", {}).get("accession"),
                    "entry_name": result.get("metadata", {}).get("name"),
                    "entry_type": result.get("metadata", {}).get("type"),
                    "description": result.get("metadata", {}).get("description"),
                }
                entries.append(entry_info)
            return pd.DataFrame(entries)
        except Exception as e:
            logger.error(f"Search failed for query '{query}': {e}")
            return pd.DataFrame()

    def get_domain_architecture(self, uniprot_id: str) -> pd.DataFrame:
        """Get the domain architecture of a protein.

        Parameters
        ----------
        uniprot_id: UniProt accession ID

        Returns
        -------
        DataFrame containing ordered domain information
        """
        df = self.get_protein_entries(uniprot_id)
        if not df.empty:
            # Filter for domain entries and sort by position
            domains = df[df["entry_type"] == "domain"].sort_values("start")
            return domains
        return pd.DataFrame()

    def get_protein_families(self, uniprot_id: str) -> pd.DataFrame:
        """Get protein family classifications.

        Parameters
        ----------
        uniprot_id: UniProt accession ID

        Returns
        -------
            DataFrame containing family information
        """
        df = self.get_protein_entries(uniprot_id)
        if not df.empty:
            # Filter for family entries
            families = df[df["entry_type"] == "family"]
            return families
        return pd.DataFrame()

    def get_protein_go_terms(self, uniprot_id: str) -> list[dict]:
        """Get GO terms associated with protein's InterPro entries.

        Parameters
        ----------
        uniprot_id: UniProt accession ID

        Returns
        -------
            List of dictionaries containing GO term information
        """
        endpoint = f"protein/uniprot/{uniprot_id}/"
        try:
            data = self._make_request(endpoint)
            data = data.get("metadata", {}).get("go_terms", [])
            data = pd.DataFrame(data)
            data["category"] = data["category"].apply(lambda x: x["name"])
            return data.sort_values("category")
        except Exception as e:
            logger.error(f"Failed to get GO terms for {uniprot_id}: {e}")
            return []
