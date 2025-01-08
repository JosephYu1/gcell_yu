import io
from pathlib import Path

import pandas as pd
import requests

from .data import organism_to_uniprot


def get_string_network(
    identifiers: list[str],
    species: str = "human",
    network_type: str = "physical",
    required_score: int = 400,
    save_image: str | None = None,
) -> pd.DataFrame:
    """Retrieve protein-protein interaction network data from STRING database.

    This function queries the STRING database API to get protein interaction networks
    for a list of proteins. It can retrieve either physical or functional interaction
    networks for various species.

    Parameters
    ----------
    identifiers (list[str]): List of protein identifiers (gene symbols or STRING IDs).
    species (str, optional): Species name for the proteins. Defaults to "human".
    Supported species can be found in organism_to_uniprot dictionary.
    network_type (str, optional): Type of interaction network. Defaults to "physical".
        Options:
            - "physical": Direct physical protein-protein interactions
            - "functional": Both direct and indirect functional associations
    required_score (int, optional): Minimum required interaction score (0-1000).
        Defaults to 400. Higher scores indicate higher confidence:
            - Low confidence: >150
            - Medium confidence: >400
            - High confidence: >700
            - Highest confidence: >900
    save_image (str | None, optional): If provided, saves the network visualization
        to the specified file path. The file extension determines the format
        (e.g., '.png', '.pdf', '.svg'). Defaults to None.

    Returns
    -------
    pd.DataFrame: DataFrame containing the interaction network data with columns:
        - nodeId: STRING protein id of partner A
        - stringId_A: STRING protein id of partner A
        - stringId_B: STRING protein id of partner B
        - score: Combined interaction score
        Additional columns depend on the interaction type and score.

    Raises
    ------
    ValueError: If an unsupported network_type or species is provided.
    requests.RequestException: If the image download fails.

    Examples
    --------
    >>> # Get physical interactions and save network image
    >>> network = get_string_network(
    ...     ["PTCH1", "SHH", "GLI1"], save_image="hedgehog_network.png"
    ... )
    >>> # Get high-confidence functional interactions with PDF image
        >>> mouse_network = get_string_network(
        ...     identifiers=["Ptch1", "Shh", "Gli1"],
        ...     species="mouse",
        ...     network_type="functional",
        ...     required_score=700,
        ...     save_image="mouse_network.pdf",
        ... )

    Notes
    -----
    - The STRING database API is queried at https://string-db.org
    - Protein identifiers can be gene symbols or STRING IDs
    - The returned interaction scores range from 0 to 1000
    - See https://string-db.org/cgi/help for more details about score computation
    - Network images are generated using STRING's image API
    """
    if network_type not in ["physical", "functional"]:
        raise ValueError(
            f"Unsupported network type: {network_type}",
            "Supported network types: physical, functional",
        )
    identifiers = "%0d".join(identifiers)
    if species not in organism_to_uniprot:
        print("Supported species: ", organism_to_uniprot.keys())
        raise ValueError(f"Unsupported species: {species}")

    species_id = organism_to_uniprot[species].split("_")[1]

    # Get network data
    response = requests.get(
        f"https://string-db.org/api/tsv/network?identifiers={identifiers}&"
        f"species={species_id}&"
        f"network_type={network_type}&"
        f"required_score={required_score}"
    )
    network_data = pd.read_csv(io.StringIO(response.content.decode("utf-8")), sep="\t")

    # Save network image if requested
    if save_image:
        image_response = requests.get(
            f"https://string-db.org/api/image/network?"
            f"identifiers={identifiers}&"
            f"species={species_id}&"
            f"network_type={network_type}&"
            f"required_score={required_score}"
        )

        if image_response.status_code == 200:
            with Path(save_image).open("wb") as f:
                f.write(image_response.content)
        else:
            raise requests.RequestException(
                f"Failed to download network image: {image_response.status_code}"
            )

    return network_data
