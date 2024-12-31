from .data import (
    get_lddt_from_gene_name,
    get_lddt_from_uniprot_id,
    get_seq_from_gene_name,
    get_seq_from_uniprot_id,
    get_uniprot_from_gene_name,
    initialize_all,
    organism_to_uniprot,
    ucsc_to_organism,
)
from .protein import Protein

# Initialize protein data when the module is imported
initialize_all()

__all__ = [
    "get_seq_from_gene_name",
    "get_seq_from_uniprot_id",
    "get_lddt_from_gene_name",
    "get_lddt_from_uniprot_id",
    "get_uniprot_from_gene_name",
    "organism_to_uniprot",
    "ucsc_to_organism",
    "Protein",
]
