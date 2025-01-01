from gcell.protein import Protein, get_seq_from_gene_name, get_uniprot_from_gene_name


def test_protein_initialization():
    """Test if Protein class can be initialized with a gene name"""
    p = Protein("PAX5")
    assert isinstance(p, Protein)


def test_protein_plotting_methods():
    """Test if plotting methods run without errors"""
    p = Protein("PAX5")
    # Test if these methods run without raising exceptions
    p.plot_plddt_manuscript()
    p.plotly_plddt()


def test_get_seq_from_gene_name():
    """Test if correct sequence is returned for PAX5"""
    expected_sequence = (
        "MDLEKNYPTPRTSRTGHGGVNQLGGVFVNGRPLPDVVRQRIVELAHQGVRPCDISRQLRVSHGCVSKILGRYYETGSIKPGVIGGSKPKVATPKVVEKIAE"
        "YKRQNPTMFAWEIRDRLLAERVCDNDTVPSVSSINRIIRTKVQQPPNQPVPASSHSIVSTGSVTQVSSVSTDSAGSSYSISGILGITSPSADTNKRKRDEGI"
        "QESPVPNGHSLPGRDFLRKQMRGDLFTQQQLEVLDRVFERQHYSDIFTTTEPIKPEQTTEYSAMASLAGGLDDMKANLASPTPADIGSSVPGPQSYPIVTGR"
        "DLASTTLPGYPPHVPPAGQGSYSAPTLTGMVPGSEFSGSPYSHPQYSSYNDSWRFPNPGLLGSPYYYSAAARGAAPPAAATAYDRH"
    )
    assert str(get_seq_from_gene_name("PAX5")) == expected_sequence


def test_get_uniprot_from_gene_name():
    """Test if correct UniProt ID is returned for PAX5"""
    assert get_uniprot_from_gene_name("PAX5") == "Q02548"
