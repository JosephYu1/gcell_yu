from unittest.mock import Mock, patch

import pandas as pd
import pytest
import requests

from gcell.protein.interpro import InterProAPI

# Sample API responses for mocking
SAMPLE_PROTEIN_ENTRY = {
    "results": [
        {
            "metadata": {
                "accession": "IPR000001",
                "name": "Kringle",
                "type": "Domain",
                "source_database": "PFAM",
                "description": "Kringle domain",
            },
            "proteins": {
                "locations": [
                    {"fragments": [{"start": 100, "end": 200}], "score": 0.95}
                ]
            },
        }
    ]
}

SAMPLE_GO_TERMS = {
    "results": [
        {"identifier": "GO:0004674", "name": "protein serine/threonine kinase activity"}
    ]
}

SAMPLE_PATHWAYS = {
    "results": [{"identifier": "KEGG:hsa04010", "name": "MAPK signaling pathway"}]
}

SAMPLE_SEARCH_RESULTS = {
    "results": [
        {
            "metadata": {
                "accession": "IPR000001",
                "name": "Kinase",
                "type": "Domain",
                "description": "Protein kinase domain",
            }
        }
    ]
}


@pytest.fixture
def interpro_api():
    return InterProAPI()


@pytest.fixture
def mock_response():
    """Create a mock response with success status"""
    mock = Mock()
    mock.raise_for_status = Mock()
    return mock


def test_init(interpro_api):
    """Test initialization of InterProAPI"""
    assert interpro_api.base_url == "https://www.ebi.ac.uk/interpro/api"


@patch("requests.get")
def test_get_protein_entries_success(mock_get, interpro_api, mock_response):
    """Test successful protein entries retrieval"""
    mock_response.json = Mock(return_value=SAMPLE_PROTEIN_ENTRY)
    mock_get.return_value = mock_response

    df = interpro_api.get_protein_entries("P12345")

    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert "entry_id" in df.columns
    assert df.iloc[0]["entry_id"] == "IPR000001"
    assert df.iloc[0]["start"] == 100
    assert df.iloc[0]["end"] == 200


@patch("requests.get")
def test_get_protein_entries_error(mock_get, interpro_api):
    """Test error handling in protein entries retrieval"""
    mock_get.side_effect = Exception("API Error")

    df = interpro_api.get_protein_entries("P12345")

    assert isinstance(df, pd.DataFrame)
    assert df.empty


@patch("requests.get")
def test_search_by_name_success(mock_get, interpro_api, mock_response):
    """Test successful search by name"""
    mock_response.json = Mock(return_value=SAMPLE_SEARCH_RESULTS)
    mock_get.return_value = mock_response

    df = interpro_api.search_by_name("kinase", entry_type="Domain")

    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert df.iloc[0]["entry_name"] == "Kinase"
    assert df.iloc[0]["entry_type"] == "Domain"


@patch("requests.get")
def test_get_domain_architecture_success(mock_get, interpro_api, mock_response):
    """Test successful domain architecture retrieval"""
    mock_response.json = Mock(return_value=SAMPLE_PROTEIN_ENTRY)
    mock_get.return_value = mock_response

    df = interpro_api.get_domain_architecture("P12345")

    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert df.iloc[0]["entry_type"] == "Domain"


@patch("requests.get")
def test_get_protein_families_success(mock_get, interpro_api, mock_response):
    """Test successful protein families retrieval"""
    sample_family_entry = {
        "results": [
            {
                "metadata": {
                    "accession": "IPR000001",
                    "name": "Kinase Family",
                    "type": "Family",
                    "description": "Protein kinase family",
                }
            }
        ]
    }
    mock_response.json = Mock(return_value=sample_family_entry)
    mock_get.return_value = mock_response

    df = interpro_api.get_protein_families("P12345")

    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert df.iloc[0]["entry_type"] == "Family"


@patch("requests.get")
def test_get_protein_go_terms_success(mock_get, interpro_api, mock_response):
    """Test successful GO terms retrieval"""
    mock_response.json = Mock(return_value=SAMPLE_GO_TERMS)
    mock_get.return_value = mock_response

    go_terms = interpro_api.get_protein_go_terms("P12345")

    assert isinstance(go_terms, list)
    assert len(go_terms) > 0
    assert go_terms[0]["identifier"] == "GO:0004674"


@patch("requests.get")
def test_get_protein_pathways_success(mock_get, interpro_api, mock_response):
    """Test successful pathways retrieval"""
    mock_response.json = Mock(return_value=SAMPLE_PATHWAYS)
    mock_get.return_value = mock_response

    pathways = interpro_api.get_protein_pathways("P12345")

    assert isinstance(pathways, list)
    assert len(pathways) > 0
    assert pathways[0]["identifier"] == "KEGG:hsa04010"


@patch("requests.get")
def test_get_entry_info_success(mock_get, interpro_api, mock_response):
    """Test successful entry info retrieval"""
    sample_entry_info = {
        "metadata": {"accession": "IPR000001", "name": "Kringle", "type": "Domain"}
    }
    mock_response.json = Mock(return_value=sample_entry_info)
    mock_get.return_value = mock_response

    entry_info = interpro_api.get_entry_info("IPR000001")

    assert isinstance(entry_info, dict)
    assert entry_info["metadata"]["accession"] == "IPR000001"


@patch("requests.get")
def test_api_error_handling(mock_get, interpro_api):
    """Test API error handling"""
    mock_get.side_effect = requests.exceptions.RequestException("API Error")

    with pytest.raises(requests.exceptions.RequestException):
        interpro_api._make_request("some/endpoint")


def test_empty_dataframe_handling(interpro_api):
    """Test handling of empty DataFrames"""
    empty_df = pd.DataFrame()

    # Test domain architecture with empty DataFrame
    with patch.object(interpro_api, "get_protein_entries", return_value=empty_df):
        result = interpro_api.get_domain_architecture("P12345")
        assert isinstance(result, pd.DataFrame)
        assert result.empty

        # Test protein families with empty DataFrame
        result = interpro_api.get_protein_families("P12345")
        assert isinstance(result, pd.DataFrame)
        assert result.empty
