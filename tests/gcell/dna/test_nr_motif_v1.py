import traceback
from pathlib import Path

import pandas as pd
import pytest

from gcell._settings import get_setting
from gcell.dna.motif import Motif, MotifCluster
from gcell.dna.nr_motif_v1 import NrMotifV1


@pytest.fixture(scope="module")
def test_data_dir():
    return Path(get_setting("annotation_dir")) / "nr_motif_v1"


def test_init_and_basic_functionality(test_data_dir):
    """Test initialization and verify that all methods can be called"""

    # Initialize the class
    nr_motif = NrMotifV1(test_data_dir)

    # Test basic attributes
    assert isinstance(nr_motif.annotations, pd.DataFrame)
    assert len(nr_motif.annotations) > 0

    # Test all methods to ensure they run without errors
    try:
        # Test getting motif list
        motif_list = nr_motif.get_motif_list()
        assert isinstance(motif_list, list)
        assert len(motif_list) > 0

        # Test getting a specific motif
        motif = nr_motif.get_motif(motif_list[0])
        assert isinstance(motif, Motif)

        # Test getting cluster list
        cluster_list = nr_motif.get_motifcluster_list()
        assert isinstance(cluster_list, list)
        assert len(cluster_list) > 0

        # Test getting cluster by name
        cluster = nr_motif.get_motif_cluster_by_name(cluster_list[0])
        assert isinstance(cluster, MotifCluster)

        # Test getting cluster by ID
        cluster_by_id = nr_motif.get_motif_cluster_by_id("AHR")
        assert isinstance(cluster_by_id, MotifCluster)

        # Test motif to cluster mapping
        assert isinstance(nr_motif.motif_to_cluster, dict)
        assert len(nr_motif.motif_to_cluster) > 0

        # Test cluster gene list
        assert isinstance(nr_motif.cluster_gene_list, dict)

        # Test scanner property
        scanner = nr_motif.scanner
        assert scanner is not None

    except Exception as e:
        pytest.fail(f"Method execution failed: {str(e)}, {traceback.format_exc()}")


def test_repr(test_data_dir):
    """Test the string representation of NrMotifV1"""
    nr_motif = NrMotifV1(test_data_dir)
    repr_str = repr(nr_motif)

    # Check that the repr contains key information
    assert "NrMotifV1" in repr_str
    assert "motifs:" in repr_str
    assert "clusters:" in repr_str
    assert "matrices:" in repr_str
    assert str(test_data_dir) in repr_str

    # Check that the numbers make sense
    assert int(repr_str.split("motifs:")[1].split("\n")[0].strip().replace(",", "")) > 0
    assert (
        int(repr_str.split("clusters:")[1].split("\n")[0].strip().replace(",", "")) > 0
    )
    assert (
        int(repr_str.split("matrices:")[1].split("\n")[0].strip().replace(",", "")) > 0
    )
