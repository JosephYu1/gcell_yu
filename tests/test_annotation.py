import pytest

import gcell
from gcell.dna.genome import Genome
from gcell.rna.gencode import Gencode


@pytest.fixture(scope="module")
def setup_paths():
    data_dir = "/home/xf2217/.gcell_data"
    gcell.update_settings(
        {
            "annotation_dir": f"{data_dir}/annotations",
            "genome_dir": f"{data_dir}/genomes",
            "cache_dir": f"{data_dir}/cache",
        }
    )


@pytest.fixture(scope="module")
def genome(setup_paths):
    return Genome("hg38")


@pytest.fixture(scope="module")
def gencode(setup_paths):
    return Gencode("hg38", version=44)


def test_genome_sequence(genome):
    expected_seq = "TGGCGCAGGGTCCGGCGGCGCCGAGGGGTGGGCGAGCCTCGGTCTCGAGCCTCTTGGCTTCCTCCGCCCGTCCCCACTCCGGTCCCGGTTTGGGCCCTGC"
    assert str(genome.genome_seq["chr1"][1000512 : 1000512 + 100]) == expected_seq


def test_chromosome_sizes(genome):
    assert genome.chrom_sizes["chr1"] == 248956422


def test_chromosome_gap_assembly(genome):
    assert genome.chrom_gap.assembly == "hg38"


def test_blacklist_end_position(genome):
    assert genome.blacklist.iloc[0].End == 45700


def test_gencode_gene_name(gencode):
    assert gencode.gtf.iloc[0].gene_name == "DDX11L1"
