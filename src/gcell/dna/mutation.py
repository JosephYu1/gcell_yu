"""
Module to handle mutations in a specific cell type that alters motifs.

Classes
-------
Mutations: Class to handle mutations
SVs: Class to handle structural variations (half-baked)

Functions
---------
read_rsid_parallel: Function to read rsids in parallel
read_rsid: Function to read rsids
read_gwas_catalog: Function to read GWAS catalog
read_vcf: Function to read VCF file
prepare_gnomad_data: Function to prepare gnomad data
fetch_rsid_data: Function to fetch rsid data
tabix_query: Function to query tabix file
tabix_index: Function to index tabix file
bgzip: Function to compress file with bgzip
"""

import concurrent.futures
import random
import subprocess
import time
from pathlib import Path
from subprocess import Popen

import pandas as pd
import requests
from tqdm import tqdm

from .._logging import get_logger
from .genome import GenomicRegionCollection
from .motif import MotifCollection
from .sequence import DNASequence, DNASequenceCollection

logger = get_logger(__name__)


class Mutations(GenomicRegionCollection):
    """Class to handle mutations. Subclass of GenomicRegionCollection.

    Parameters
    ----------
    genome: Genome
    df: DataFrame containing mutations


    """

    def __init__(self, genome, df):
        super().__init__(genome, df)
        if df is not None:
            self.collect_ref_sequence(30, 30)
            self.collect_alt_sequence(30, 30)
        return

    def collect_ref_sequence(self, upstream=0, downstream=0):
        """Collect reference sequences centered at the mutation sites.
        Adds Ref_seq attribute to the class.

        Parameters
        ----------
        upstream: Upstream length
        downstream: Downstream length
        """
        self.Ref_seq = [
            s.seq
            for s in super()
            .collect_sequence(upstream=upstream, downstream=downstream)
            .sequences
        ]

    def collect_alt_sequence(self, upstream=0, downstream=0):
        """Collect alternative sequences centered at the mutation sites.
        Adds Alt_seq attribute to the class.

        Parameters
        ----------
        upstream: Upstream length
        downstream: Downstream length
        """
        if self.Ref_seq is None:
            self.collect_ref_sequence(upstream, downstream)
        n_mut = len(self.Ref_seq)
        Alt_seq = DNASequenceCollection([DNASequence(s) for s in self.Ref_seq.values])
        Alt_seq = Alt_seq.mutate([upstream] * n_mut, self.Alt.values)
        Alt_seq = [s.seq for s in Alt_seq.sequences]
        self.Alt_seq = Alt_seq

    def get_motif_diff(self, motif: MotifCollection) -> dict:
        """Get motif difference between reference and alternative sequences.

        Parameters
        ----------
        motif: Motif

        Returns
        -------
        dict: Dictionary containing motif difference between reference and alternative sequences
        """
        Alt_seq = DNASequenceCollection(
            [
                DNASequence(row.Alt_seq, row.RSID + "_" + row.Alt)
                for i, row in self.df.iterrows()
            ]
        )
        Ref_seq = DNASequenceCollection(
            [
                DNASequence(row.Ref_seq, row.RSID + "_" + row.Ref)
                for i, row in self.df.iterrows()
            ]
        )
        return {"Alt": Alt_seq.scan_motif(motif), "Ref": Ref_seq.scan_motif(motif)}


class SVs:
    """Class to handle SVs

    Note:
        This is still a work in progress
    """

    def __init__(self, bedpe_file, genome):
        self.genome = genome
        self.bedpe_file = bedpe_file
        self.bedpe = self.read_bedpe()
        return

    def read_bedpe(self):
        """Read bedpe file"""
        pd.read_csv(self.bedpe_file, sep="\t", header=None)
        return


def bgzip(filename):
    """Call bgzip to compress a file."""
    Popen(["bgzip", "-f", filename])


def tabix_index(filename, preset="gff", chrom=1, start=4, end=5, skip=0, comment="#"):
    """Call tabix to create an index for a bgzip-compressed file."""
    Popen(
        [
            "tabix",
            "-p",
            preset,
            "-s",
            chrom,
            "-b",
            start,
            "-e",
            end,
            "-S",
            skip,
            "-c",
            comment,
        ]
    )


def tabix_query(
    filename: str,
    chrom: str,
    start: int,
    end: int,
    output_file: str,
    with_header: bool = True,
) -> str:
    """
    Calls tabix to query a VCF file and saves the output to a file.

    Parameters
    ----------
    filename (str): The path to the VCF file.
    chrom (str): The chromosome.
    start (int): The start position of the query.
    end (int): The end position of the query.
    output_file (str): The file to save the output.

    Returns
    -------
    str: The path to the output file.
    """
    query = f"{chrom}:{start}-{end}"

    # Construct the tabix command
    command = ["tabix", filename, query]

    if with_header:
        command.append("-h")

    # output_file = os.path.join(base_dir, output_file)

    try:
        with Path.open(output_file, "w") as f:
            # Execute the command and redirect the output to the file
            subprocess.run(command, stdout=f, check=True)

        # compress the output file
        bgzip(output_file)
        # index the output file
        tabix_index(output_file + ".gz")

    except subprocess.CalledProcessError as e:
        logger.error(f"An error occurred while querying: {e}")
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")

    return output_file + ".gz"


def read_gwas_catalog(genome, gwas_catalog_csv_file):
    """Read GWAS catalog for mutations

    Parameters
    ----------
    gwas_catalog_csv_file: GWAS catalog file in csv format.

    Returns
    -------
    Mutations: Class to handle mutations
    """
    gwas = pd.read_csv(gwas_catalog_csv_file, sep="\t")
    chrom = gwas["CHR_ID"].astype(str).apply(lambda x: "chr" + x)
    risk_allele = gwas["STRONGEST SNP-RISK ALLELE"].apply(lambda x: x.split("-")[1])
    variants = pd.DataFrame(
        {
            "Chromosome": chrom,
            "Start": gwas["CHR_POS"] - 1,
            "End": gwas["CHR_POS"],
            "Alt": risk_allele,
            "RSID": gwas["SNPS"],
        }
    )
    variants = variants.drop_duplicates()
    grc = GenomicRegionCollection(genome, variants)
    variants["Ref"] = [s.seq for s in grc.collect_sequence().sequences]
    # filter out variants with same risk allele and reference allele
    variants = variants.query("Ref != Alt").reset_index(drop=True)
    return Mutations(genome, variants)


def read_vcf(vcf_file):
    """Read VCF file for mutations

    Parameters
    ----------
    vcf_file: VCF file

    Returns
    -------
    Mutations: Class to handle mutations
    """
    pd.read_csv(vcf_file, sep="\t", header=None)

    return


# e.g. https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.chr1.vcf.bgz
def prepare_gnomad_data(
    gnomad_path="/pmglocal/alb2281/get_data/gnomad",
    gnomad_base_url="https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.0/vcf/genomes/",
):
    """
    Download tabix index for gnomad data.

    Parameters
    ----------
    gnomad_path: Path to store gnomad data
    gnomad_base_url: Base URL for gnomad data
    """
    for chrom in list(range(1, 23)) + ["X", "Y"]:
        chrom = f"chr{chrom}"
        print(f"Downloading {chrom}...")
        r = requests.get(
            gnomad_base_url + f"gnomad.genomes.v4.0.sites.{chrom}.vcf.bgz.tbi"
        )
        with Path.open(gnomad_path + chrom + ".vcf.bgz.tbi", "wb") as f:
            f.write(r.content)

    return


def fetch_rsid_data(server, rsid, max_retries=5):
    """Fetch RSID data with retry mechanism for rate limiting.

    Parameters
    ----------
    server: Ensembl server
    rsid: RSID
    max_retries: Maximum number of retries
    """
    ext = f"/variation/human/{rsid}?"
    for i in range(max_retries):
        try:
            r = requests.get(server + ext, headers={"Content-Type": "application/json"})
            r.raise_for_status()
            decoded = pd.DataFrame(r.json()["mappings"])
            decoded["RSID"] = rsid
            return decoded
        except requests.exceptions.HTTPError as err:
            if r.status_code == 429 and i < max_retries - 1:  # Too Many Requests
                wait_time = (2**i) + random.random()
                time.sleep(wait_time)
            else:
                raise err


def read_rsid_parallel(
    genome, rsid_list, num_workers=10
) -> tuple[Mutations, list, list]:
    """Read VCF file, only support hg38

    Parameters
    ----------
    genome: Genome
    rsid_list: List of RSIDs
    num_workers: Number of workers

    Returns
    -------
    Mutations: Class to handle mutations
    """
    server = "http://rest.ensembl.org"
    df = []

    processed_rsids = []
    failed_rsids = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        future_to_rsid = {
            executor.submit(fetch_rsid_data, server, rsid): rsid
            for rsid in tqdm(rsid_list)
        }
        for future in concurrent.futures.as_completed(future_to_rsid):
            try:
                df.append(future.result())
                processed_rsids.append(future_to_rsid[future])
            except Exception as e:
                failed_rsids.append(future_to_rsid[future])
                logger.error(f"Error processing {future_to_rsid[future]}: {e}")

    if len(df) > 0:
        df = (
            pd.concat(df)
            .query('~location.str.contains("CHR")')
            .query('assembly_name=="GRCh38"')
        )
        df["Start"] = df["start"] - 1
        df["End"] = df["start"]
        df["Chromosome"] = df.seq_region_name.apply(lambda x: "chr" + x)
        df["Ref"] = df.allele_string.apply(lambda x: x.split("/")[0])
        df["Alt"] = df.allele_string.apply(lambda x: x.split("/")[1])
        return (
            Mutations(genome, df[["Chromosome", "Start", "End", "Ref", "Alt", "RSID"]]),
            processed_rsids,
            failed_rsids,
        )
    else:
        return Mutations(genome, None), processed_rsids, failed_rsids


def read_rsid(genome, rsid_file) -> Mutations:
    """Read VCF file, only support hg38

    Parameters
    ----------
    genome: Genome
    rsid_file: File containing RSIDs

    Returns
    -------
    Mutations: Class to handle mutations
    """
    import numpy as np

    rsid_list = np.loadtxt(rsid_file, dtype=str)
    server = "http://rest.ensembl.org"
    df = []
    for rsid in tqdm(rsid_list):
        ext = f"/variation/human/{rsid}?"
        r = requests.get(server + ext, headers={"Content-Type": "application/json"})
        if not r.ok:
            r.raise_for_status()
        decoded = pd.DataFrame(r.json()["mappings"])
        decoded["RSID"] = rsid
        df.append(decoded)
    df = (
        pd.concat(df)
        .query('~location.str.contains("CHR")')
        .query('assembly_name=="GRCh38"')
    )
    df["Start"] = df["start"] - 1
    df["End"] = df["start"]
    df["Chromosome"] = df.seq_region_name.apply(lambda x: "chr" + x)
    df["Ref"] = df.allele_string.apply(lambda x: x.split("/")[0])
    df["Alt"] = df.allele_string.apply(lambda x: x.split("/")[1])

    return Mutations(genome, df[["Chromosome", "Start", "End", "Ref", "Alt", "RSID"]])
