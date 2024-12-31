from pathlib import Path

import numpy as np
import pandas as pd
from pyranges import PyRanges as pr
from pyranges import read_gtf
from tqdm import tqdm

from .gene import Gene, GeneSets


class GTF:
    """Base class for handling GTF annotation files"""

    def __init__(self, gtf_path, exclude_chrs=["chrM", "chrY"]):
        """Initialize GTF reader with a path to GTF file

        Args:
            gtf_path: Path to GTF file (can be gzipped)
            exclude_chrs: List of chromosomes to exclude
        """
        self.gtf_path = Path(gtf_path)

        if not self.gtf_path.exists():
            raise FileNotFoundError(f"GTF file not found: {self.gtf_path}")

        self.gtf = self._load_gtf(exclude_chrs)

    def _load_gtf(self, exclude_chrs):
        """Load and process GTF file into standardized format"""
        gtf_df = read_gtf(str(self.gtf_path)).as_df()

        # Process transcripts for positive and negative strands
        positive = gtf_df[(gtf_df.Feature == "transcript") & (gtf_df.Strand == "+")][
            [
                "Chromosome",
                "Start",
                "Start",
                "Strand",
                "gene_name",
                "gene_id",
                "gene_type",
            ]
        ]
        negative = gtf_df[(gtf_df.Feature == "transcript") & (gtf_df.Strand == "-")][
            [
                "Chromosome",
                "End",
                "End",
                "Strand",
                "gene_name",
                "gene_id",
                "gene_type",
            ]
        ]

        columns = [
            "Chromosome",
            "Start",
            "End",
            "Strand",
            "gene_name",
            "gene_id",
            "gene_type",
        ]
        positive.columns = columns
        negative.columns = columns

        gtf_df = pd.concat([positive, negative], axis=0).drop_duplicates().reset_index()
        gtf_df["gene_id"] = gtf_df.gene_id.str.split(".").str[0]

        # Filter excluded chromosomes
        gtf_df = gtf_df[~gtf_df.Chromosome.isin(exclude_chrs)]

        # Remove genes with multiple chromosomes
        gtf_df["chrom_count"] = gtf_df.groupby("gene_name")["Chromosome"].transform(
            "nunique"
        )
        gtf_df = gtf_df[gtf_df.chrom_count == 1]

        return gtf_df

    def get_gene(self, gene_name):
        df = self.gtf[self.gtf.gene_name == gene_name]
        return Gene(
            name=df.gene_name.iloc[0],
            id=df.gene_id.iloc[0],
            chrom=df.Chromosome.iloc[0],
            strand=df.Strand.iloc[0],
            tss_list=df[
                [
                    "Chromosome",
                    "Start",
                    "End",
                    "Strand",
                    "gene_name",
                    "gene_id",
                    "gene_type",
                ]
            ],
        )

    def get_genes(self, gene_names):
        gene_names = np.intersect1d(gene_names, np.unique(self.gtf.gene_name.values))
        return GeneSets([self.get_gene(gene_name) for gene_name in tqdm(gene_names)])

    def get_gene_id(self, gene_id):
        df = self.gtf[self.gtf.gene_id.str.startswith(gene_id)]
        return Gene(
            name=df.gene_name.iloc[0],
            id=df.gene_id.iloc[0],
            chrom=df.Chromosome.iloc[0],
            strand=df.Strand.iloc[0],
            tss_list=df[
                [
                    "Chromosome",
                    "Start",
                    "End",
                    "Strand",
                    "gene_name",
                    "gene_id",
                    "gene_type",
                ]
            ],
        )

    def get_genebodies(self, gene_names=None):
        """If none, return all gene bodies in a pandas dataframe"""
        genebodies = self.gtf.query('gene_type == "protein_coding"')
        genebodies["Start"] = genebodies.groupby(["Chromosome", "gene_name"])[
            "Start"
        ].transform("min")
        genebodies["End"] = genebodies.groupby(["Chromosome", "gene_name"])[
            "End"
        ].transform("max")
        genebodies = genebodies.drop_duplicates(subset=["Chromosome", "gene_name"])
        genebodies = genebodies[
            [
                "Chromosome",
                "Start",
                "End",
                "Strand",
                "gene_name",
                "gene_id",
                "gene_type",
            ]
        ]
        if gene_names is not None:
            genebodies = genebodies[genebodies.gene_name.isin(gene_names)]
        return genebodies

    def get_exp_feather(self, peaks, extend_bp=300):
        exp = (
            pr(peaks, int64=True)
            .join(pr(self.gtf, int64=True).extend(extend_bp), how="left")
            .as_df()
        )
        return exp.reset_index(drop=True)

    def query_region(self, chrom, start, end, strand=None):
        result = self.gtf.query(
            f'Chromosome == "{chrom}" & Start > {start} & End < {end}'
        )
        if strand is not None:
            result = result.query(f'Strand == "{strand}"')
        return result
