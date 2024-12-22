from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from pyranges import PyRanges as pr
from pyranges import read_gtf
from tqdm import tqdm

from .gene import Gene, GeneSets


class Gencode:
    """Read gencode gene annotation given genome assembly and version,
    returns a pandas dataframe"""

    @classmethod
    def from_config(cls, config):
        return cls(
            assembly=config.get("assembly"),
            version=config.get("gencode_version"),
            gtf_dir=config.get("annotation_dir"),
        )

    def __init__(
        self, assembly="hg38", version=40, gtf_dir=".", exclude_chrs=["chrM", "chrY"]
    ):
        super().__init__()

        self.assembly = assembly
        self.gtf_dir = Path(gtf_dir)  # Convert to Path object

        if self.assembly == "hg38":
            self.url = f"http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{version}/gencode.v{version}.basic.annotation.gtf.gz"
            self.gtf = self.gtf_dir / f"gencode.v{str(version)}.basic.annotation.gtf.gz"
        elif self.assembly == "mm10":
            self.url = f"http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_{version}/gencode.v{version}.basic.annotation.gtf.gz"
            self.gtf = (
                self.gtf_dir
                / f"gencode.{self.assembly}.v{str(version)}.basic.annotation.gtf.gz"
            )
        elif self.assembly == "hg19":
            self.url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{version}/GRCh37_mapping/gencode.v{version}lift37.basic.annotation.gtf.gz"
            self.gtf = (
                self.gtf_dir / f"gencode.v{str(version)}lift37.basic.annotation.gtf.gz"
            )

        feather_path = self.gtf_dir / f"gencode.v{str(version)}.{self.assembly}.feather"

        if feather_path.exists():
            self.gtf = pd.read_feather(feather_path)
            self.feather_file = feather_path
        else:
            if self.gtf.exists():
                self.gtf = read_gtf(str(self.gtf)).as_df()
            else:
                # download gtf to the specified directory
                os.system(f"wget -P {self.gtf_dir} {self.url} -O {self.gtf}")
                self.gtf = read_gtf(str(self.gtf)).as_df()

            positive = self.gtf[
                (self.gtf.Feature == "transcript") & (self.gtf.Strand == "+")
            ][
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
            negative = self.gtf[
                (self.gtf.Feature == "transcript") & (self.gtf.Strand == "-")
            ][
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

            positive.columns = [
                "Chromosome",
                "Start",
                "End",
                "Strand",
                "gene_name",
                "gene_id",
                "gene_type",
            ]
            negative.columns = [
                "Chromosome",
                "Start",
                "End",
                "Strand",
                "gene_name",
                "gene_id",
                "gene_type",
            ]

            self.gtf = (
                pd.concat([positive, negative], axis=0).drop_duplicates().reset_index()
            )
            self.gtf["gene_id"] = self.gtf.gene_id.str.split(".").str[0]
            self.gtf.to_feather(feather_path)
            self.feather_file = feather_path

        self.gtf = self.gtf[~self.gtf.Chromosome.isin(exclude_chrs)]
        # count number of chromosomes per gene and remove the genes with multiple chromosomes
        self.gtf["chrom_count"] = self.gtf.groupby("gene_name")["Chromosome"].transform(
            "nunique"
        )
        self.gtf = self.gtf[self.gtf.chrom_count == 1]
        return

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
        # get genebody for each gene in gencode using the longest transcript start and end
        genebodies = self.gtf.query(
            'gene_type == "protein_coding" & Chromosome != "chrM" & Chromosome != "chrY"'
        )
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
