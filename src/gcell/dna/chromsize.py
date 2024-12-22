from __future__ import annotations

import gzip
from io import StringIO
from pathlib import Path

import pandas as pd
import requests


class ChromSize:
    def __init__(self, assembly=None, annotation_dir=None):
        self.assembly = assembly
        self.annotation_dir = Path(annotation_dir) if annotation_dir else None
        if self.assembly is None:
            raise ValueError("assembly is not specified")
        if self.annotation_dir is None:
            raise ValueError("annotation_dir is not specified")

        self.chrom_sizes = self.parse_or_download_chrom_sizes()

    def _download_chrom_sizes(self):
        url = f"http://hgdownload.soe.ucsc.edu/goldenPath/{self.assembly}/bigZips/{self.assembly}.chrom.sizes"
        response = requests.get(url)
        if response.status_code != 200:
            raise ConnectionError("Failed to download chromosome data")
        return self._parse_chrom_data(response.text)

    def _parse_chrom_data(self, data):
        chrom_sizes = {}
        lines = data.strip().split("\n")
        for line in lines:
            parts = line.split("\t")
            if len(parts) == 2:
                chrom, length = parts
                chrom_sizes[chrom] = int(length)
        return chrom_sizes

    def get_dict(self, chr_included=None):
        if chr_included is None:
            return self.chrom_sizes
        else:
            return {chr: self.chrom_sizes.get(chr, None) for chr in chr_included}

    @property
    def dict(self):
        return self.chrom_sizes

    def save_chrom_sizes(self):
        filepath = self.annotation_dir / f"{self.assembly}_chrom_sizes.txt"
        filepath.write_text(
            "\n".join(
                f"{chrom}\t{length}" for chrom, length in self.chrom_sizes.items()
            )
        )

    def parse_or_download_chrom_sizes(self):
        filepath = self.annotation_dir / f"{self.assembly}_chrom_sizes.txt"
        if filepath.exists():
            return self._parse_chrom_data(filepath.read_text())
        else:
            return self._download_chrom_sizes()

    def as_pyranges(self):
        try:
            import pyranges as pr

            cs = pd.DataFrame(
                {
                    "Chromosome": list(self.chrom_sizes.keys()),
                    "Start": 0,
                    "End": list(self.chrom_sizes.values()),
                }
            ).sort_values(by=["Chromosome", "Start", "End"])
            return pr.PyRanges(cs, int64=True)
        except ImportError:
            raise ImportError("pyranges is not installed")

    def __repr__(self) -> str:
        return (
            f"ChromSize(assembly={self.assembly}, annotation_dir={self.annotation_dir})"
        )


class ChromGap:
    def __init__(self, assembly=None, annotation_dir=None, config=None):
        if config is not None:
            self.assembly = config.get("assembly")
            self.annotation_dir = Path(config.get("annotation_dir"))
        else:
            self.assembly = assembly
            self.annotation_dir = Path(annotation_dir) if annotation_dir else None

        if self.assembly is None:
            raise ValueError("assembly is not specified")
        if self.annotation_dir is None:
            raise ValueError("annotation_dir is not specified")

        self.agp_data = self.parse_or_download_agp()

    def _download_agp(self):
        url = f"https://hgdownload.soe.ucsc.edu/goldenPath/{self.assembly}/bigZips/{self.assembly}.agp.gz"
        response = requests.get(url)
        if response.status_code != 200:
            raise ConnectionError("Failed to download AGP data")
        return gzip.decompress(response.content).decode("utf-8")

    def _parse_agp_data(self, data):
        columns = [
            "chrom",
            "start",
            "end",
            "part_number",
            "component_type",
            "component_id",
            "component_start",
            "component_end",
            "orientation",
        ]
        return pd.read_csv(StringIO(data), sep="\t", comment="#", names=columns)

    def get_telomeres(self, return_tabix=False):
        df = self.agp_data[self.agp_data["component_start"] == "telomere"]
        if return_tabix:
            return pandas_to_tabix_region(df)
        return df

    def get_heterochromatin(self, return_tabix=False):
        df = self.agp_data[
            self.agp_data["component_start"].isin(["heterochromatin", "centromere"])
        ]
        if return_tabix:
            return pandas_to_tabix_region(df)
        return df

    def save_agp_data(self):
        filepath = self.annotation_dir / f"{self.assembly}_agp.txt"
        self.agp_data.to_csv(filepath, sep="\t", index=False)

    def parse_or_download_agp(self):
        filepath = self.annotation_dir / f"{self.assembly}_agp.txt"
        if filepath.exists():
            return pd.read_csv(filepath, sep="\t")
        else:
            data = self._download_agp()
            agp_data = self._parse_agp_data(data)
            self.agp_data = agp_data
            self.save_agp_data()
            return agp_data

    def __repr__(self) -> str:
        return (
            f"ChromGap(assembly={self.assembly}, annotation_dir={self.annotation_dir})"
        )


def pandas_to_tabix_region(df, chrom_col="chrom", start_col="start", end_col="end"):
    return " ".join(
        df.apply(
            lambda x: f"{x[chrom_col]}:{x[start_col]}-{x[end_col]}", axis=1
        ).tolist()
    )
