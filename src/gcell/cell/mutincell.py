"""
# TODO: This is still a work in progress
Class to handle mutations in a specific cell type that alters motifs.

Classes
-------
MutationsInCellType: Base class to handle mutations in a specific cell type
CellMutCollection: Base class to handle mutations in a collection of cell types
GETHydraCellMutCollection: Class to handle mutations in a collection of GETHydraCellType
"""

import concurrent.futures
import contextlib
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from omegaconf import DictConfig
from pyranges import PyRanges as pr

from .._logging import get_logger
from .._settings import get_setting
from ..cell.celltype import CellType, GETCellType, GETHydraCellType
from ..config.config import load_config
from ..dna.genome import Genome, GenomicRegionCollection
from ..dna.mutation import Mutations, read_rsid_parallel
from ..dna.nr_motif_v1 import NrMotifV1
from ..dna.sequence import DNASequence, DNASequenceCollection

annotation_dir = get_setting("annotation_dir")
logger = get_logger(__name__)


class MutationsInCellType:
    """Class to handle mutations in a specific cell type

    Parameters
    ----------
    genome : Genome
        The genome object.
    df : pd.DataFrame
        The mutation data.
    cell_type : CellType
        The cell type object.

    """

    def __init__(self, genome: Genome, df: pd.DataFrame, cell_type: CellType):
        self.celltype = cell_type
        df = pr(cell_type.peak_annot).join(pr(df)).df
        # keep only variants with one base change
        df = df.query("Ref.str.len()==1 & Alt.str.len()==1")
        df["upstream"] = df.Start_b - df.Start
        df["downstream"] = df.End - df.End_b
        self.mut = GenomicRegionCollection(genome, df)
        self.upstream = self.mut.df.upstream.values
        self.downstream = self.mut.df.downstream.values
        self.Alt = self.mut.df.Alt.values
        self.Ref = self.mut.df.Ref.values
        # self.get_original_input()
        # self.get_altered_input()

    def get_original_input(self, motif):
        """
        Get the original input for the mutation.
        """
        self.Ref_peak_seq = [
            s.seq for s in self.mut.collect_sequence(upstream=0, downstream=0).sequences
        ]
        Ref_peak_seq = DNASequenceCollection(
            [
                DNASequence(self.Ref_peak_seq[i], row.RSID + "_" + row.Ref)
                for i, row in self.mut.df.iterrows()
            ]
        )
        self.Ref_input = Ref_peak_seq.scan_motif(motif)

    def get_altered_input(self, motif):
        """
        Get the altered input for the mutation.
        """
        if self.Ref_peak_seq is None:
            self.get_original_input()
        Alt_peak_seq = DNASequenceCollection(
            [DNASequence(s) for s in self.Ref_peak_seq]
        )
        Alt_peak_seq = Alt_peak_seq.mutate(list(self.upstream), self.Alt)
        Alt_peak_seq = [s.seq for s in Alt_peak_seq.sequences]
        self.Alt_peak_seq = Alt_peak_seq
        Alt_peak_seq = DNASequenceCollection(
            [
                DNASequence(self.Alt_peak_seq[i], row.RSID + "_" + row.Alt)
                for i, row in self.mut.df.iterrows()
            ]
        )
        self.Alt_input = Alt_peak_seq.scan_motif(motif)


class CellMutCollection:
    """
    Class to handle mutations in a collection of cell types.

    Parameters
    ----------
    get_config_path : str
        The path to the configuration file.
    genome_path : str
        The path to the genome file.
    motif_path : str
        The path to the motif file.
    celltype_list : list[str]
        The list of cell types.
    variant_list : list[str]
        The list of variants.

    """

    def __init__(
        self,
        get_config_path,
        genome_path,
        motif_path,
        celltype_list,
        variant_list,
        variant_to_genes,
        output_dir,
        num_workers,
        *,
        debug: bool = False,
        celltype_annot_path: str | None = None,
        cell_obj_dict: dict[str, GETCellType] | None = None,
    ):
        self.output_dir = output_dir
        if celltype_annot_path is not None:
            celltype_annot = pd.read_csv(celltype_annot_path)
            self.celltype_annot_dict = celltype_annot.set_index("id").celltype.to_dict()
        else:
            self.celltype_annot_dict = {cell_id: cell_id for cell_id in celltype_list}

        self.num_workers = num_workers
        self.get_config_path = get_config_path

        self.genome = Genome("hg38", genome_path)
        self.motif = NrMotifV1.load_from_pickle(motif_path)
        self.debug = debug

        self.get_config = load_config(get_config_path)
        self.get_config.celltype.jacob = True
        self.get_config.celltype.num_cls = 2
        self.get_config.celltype.input = True
        self.get_config.celltype.embed = False
        self.get_config.celltype.assets_dir = ""
        self.get_config.s3_file_sys = ""
        self.get_config.celltype.data_dir = "/manitou/pmg/users/xf2217/pretrain_human_bingren_shendure_apr2023/fetal_adult/"
        self.get_config.celltype.interpret_dir = (
            "/pmglocal/xf2217/Interpretation_all_hg38_allembed_v4_natac"
        )

        if self.debug:
            variant_list = variant_list[:2]
            celltype_list = celltype_list[:2]

        if not Path(self.output_dir).exists():
            Path(self.output_dir).mkdir(parents=True)
            Path(f"{self.output_dir}/csv").mkdir(parents=True)
            Path(f"{self.output_dir}/feather").mkdir(parents=True)
            Path(f"{self.output_dir}/logs").mkdir(parents=True)

        self.celltype_list = celltype_list
        self.celltype_cache = {}
        if cell_obj_dict is not None:
            for cell_id, cell_obj in cell_obj_dict.items():
                self.celltype_cache[cell_id] = cell_obj
        self.jacobian_cache = {}
        self.variant_muts, self.variant_list, self.failed_rsids = read_rsid_parallel(
            self.genome, variant_list, num_workers=self.num_workers
        )
        self.variant_to_genes = self.filter_variant_to_genes_map(variant_to_genes)
        self.variant_to_normal_variants = {}

        all_variant_mut_df_col = [self.variant_muts.df]
        all_failed_variant_list = []
        for rsid in self.variant_list:
            normal_variants_muts, processed_normal_variants, failed_normal_variants = (
                self.get_nearby_variants(rsid)
            )
            all_variant_mut_df_col.append(normal_variants_muts.df)
            self.variant_to_normal_variants[rsid] = processed_normal_variants
            all_failed_variant_list += failed_normal_variants

        self.all_variant_mut_df = pd.concat(
            all_variant_mut_df_col, ignore_index=True
        ).drop_duplicates()
        self.all_failed_variant_list = list(set(all_failed_variant_list))
        self.motif_diff_df = self.generate_motif_diff_df(save_motif_df=True)

        with Path(f"{self.output_dir}/logs/failed_rsid_not_in_ref.txt").open("w") as f:
            for item in self.all_failed_variant_list:
                f.write(item)
                f.write("\n")

    def filter_variant_to_genes_map(
        self, variant_to_genes
    ) -> dict[tuple[str, str], list[str]]:
        variant_to_genes = {
            rsid: gene
            for rsid, gene in variant_to_genes.items()
            if rsid in self.variant_list
        }
        celltype_specific_variant_to_genes = {}
        for variant in variant_to_genes:
            for cell_id in self.celltype_list:
                nearby_genes = self.get_nearby_genes(variant, cell_id)
                celltype_specific_variant_to_genes[(variant, cell_id)] = nearby_genes
        return celltype_specific_variant_to_genes

    def load_normal_filter_normal_variants(self, normal_variants_path) -> pd.DataFrame:
        """
        Load the normal variants and filter them.

        Parameters
        ----------
        normal_variants_path : str
            The path to the normal variants file.

        Returns
        -------
        pd.DataFrame
            The filtered normal variants.
        """
        normal_variants = pd.read_csv(
            normal_variants_path, sep="\t", comment="#", header=None
        )
        normal_variants.columns = [
            "Chromosome",
            "Start",
            "RSID",
            "Ref",
            "Alt",
            "Qual",
            "Filter",
            "Info",
        ]
        normal_variants["End"] = normal_variants.Start
        normal_variants["Start"] = normal_variants.Start - 1
        normal_variants = normal_variants[
            [
                "Chromosome",
                "Start",
                "End",
                "RSID",
                "Ref",
                "Alt",
                "Qual",
                "Filter",
                "Info",
            ]
        ]
        normal_variants = normal_variants.query("Ref.str.len()==1 & Alt.str.len()==1")
        normal_variants["AF"] = normal_variants.Info.transform(
            lambda x: float(re.findall(r"AF=([0-9e\-\.]+)", x)[0])
        )
        return normal_variants

    def generate_motif_diff_df(self, *, save_motif_df: bool = True):
        """
        Generate the motif difference dataframe.

        Parameters
        ----------
        save_motif_df : bool, optional
            Whether to save the motif difference dataframe, by default True

        Returns
        -------
        pd.DataFrame
            The motif difference dataframe.
        """
        variants_rsid = self.all_variant_mut_df.copy()
        variants_rsid = variants_rsid.dropna()
        variants_rsid = variants_rsid.drop_duplicates(subset="RSID", keep="first")
        variants_rsid = Mutations(self.genome, variants_rsid)
        motif_diff = variants_rsid.get_motif_diff(self.motif)
        motif_diff_df = pd.DataFrame(
            (motif_diff["Alt"].values - motif_diff["Ref"].values),
            index=variants_rsid.df.RSID.values,
            columns=self.motif.cluster_names,
        )

        if save_motif_df:
            motif_diff_df.to_csv(Path(self.output_dir) / "motif_diff_df.csv")
        self.motif_diff_df = motif_diff_df
        return motif_diff_df

    def get_variant_score(self, args_tuple) -> pd.DataFrame:
        """
        Get the score for a variant.

        Parameters
        ----------
        args_tuple : tuple
            The tuple of arguments.

        Returns
        -------
        pd.DataFrame
            The score for the variant.
        """
        variant, gene, cell_id = args_tuple
        variant_df = self.all_variant_mut_df[
            self.all_variant_mut_df["RSID"] == variant
        ].iloc[0]
        motif_diff_score = self.motif_diff_df.loc[variant]
        # Use cached cell type
        if cell_id in self.celltype_cache:
            cell = self.celltype_cache[cell_id]
        else:
            cell = GETCellType(cell_id, self.get_config)
            self.celltype_cache[cell_id] = cell
        # Use cached jacobian
        if (cell_id, gene) in self.jacobian_cache:
            motif_importance = self.jacobian_cache[(cell_id, gene)]
        else:
            motif_importance = cell.get_gene_jacobian_summary(gene, "motif")[0:-1]
            self.jacobian_cache[(cell_id, gene)] = motif_importance

        diff = motif_diff_score.copy().values
        combined_score = diff * motif_importance.values
        combined_score = pd.Series(combined_score, index=motif_diff_score.index.values)
        combined_score = pd.DataFrame(combined_score, columns=["score"])
        combined_score["gene"] = gene
        combined_score["variant"] = variant_df.RSID
        combined_score["chrom"] = variant_df.Chromosome
        combined_score["pos"] = variant_df.Start
        combined_score["ref"] = variant_df.Ref
        combined_score["alt"] = variant_df.Alt
        combined_score["celltype"] = self.celltype_annot_dict[cell.celltype]
        combined_score["diff"] = diff
        combined_score["motif_importance"] = motif_importance.values
        return combined_score

    def get_scores_for_single_risk_variant(self, variant) -> pd.DataFrame:
        """
        Get the scores for a single risk variant.

        Parameters
        ----------
        variant : str
            The variant.

        Returns
        -------
        pd.DataFrame
            The scores for the variant.
        """
        variants_to_run = [variant] + self.variant_to_normal_variants[variant]

        scores = []
        failed_args = []
        for cur_celltype in self.celltype_list:
            gene_set = self.variant_to_genes[(variant, cur_celltype)]
            for cur_variant in variants_to_run:
                for cur_gene in gene_set:
                    try:
                        scores.append(
                            self.get_variant_score(
                                (cur_variant, cur_gene, cur_celltype)
                            )
                        )
                    except Exception:
                        failed_args.append((cur_variant, cur_gene, cur_celltype))

        if scores:
            scores = pd.concat(scores, axis=0)
            scores.reset_index().to_feather(
                f"{self.output_dir}/feather/{variant}.feather"
            )
            scores.reset_index().to_csv(f"{self.output_dir}/csv/{variant}.csv")
        return failed_args

    def get_all_variant_scores(self) -> pd.DataFrame:
        """
        Get the scores for all variants.

        Returns
        -------
        pd.DataFrame
            The scores for all variants.
        """
        scores = []
        failed_args = []

        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self.num_workers
        ) as executor:
            future_to_score = {
                executor.submit(
                    self.get_scores_for_single_risk_variant, variant
                ): variant
                for variant in self.variant_list
            }
            for future in concurrent.futures.as_completed(future_to_score):
                failed_args.append(future.result())

        flat_failed_args = [
            failed_tuple
            for failed_tuple_subgroup in failed_args
            for failed_tuple in failed_tuple_subgroup
        ]
        with Path(f"{self.output_dir}/logs/failed_args_no_overlap_with_annot.txt").open(
            "w"
        ) as f:
            for item in flat_failed_args:
                f.write(str(item))
                f.write("\n")
        # read all scores and concat
        for variant in self.variant_list:
            with contextlib.suppress(Exception):
                scores.append(
                    pd.read_feather(f"{self.output_dir}/feather/{variant}.feather")
                )
        scores = pd.concat(scores, axis=0)
        return scores

    def get_nearby_variants(
        self, variant, distance=2000
    ) -> tuple[Mutations, list[str], list[str]]:
        """
        Get the nearby variants for a given variant.

        Parameters
        ----------
        variant : str
            The variant.
        distance : int, optional
            The distance to search for nearby variants, by default 2000

        Returns
        -------
        tuple[Mutations, list[str], list[str]]
            The nearby variants.
        """
        chrom = self.variant_muts.df.query(f'RSID=="{variant}"')["Chromosome"].values[0]
        start = (
            self.variant_muts.df.query(f'RSID=="{variant}"')["Start"].values[0]
            - distance
        )
        end = (
            self.variant_muts.df.query(f'RSID=="{variant}"')["End"].values[0] + distance
        )
        filename = f"https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.{chrom}.vcf.bgz"
        query = f"{chrom}:{start}-{end}"
        command = ["tabix", filename, query]
        result = subprocess.run(command, stdout=subprocess.PIPE)
        result_lines = result.stdout.decode("utf-8").strip().split("\n")

        processed_rsids, failed_rsids = [], []
        chrom_list, start_list, end_list, ref_list, alt_list, rsid_list, af_list = (
            [],
            [],
            [],
            [],
            [],
            [],
            [],
        )
        for line in result_lines:
            chrom, pos, normal_rsid, ref, alt, qual, filter, info = line.split(
                "\t", maxsplit=7
            )

            # Parse out allele frequency
            af = info.split(";")[2]
            af = float(af.split("=")[1]) if af.startswith("AF=") else None

            # Filter by PASS and allele frequency
            if filter == "PASS" and af is not None and af > 1e-2:
                normal_rsid = normal_rsid.split(";")[0]
                processed_rsids.append(normal_rsid)
                chrom_list.append(chrom)
                start_list.append(int(pos) - 1)
                end_list.append(int(pos))
                ref_list.append(ref)
                alt_list.append(alt.split(",")[0])
                rsid_list.append(normal_rsid)
                af_list.append(af)

        df = pd.DataFrame.from_dict(
            {
                "Chromosome": chrom_list,
                "Start": start_list,
                "End": end_list,
                "Ref": ref_list,
                "Alt": alt_list,
                "RSID": rsid_list,
            }
        )
        if self.debug:
            processed_rsids = processed_rsids[:2]
            df = df.iloc[:2]

        if len(df) > 0:
            return [
                Mutations(
                    self.genome,
                    df[["Chromosome", "Start", "End", "Ref", "Alt", "RSID"]],
                ),
                processed_rsids,
                failed_rsids,
            ]
        else:
            return [Mutations(self.genome, None), processed_rsids, failed_rsids]

    def get_nearby_genes(self, variant, cell_id, distance=2000000) -> list[str]:
        """
        Get the nearby genes for a given variant.

        Parameters
        ----------
        variant : str
            The variant.
        cell_id : str
            The cell type.
        distance : int, optional
            The distance to search for nearby genes, by default 2000000

        Returns
        -------
        list[str]
            The nearby genes.
        """
        if cell_id in self.celltype_cache:
            cell = self.celltype_cache[cell_id]
        else:
            cell = GETCellType(cell_id, self.get_config)
            self.celltype_cache[cell_id] = cell

        variant_row = self.variant_muts.df[
            self.variant_muts.df["RSID"] == variant
        ].iloc[0]
        chrom = variant_row["Chromosome"]
        pos = variant_row["Start"]

        # Filter genes directly using boolean masks
        genes = cell.gene_annot[
            (cell.gene_annot["Chromosome"] == chrom)
            & (cell.gene_annot["Start"] > pos - distance)
            & (cell.gene_annot["Start"] < pos + distance)
        ]
        return np.unique(genes.gene_name.values).tolist()


class GETHydraCellMutCollection(CellMutCollection):
    def __init__(self, cfg: DictConfig):
        self.cfg = cfg
        self.output_dir = (
            Path(cfg.machine.output_dir)
            / cfg.run.project_name
            / cfg.run.run_name
            / "variant_analysis"
        )
        self.get_config_path = "local_interpret"
        self.genome_path = cfg.machine.fasta_path
        self.num_workers = cfg.machine.num_workers
        self.debug = False
        hydra_celltype = GETHydraCellType.from_config(cfg)
        self.celltype_cache = {hydra_celltype.celltype: hydra_celltype}
        self.celltype_list = [hydra_celltype.celltype]
        self.celltype_annot_dict = {hydra_celltype.celltype: hydra_celltype.celltype}
        self.setup_directories()
        self.setup_genome_and_motif()
        self.setup_variants()

        self.jacobian_cache = {}

    def setup_directories(self):
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        Path(f"{self.output_dir}/csv").mkdir(parents=True, exist_ok=True)
        Path(f"{self.output_dir}/feather").mkdir(parents=True, exist_ok=True)
        Path(f"{self.output_dir}/logs").mkdir(parents=True, exist_ok=True)

    def setup_genome_and_motif(self):
        self.genome = Genome("hg38")
        self.motif = NrMotifV1.load_from_pickle()

    def setup_variants(self):
        genes_list = self.cfg.task.gene_list.split(",")
        self.variant_list = self.cfg.task.mutations.split(",")
        self.variant_to_genes = {rsid: genes_list for rsid in self.variant_list}

        self.variant_muts, self.variant_list, self.failed_rsids = read_rsid_parallel(
            self.genome, self.variant_list, num_workers=self.num_workers
        )
        self.variant_to_genes = self.filter_variant_to_genes_map(self.variant_to_genes)

        self.setup_normal_variants()
        self.generate_motif_diff_df(save_motif_df=True)

    def setup_normal_variants(self):
        self.variant_to_normal_variants = {}
        all_variant_mut_df_col = [self.variant_muts.df]
        all_failed_variant_list = []

        for rsid in self.variant_list:
            normal_variants_muts, processed_normal_variants, failed_normal_variants = (
                self.get_nearby_variants(rsid)
            )
            all_variant_mut_df_col.append(normal_variants_muts.df)
            self.variant_to_normal_variants[rsid] = processed_normal_variants
            all_failed_variant_list += failed_normal_variants

        self.all_variant_mut_df = pd.concat(
            all_variant_mut_df_col, ignore_index=True
        ).drop_duplicates()
        self.all_failed_variant_list = list(set(all_failed_variant_list))

        with Path(f"{self.output_dir}/logs/failed_rsid_not_in_ref.txt").open("w") as f:
            for item in self.all_failed_variant_list:
                f.write(f"{item}\n")
