"""Define Sample class and subclasses specific to various input formats.

Classes defined herein include:
* Sample
* TextSample
* VCFSample

"""

from __future__ import annotations

import argparse
import contextlib
import csv
import gzip
import logging
import os
import re
from collections import defaultdict
from operator import attrgetter
from typing import Literal, cast

import pandas as pd

from yhaplo import node as node_module  # noqa F401
from yhaplo import snp as snp_module  # noqa F401
from yhaplo import tree as tree_module  # noqa F401
from yhaplo.config import IID_TYPE, Config
from yhaplo.utils.optional_dependencies import (
    check_vcf_dependencies,
)
from yhaplo.utils.vcf import check_vcf_index

with contextlib.suppress(ImportError):
    from pysam import VariantFile

logger = logging.getLogger(__name__)


def call_haplogroups_from_config(config: Config) -> pd.DataFrame:
    """Call haplogroups from a Config instance.

    Parameters
    ----------
    config : Config
        Yhaplo Config instance.

    Returns
    -------
    haplogroup_df : pd.DataFrame
        DataFrame of haplogroup calling results.
        Index: Individual identifier.
        Columns:
        - hg_snp_obs: Haplogroup using a variant of representative-SNP form.
               Rather than using one representative SNP per haplogroup,
               use the most highly ranked SNP this individual was observed
               to carry in the derived state.
        - hg_snp: Haplogroup in representative-SNP form (e.g., "Q-M3").
        - ycc_haplogroup: Haplogroup using YCC nomenclature (e.g., "Q1a2a1a1").

    """
    if config.run_from_sample_major_txt:
        TextSample.call_haplogroups(config)
    elif config.run_from_vcf:
        VCFSample.call_haplogroups(config)
    else:
        logger.info("Mode: No input data\n")
        tree_module.Tree(config)

    haplogroup_df = Sample.haplogroup_df()

    return haplogroup_df


class Sample:
    """Class representing an individual.

    Attributes
    ----------
    iid : IID_TYPE
        Individual identifier.
    haplogroup_node : Node | None
        Node corresponding to haplogroup, once called.
    most_derived_snp : SNP | None
        Most derived SNP observed
    der_snp_list : list[SNP]
        List of SNPs observed in the derived state.
    anc_snp_list: list[SNP]
        List of SNPs observed in the ancestral state.
    anc_der_count_tuples : list[tuple[Node, int, int]]
        List of tuples, each of which includes a node and counts of SNPs observed
        in the ancestral and derived states.

    """

    config: Config
    args: argparse.Namespace
    tree: tree_module.Tree

    num_assigned = 0
    num_root_calls = 0
    sample_list: list[Sample] = []

    def __init__(self, iid: IID_TYPE):
        """Instantiate Sample.

        Parameters
        ----------
        iid : IID_TYPE
            Individual identifier.

        """
        self.iid = iid
        self.haplogroup_node: node_module.Node | None = None
        self.most_derived_snp: snp_module.SNP | None = None
        self.der_snp_list: list[snp_module.SNP]
        self.anc_snp_list: list[snp_module.SNP]
        self.anc_der_count_tuples: list[tuple[node_module.Node, int, int]]

        type(self).sample_list.append(self)

    def __repr__(self) -> str:
        """Return string representation."""

        return (
            f"<{__name__}.{self.__class__.__name__}: "
            f'iid={str(self.iid)}, hg_snp_obs="{self.hg_snp_obs}", '
            f'hg_snp="{self.hg_snp}", haplogroup="{self.haplogroup}">'
        )

    def __str__(self) -> str:
        """Return printable string representation."""

        return (
            f"{str(self.iid):8s} {self.hg_snp_obs:15s} "
            f"{self.hg_snp:15s} {self.haplogroup:25s}"
        )

    # Haplogroup calling
    # ----------------------------------------------------------------------
    def call_haplogroup(self) -> None:
        """Call haplogroup."""

        type(self).num_assigned += 1
        tree = type(self).tree
        (
            path,
            self.anc_snp_list,
            self.anc_der_count_tuples,
        ) = tree.identify_phylogenetic_path(self)
        self.der_snp_list = path.der_snp_list
        self.most_derived_snp = path.most_derived_snp

        if self.most_derived_snp:
            self.haplogroup_node = self.most_derived_snp.node
        else:
            type(self).num_root_calls += 1
            self.haplogroup_node = type(self).tree.root

        with contextlib.suppress(NotImplementedError):
            self.fix_haplogroup_if_artifact()

        self.write_real_time_output()
        self.purge_data()

    def get_genotype(self, position: int) -> str:
        """Return consensus genotype for position. Subclasses must override."""

        raise NotImplementedError

    def fix_haplogroup_if_artifact(self) -> None:
        """Fix artifactual haplogroup assignments. Subclasses may override."""

        raise NotImplementedError

    def write_real_time_output(self) -> None:
        """Write real-time output if requested."""

        args, config = type(self).args, type(self).config

        if args.write_haplogroups_real_time:
            config.haplogroup_real_time_file.write(
                f"{str(self)} {self.haplogroup_dfs_rank:5d}\n"
            )

        if args.haplogroup_to_list_genotypes_for:
            config.hg_genos_file.write(f"{self.str_compressed}\n\n")

    def purge_data(self) -> None:
        """Clear data structures if no longer needed."""

        args = type(self).args

        if not (
            args.write_der_snps
            or args.write_der_snps_detail
            or args.write_haplogroup_paths
            or args.write_haplogroup_paths_detail
        ):
            self.der_snp_list.clear()

        if not (args.write_anc_snps or args.write_anc_snps_detail):
            self.anc_snp_list.clear()

        if not self.args.write_anc_der_counts:
            self.anc_der_count_tuples.clear()

    # Haplogroup properties
    # ----------------------------------------------------------------------
    @property
    def haplogroup(self) -> str:
        """Return haplogroup using YCC nomenclature (e.g., "Q1a2a1a1")."""

        if self.haplogroup_node is None:
            raise RuntimeError(f"Haplogroup not yet computed for {self.iid}")

        return self.haplogroup_node.haplogroup

    @property
    def hg_snp(self) -> str:
        """Return haplogroup in representative-SNP form (e.g., "Q-M3")."""

        if self.haplogroup_node is None:
            raise RuntimeError(f"Haplogroup not yet computed for {self.iid}")

        return self.haplogroup_node.hg_snp

    @property
    def hg_trunc(self) -> str:
        """Return haplogroup in truncated form (e.g., "Q")."""

        if self.haplogroup_node is None:
            raise RuntimeError(f"Haplogroup not yet computed for {self.iid}")

        return self.haplogroup_node.hg_trunc

    @property
    def hg_snp_obs(self) -> str:
        """Return haplogroup using a variant of representative-SNP form.

        Rather than using one representative SNP per haplogroup,
        use the most highly ranked SNP this individual was observed to carry
        in the derived state.

        """
        if self.haplogroup_node is None:
            raise RuntimeError(f"Haplogroup not yet computed for {self.iid}")

        if self.most_derived_snp:
            hg_snp_obs = self.most_derived_snp.hg_snp
        elif self.haplogroup_node:
            hg_snp_obs = type(self).tree.root.haplogroup
        else:
            raise RuntimeError(f"Haplogroup not yet computed for {self.iid}")

        return hg_snp_obs

    @property
    def haplogroup_dfs_rank(self) -> int:
        """Return depth-first-search ranking of haplogroup node."""

        if self.haplogroup_node is None:
            raise RuntimeError(f"Haplogroup not yet computed for {self.iid}")

        haplogroup_dfs_rank = self.haplogroup_node.dfs_rank
        return haplogroup_dfs_rank

    @property
    def haplogroup_dict(self) -> dict[str, str | int]:
        """Return dictionary with various representations of the haplogroup call.

        Returns
        -------
        haplogroup_dict : dict[str, str | int]
            Keys:
            - "iid": Individual identifier.
            - "hg_snp_obs": Haplogroup using a variant of representative-SNP form.
                   Rather than using one representative SNP per haplogroup,
                   use the most highly ranked SNP this individual was observed
                   to carry in the derived state.
            - "hg_snp": Haplogroup in representative-SNP form (e.g., "Q-M3").
            - "ycc_haplogroup": Haplogroup using YCC nomenclature (e.g., "Q1a2a1a1").

        """
        haplogroup_dict = {
            "iid": self.iid,
            "hg_snp_obs": self.hg_snp_obs,
            "hg_snp": self.hg_snp,
            "ycc_haplogroup": self.haplogroup,
        }

        return haplogroup_dict

    # String-representation properties and methods
    # ----------------------------------------------------------------------
    @property
    def str_compressed(self) -> str:
        """Return compressed string representation."""

        str_compressed = re.sub(r"\s+", " ", str(self))
        return str_compressed

    @property
    def str_simple(self) -> str:
        """Return string representation with just iid, haplogroup, and hg_snp."""

        return f"{str(self.iid):8s} {self.haplogroup:25s} {self.hg_snp:15s}"

    @property
    def str_for_counts(self) -> str:
        """Return string representation for ancestral/derived counts output."""

        left_part = f"{str(self.iid):8s} {self.haplogroup}"
        right_part = f"{self.hg_snp_obs} {self.hg_snp}"
        str_for_counts = f"{left_part} | {right_part}"

        return str_for_counts

    # String-representation methods
    # ----------------------------------------------------------------------
    def str_snps(
        self,
        allele_state: Literal["derived", "ancestral"] = "derived",
    ) -> str:
        """Return string representation with derived or ancestral SNPs."""

        if allele_state == "derived":
            snp_list = self.der_snp_list
        elif allele_state == "ancestral":
            snp_list = self.anc_snp_list
        else:
            raise ValueError(
                f'allele_state must be "ancestral" or "derived", not "{allele_state}"'
            )

        snp_list_string = " ".join(snp.str_short for snp in snp_list)
        str_snps = f"{self.str_simple} | {snp_list_string}"

        return str_snps

    def str_haplogroup_path(
        self,
        include_snps: bool = False,
    ) -> str:
        """Return string representation with haplogroup path."""

        if self.most_derived_snp:
            snp_label_list_dict = defaultdict(list)
            for snp in self.der_snp_list:
                snp_label_list_dict[snp.node.haplogroup].append(snp.label_cleaned)

            path_info_list = []
            for node in self.most_derived_snp.back_trace_path():
                if node.haplogroup in snp_label_list_dict:
                    snp_label_list = snp_label_list_dict[node.haplogroup]
                    num_snps = len(snp_label_list)
                    path_info = f"{node.haplogroup}:{num_snps}"
                    if include_snps:
                        path_info = f"{path_info}:{','.join(snp_label_list)}"

                    path_info_list.append(path_info)

            haplogroup_path = " ".join(path_info_list)
        else:
            haplogroup_path = ""

        str_haplogroup_path = f"{self.str_simple} | {haplogroup_path}"
        return str_haplogroup_path

    # Class methods: results
    # ----------------------------------------------------------------------
    @classmethod
    def haplogroup_df(cls) -> pd.DataFrame:
        """Return DataFrame of haplogroup calling results.

        Returns
        -------
        haplogroup_df : pd.DataFrame
            DataFrame of haplogroup calling results.
            Index: Individual identifier.
            Columns:
            - hg_snp_obs: Haplogroup using a variant of representative-SNP form.
                   Rather than using one representative SNP per haplogroup,
                   use the most highly ranked SNP this individual was observed
                   to carry in the derived state.
            - hg_snp: Haplogroup in representative-SNP form (e.g., "Q-M3").
            - ycc_haplogroup: Haplogroup using YCC nomenclature (e.g., "Q1a2a1a1").

        """
        haplogroup_df = pd.DataFrame(
            [
                (sample.iid, sample.hg_snp_obs, sample.hg_snp, sample.haplogroup)
                for sample in cls.sample_list
            ],
            columns=["iid", "hg_snp_obs", "hg_snp", "ycc_haplogroup"],
        ).set_index("iid")

        return haplogroup_df

    # Class methods: configuration
    # ----------------------------------------------------------------------
    @classmethod
    def configure(cls, config: Config) -> None:
        """Configure class.

        This method will:
        - Store a Config instance
        - Store command-line arguments
        - Build and store a Tree instance, if not already present
        - Clear the sample list
        - Run some error checking
        - Log some information

        """
        cls.config = config
        cls.args = cls.config.args
        if not hasattr(cls, "tree"):
            cls.tree = tree_module.Tree(config)
        else:
            logger.info("Tree: Previously constructed\n")
            logger.info("Variants: Previously loaded\n")

        cls.num_assigned = 0
        cls.num_root_calls = 0
        cls.sample_list.clear()
        cls.check_number_of_run_modes()

        if cls.args.write_haplogroups_real_time:
            logger.info(
                f"\nWill write haplogroups as they are called:\n"
                f"    {cls.config.haplogroup_real_time_fp}\n"
                "Note: This file includes DFS rank, so it can be sorted "
                "ex post facto with:\n"
                f"    sort -nk5 {cls.config.haplogroup_real_time_fp}\n"
            )

        if config.run_from_sample_major_txt:
            input_description = f"sample-major text file:\n    {cls.args.data_fp}"
        elif config.run_from_vcf:
            input_description = f"variant-major VCF/BCF file:\n    {cls.args.data_fp}"
        else:
            assert config.iid_to_ablock is not None
            num_ablocks = len(config.iid_to_ablock)
            plural_s = "s" if num_ablocks > 1 else ""
            input_description = f"[{num_ablocks}] 23andMe ablock{plural_s}..."

        logger.info(f"\nGenotypes\n\nLoading genotypes from {input_description}\n")

    @classmethod
    def check_number_of_run_modes(cls) -> None:
        """Check the number of run modes.

        Raises
        ------
        ValueError
            When more then one run mode have been selected.

        """
        number_of_run_modes_selected = (
            cls.config.run_from_sample_major_txt
            + cls.config.run_from_vcf
            + cls.config.run_from_ablocks
        )
        if number_of_run_modes_selected > 1:
            raise ValueError(
                "Expecting no more than one run mode\n"
                f"    {number_of_run_modes_selected} selected\n"
            )

    # Class methods: class variable mutators
    # ----------------------------------------------------------------------
    @classmethod
    def sort_sample_list(cls) -> None:
        """Sort sample list by haplogroup, then by iid."""

        cls.sample_list.sort(key=attrgetter("iid"))
        if cls.sample_list[0].haplogroup_node:
            cls.sample_list.sort(key=attrgetter("haplogroup_dfs_rank"))

    # Class methods: output writers
    # ----------------------------------------------------------------------
    @classmethod
    def write_results(cls) -> None:
        """Sort samples, write results, and close optional real-time output files."""

        cls.sort_sample_list()

        logger.info(
            "\nHaplogroups\n\n"
            f"    {cls.num_assigned:8d} assigned\n"
            f"    {cls.num_root_calls:8d} assigned to root haplogroup: "
            f"{cls.tree.root.haplogroup}\n"
        )
        if cls.num_root_calls > 0:
            logger.warning(
                "WARNING. If the dataset does not include fixed reference sites,\n"
                "         re-run with alternative root (e.g., with: -r A0-T).\n\n\n"
            )

        if not cls.config.suppress_output:  # Use str(sample)
            logger.info("\nOutput\n")
            cls.write_haplogroups()

        if cls.args.write_anc_der_counts:  # Use sample.str_for_counts
            cls.write_anc_der_counts()

        if cls.args.write_haplogroup_paths_detail:  # Use sample.str_haplogroup_path()
            cls.write_haplogroup_paths(include_snps=True)
        elif cls.args.write_haplogroup_paths:  # Use sample.str_haplogroup_path()
            cls.write_haplogroup_paths()

        if cls.args.write_der_snps:  # Use sample.str_snps()
            cls.write_snps(allele_state="derived")

        if cls.args.write_der_snps_detail:  # Use sample.str_compressed
            cls.write_snps_detail()

        if cls.args.write_anc_snps:  # Use sample.str_snps()
            cls.write_snps(allele_state="ancestral")

        if cls.args.write_anc_snps_detail:  # Use sample.str_compressed
            cls.write_snps_detail(ancestral=True)

        cls.config.close_real_time_output_files()

    @classmethod
    def write_haplogroups(cls) -> None:
        """Write haplogroup of each sample."""

        with open(cls.config.haplogroup_calls_fp, "w") as haplogroup_calls_file:
            for sample in cls.sample_list:
                haplogroup_calls_file.write(f"{sample}\n")

        logger.info(
            f"Wrote called haplogroups:\n    {cls.config.haplogroup_calls_fp}\n"
        )

    @classmethod
    def write_anc_der_counts(cls) -> None:
        """Write counts of ancestral and derived alleles encountered.

        This includes each visited node,
        other than those with no ancestral or derived alleles.

        """
        with open(cls.config.counts_anc_der_fp, "w") as counts_anc_der_file:
            for sample in cls.sample_list:
                for node, num_ancestral, num_derived in sample.anc_der_count_tuples:
                    if num_ancestral > 0 or num_derived > 0:
                        counts_anc_der_file.write(
                            f"{str(sample.iid):8s} {node.label:25s} "
                            f"{num_ancestral:3d} {num_derived:3d}\n"
                        )

                counts_anc_der_file.write(f"{sample.str_for_counts}\n\n")

        logger.info(
            "Wrote counts of ancestral and derived alleles encountered:\n"
            f"    {cls.config.counts_anc_der_fp}\n"
        )

    @classmethod
    def write_haplogroup_paths(
        cls,
        include_snps: bool = False,
    ) -> None:
        """Write haplogroup path for each sample."""

        with open(cls.config.haplogroup_paths_fp, "w") as haplogroup_paths_file:
            for sample in cls.sample_list:
                path = sample.str_haplogroup_path(include_snps)
                haplogroup_paths_file.write(f"{path}\n")

        logger.info(
            "Wrote paths with counts of derived SNPs observed:\n"
            f"    {cls.config.haplogroup_paths_fp}\n"
        )

    @classmethod
    def write_snps(
        cls,
        allele_state: Literal["derived", "ancestral"] = "derived",
    ) -> None:
        """Write list of derived or ancestral alleles encountered.

        Repeat for each sample.

        """
        if allele_state == "derived":
            snp_fp = cls.config.der_snps_fp
            type_of_snps = "derived SNPs on path"
        elif allele_state == "ancestral":
            snp_fp = cls.config.anc_snps_fp
            type_of_snps = "ancestral SNPs encountered in search"
        else:
            raise ValueError(
                f'allele_state must be "ancestral" or "derived", not "{allele_state}"'
            )

        with open(snp_fp, "w") as snp_file:
            for sample in cls.sample_list:
                snp_file.write(f"{sample.str_snps(allele_state)}\n")

        logger.info(f"Wrote lists of {type_of_snps}:\n    {snp_fp}\n")

    @classmethod
    def write_snps_detail(
        cls,
        ancestral: bool = False,
    ) -> None:
        """Write detailed information about derived or ancestral alleles observed.

        Repeat for each sample.

        """
        if ancestral:
            snp_detail_fp = cls.config.anc_snps_detail_fp
            type_of_snps = "ancestral SNP encountered"
        else:
            snp_detail_fp = cls.config.der_snps_detail_fp
            type_of_snps = "derived SNP on path"

        with open(snp_detail_fp, "w") as snp_detail_file:
            for sample in cls.sample_list:
                snp_detail_file.write(f"{sample.str_compressed}\n")
                snp_list = sample.anc_snp_list if ancestral else sample.der_snp_list
                for snp in snp_list:
                    snp_detail_file.write(f"{str(sample.iid):8s} {snp}\n")

                snp_detail_file.write("\n")

        logger.info(
            f"Wrote detailed information about each {type_of_snps}:\n"
            f"    {snp_detail_fp}\n"
        )


class TextSample(Sample):
    """Class representing an individual whose data are in a sample-major text file.

    Expected input format:
    - Row 1: Physical coordinates
    - Column 1: Individual identifiers

    Attributes
    ----------
    genotypes : list[str]
        List of genotypes.

    """

    position_to_column_index: dict[int, int]

    def __init__(
        self,
        iid: IID_TYPE,
        genotypes: list[str],
    ):
        """Instantiate TextSample.

        Parameters
        ----------
        iid : IID_TYPE
            Individual identifier.
        genotypes : list[str]
            List of genotypes.

        """
        super().__init__(iid)
        self.genotypes = genotypes

    def get_genotype(self, position: int) -> str:
        """Return genotype for position."""

        try:
            genotype = self.genotypes[type(self).position_to_column_index[position]]
        except KeyError:
            genotype = Config.missing_genotype

        return genotype

    def purge_data(self) -> None:
        """Clear genotype data and other data structures if no longer needed."""

        super().purge_data()
        self.genotypes.clear()

    @classmethod
    def call_haplogroups(cls, config: Config) -> None:
        """Call haplogroups from sample-major text file."""

        logger.info("Mode: Sample-major text\n")
        cls.configure(config)
        geno_file = (
            open(cls.args.data_fp)  # noqa SIM115
            if os.path.splitext(cls.args.data_fp)[1] != ".gz"
            else gzip.open(cls.args.data_fp, "rt")  # noqa SIM115
        )
        geno_reader = csv.reader(geno_file, delimiter="\t")
        cls.position_to_column_index = {
            position: column_index
            for column_index, position in enumerate(map(int, next(geno_reader)[1:]))
            if position in cls.tree.snp_pos_set
        }
        for genotypes in geno_reader:
            iid, genotypes = genotypes[0], genotypes[1:]
            if cls.args.single_sample_id is None or iid == cls.args.single_sample_id:
                text_sample = cls(iid, genotypes)
                text_sample.call_haplogroup()

        geno_file.close()

        cls.write_results()


class VCFSample(Sample):
    """Class representing an individual whose data are in a VCF/BCF file.

    Attributes
    ----------
    position_to_genotype : dict[int, str]
        Maps position to genotype.

    """

    def __init__(self, iid: IID_TYPE):
        """Instantiate VCFSample.

        Parameters
        ----------
        iid : IID_TYPE
            Individual identifier.

        """
        super().__init__(iid)
        self.position_to_genotype: dict[int, str] = {}

    def put_genotype(self, position: int, genotype: str) -> None:
        """Store one genotype."""

        self.position_to_genotype[position] = genotype

    def get_genotype(self, position: int) -> str:
        """Return genotype for position."""

        genotype = self.position_to_genotype.get(position, Config.missing_genotype)

        return genotype

    @classmethod
    def call_haplogroups(cls, config: Config) -> None:
        """Call haplogroups from variant-major VCF/BCF file."""

        logger.info("Mode: VCF\n")
        cls.configure(config)
        cls.load_data_from_vcf()

        if (
            cls.args.write_haplogroups_real_time
            or cls.args.haplogroup_to_list_genotypes_for
        ):
            cls.sort_sample_list()

        for vcf_sample in cls.sample_list:
            vcf_sample.call_haplogroup()

        cls.write_results()

    @classmethod
    def load_data_from_vcf(cls) -> None:
        """Load data from VCF."""

        check_vcf_dependencies()
        check_vcf_index(cls.args.data_fp)

        with VariantFile(cls.args.data_fp) as variant_file:
            iids = list(variant_file.header.samples)
            if cls.args.single_sample_id is None:
                for iid in iids:
                    cls(iid)
            else:
                if cls.args.single_sample_id in iids:
                    cls(cls.args.single_sample_id)
                else:
                    raise ValueError(
                        f"{cls.args.single_sample_id} not present in {cls.args.data_fp}"
                    )

            chromosome_set = set(variant_file.header.contigs.keys())
            y_chromosome_set = chromosome_set.intersection(Config.vcf_chrom_label_set)
            if len(y_chromosome_set) == 1:
                chromosome = y_chromosome_set.pop()
            else:
                raise ValueError(
                    f"VCF must include exactly one contig with a label in: "
                    f"{sorted(Config.vcf_chrom_label_set)}\n"
                    f"Observed: {chromosome_set}"
                )

            for variant_record in variant_file.fetch(chromosome):
                if variant_record.pos not in cls.tree.snp_pos_set:
                    continue

                for vcf_sample in cast(list[VCFSample], cls.sample_list):
                    variant_record_sample = variant_record.samples[vcf_sample.iid]
                    unique_alleles_set = set(variant_record_sample.alleles) - {None}
                    if len(unique_alleles_set) == 1:
                        vcf_sample.put_genotype(
                            variant_record.pos,
                            unique_alleles_set.pop(),
                        )
