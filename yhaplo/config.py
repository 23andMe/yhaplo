"""Define Config class, which includes command-line arguments."""

import argparse
import logging
import os
import sys
from collections.abc import Iterable, Mapping
from typing import Optional, Union

import numpy as np
from numpy.typing import NDArray

from yhaplo import __version__
from yhaplo.api.command_line_args import get_command_line_arg_defaults
from yhaplo.utils.loaders import DataFile

DASHED_LINE = "-" * 72 + "\n"
IID_TYPE = Union[int, str]
ABLOCK_TYPE = NDArray[np.uint8]

logger = logging.getLogger(__name__)


class Config:
    """Yhaplo configuration class.

    This class is a container for parameters, constants, filenames, etc.

    """

    # Constants
    # --------------------------------------------------------------------

    # Independent constants
    isogg_date = "2016.01.04"  # Date ISOGG website scraped
    root_haplogroup = "A"  # Haplogroup to associate with root node
    missing_genotype = "."  # For text input
    missing_haplogroup = "."  # For output
    vcf_chrom_label_set = {"Y", "chrY", "24"}  # To restrict VCF input to chrY

    newick_semantic_token_string = "(),:;"  # Used in regex
    alleles_string = "A C G T D I"  # Allowable alleles
    snp_label_letters_rank_string = "M P V L CTS AF B F Page U PF Z SK"
    superfluous_snp_text_list = ["IMS-", "-null"]  # Stripped out of SNP names
    multi_char_hg_trunc_string = (
        "A00 A0-T A0 A1 A1a A1b A1b1 BT CT DE CF GHIJK HIJK IJK IJ LT NO"
    )

    # Derived constants
    newick_semantic_token_set = set(newick_semantic_token_string)
    multi_char_hg_trunc_set = set(multi_char_hg_trunc_string.split())
    multi_char_hg_trunc_max_len = max([len(elem) for elem in multi_char_hg_trunc_set])
    allele_set = set(alleles_string.split())
    homozygous_genotype_set = {f"{allele}{allele}" for allele in allele_set}
    snp_label_letters_rank_dict = {
        letters: rank
        for rank, letters in enumerate(snp_label_letters_rank_string.split())
    }

    # 23andMe-specific constants
    ttam_hg_call_replacement_dict = {"BT": "B"}  # Prevents artifactual calls
    calling_progress_early_set = {100, 500, 1000, 5000}  # For progress messages
    calling_progress_interval = 10_000  # For progress messages
    platforms = ["v1", "v2", "v3", "v4", "v5"]

    # Data files
    # ----------------------------------------------------------------------
    # Tree data
    primary_tree_data_file = DataFile(
        "tree",
        f"y.tree.primary.{isogg_date}.nwk",
        "primary tree topology",
    )

    # Variant data
    isogg_data_file = DataFile("variants", f"isogg.{isogg_date}.txt", "ISOGG SNP data")
    isogg_correct_name_data_file = DataFile("variants", "isogg.correct.name.txt")
    isogg_corrections_data_file_dict = {
        "mutation": DataFile("variants", "isogg.correct.mutation.txt"),
        "position": DataFile("variants", "isogg.correct.position.txt"),
    }
    isogg_omit_name_data_file = DataFile("variants", "isogg.omit.name.txt")
    isogg_omit_pos_str_mutation_data_files = [
        DataFile("variants", "isogg.omit.bad.txt"),
        DataFile("variants", "isogg.omit.bad.23andMe.txt"),
        DataFile("variants", "isogg.omit.branch.conflict.txt"),
        DataFile("variants", "isogg.omit.branch.conflict.23andMe.v5.txt"),
    ]
    isogg_multiallelic_data_file = DataFile("variants", "isogg.multiallelic.txt")
    isogg_rep_snp_data_file = DataFile(
        "variants", "representative.SNPs.isogg.2015_tree.txt"
    )
    other_rep_snp_data_file = DataFile("variants", "representative.SNPs.additional.txt")
    preferred_snp_names_data_file = DataFile("variants", "preferred.snp_names.txt")

    # 23andMe: Block data
    pos_to_block_indexes_data_file = DataFile(
        "block",
        "pos_to_block_indexes.yaml",
        "Position-to-block-index",
        ttam_only=True,
    )

    # 23andMe: Platform positions
    platform_pos_data_subdir = "platform"
    platform_pos_fn_tp = "{platform}.b37.positions.txt"
    platform_qc_exclude_fn_tp = "{platform}.b37.qc.exclude.txt"

    # Example input files
    # ----------------------------------------------------------------------
    repo_root_dir = os.path.dirname(
        os.path.dirname(os.path.realpath(__file__))
    ).removesuffix("/ttam")
    example_genotype_dir = os.path.join(
        repo_root_dir,
        "tests",
        "fixtures",
        "input",
    )
    kgp_subset_fp = os.path.join(example_genotype_dir, "1000Y.subset.genos.txt")
    kgp_single_sample_vcf_fp = os.path.join(example_genotype_dir, "HG01938.vcf.gz")

    # Output files
    # ----------------------------------------------------------------------
    default_out_dir = "output"

    # Phylogenetic info
    aligned_primary_tree_fn = f"y.tree.primary.aligned.ycc.{isogg_date}.nwk"
    ycc_tree_fn = f"y.tree.ycc.{isogg_date}.nwk"
    hg_snp_tree_fn = f"y.tree.hg_snp.{isogg_date}.nwk"
    aligned_ycc_tree_fn = f"y.tree.aligned.ycc.{isogg_date}.nwk"
    aligned_hg_snp_tree_fn = f"y.tree.aligned.hg_snp.{isogg_date}.nwk"
    platform_ycc_tree_fn_tp = f"y.tree.platform.{{platform}}.ycc.{isogg_date}.nwk"
    platform_hg_snp_tree_fn_tp = f"y.tree.platform.{{platform}}.hg_snp.{isogg_date}.nwk"
    bf_tree_fn = f"y.tree.bf.traversal.{isogg_date}.txt"
    df_tree_fn = f"y.tree.df.traversal.{isogg_date}.txt"
    tree_table_fn = f"y.tree.table.{isogg_date}.tsv"
    bf_primary_tree_fn = f"y.tree.primary.bf.traversal.{isogg_date}.txt"
    df_primary_tree_fn = f"y.tree.primary.df.traversal.{isogg_date}.txt"
    cleaned_isogg_fn = f"isogg.snps.cleaned.{isogg_date}.txt"
    unique_isogg_fn = f"isogg.snps.unique.{isogg_date}.txt"
    dropped_isogg_fn = f"isogg.snps.dropped.{isogg_date}.txt"
    multiallelic_found_fn = "multiallelic.pos"

    # Haplogroup calls, log, optional files, 23andMe auxiliary files
    log_fn_tp = "log.{fn_label}.txt"
    haplogroup_calls_fn_tp = "haplogroups.{fn_label}.txt"
    haplogroup_real_time_fn_tp = "haplogroups.real_time.{fn_label}.txt"
    counts_anc_der_fn_tp = "counts.anc_der.{fn_label}.txt"
    haplogroup_paths_fn_tp = "paths.{fn_label}.txt"
    der_snps_fn_tp = "derived.snps.{fn_label}.txt"
    der_snps_detail_fn_tp = "derived.snps.detail.{fn_label}.txt"
    anc_snps_fn_tp = "ancestral.snps.{fn_label}.txt"
    anc_snps_detail_fn_tp = "ancestral.snps.detail.{fn_label}.txt"
    hg_genos_fn_tp = "hg.{{haplogroup}}.genotypes.{fn_label}.txt"

    def __init__(
        self,
        command_line_args: Optional[argparse.Namespace] = None,
        iid_to_ablock: Optional[
            Mapping[
                IID_TYPE,
                Union[bytes, NDArray[np.uint8]],
            ]
        ] = None,
        iid_to_platforms: Optional[Mapping[IID_TYPE, Union[str, Iterable[str]]]] = None,
        suppress_output: bool = False,
        out_dir: Optional[str] = None,
        all_aux_output: bool = False,
        root_logger: Optional[logging.Logger] = None,
    ):
        """Instantiate Config.

        Parameters
        ----------
        command_line_args : argparse.Namespace | None, optional
            Command-line arguments.
        iid_to_ablock : Mapping[IID_TYPE, bytes | ABLOCK_TYPE] | None, optional
            Mapping of individual identifier to 23andMe ablock.
        iid_to_platforms : Mapping[IID_TYPE, str | Iterable[str]] | None, optional
            Mapping of individual identifier to 23andMe genotyping platforms,
            each starting with "v". Values can be a comma-separated string
            or an iterable of strings.
        suppress_output : bool, optional
            When True, do not generate output files.
        out_dir : str | None, optional
            Output directory.
            When not None, override command_line_args value.
        all_aux_output : bool = False, optional
            Generate all auxiliary output.
            When True, override command_line_args value.
        root_logger : logging.Logger | None, optional
            If supplied, add a file handler.
            To populate the log file, set the logging level to INFO or lower.

        """
        self.args = (
            command_line_args
            if command_line_args is not None
            else get_command_line_arg_defaults()
        )
        self.invoked_from_command_line = command_line_args is not None

        self.iid_to_ablock = iid_to_ablock
        self.iid_to_platforms = iid_to_platforms
        if iid_to_ablock is not None and iid_to_platforms is None:
            raise ValueError("Calling from ablocks requires iid_to_platforms")

        self.suppress_output = suppress_output
        self.set_params_general(out_dir, all_aux_output)
        self.set_params_based_on_input_type()
        self.make_output_directories()
        self.set_output_file_paths_and_open_some(root_logger)

        self.log_welcome_message()
        if self.suppress_output:
            self.override_output_generating_args()

    def __repr__(self) -> str:
        """Return string representation."""

        return f"<{__name__}.{self.__class__.__name__}: command_line_args={self.args}>"

    def set_params_general(
        self,
        out_dir: Optional[str],
        all_aux_output: bool,
    ) -> None:
        """Set general parameters."""

        # Run type: Zero or one of these will be set True by downstream methods.
        self.run_from_ablocks = False
        self.run_from_sample_major_txt = False
        self.run_from_vcf = False

        # Output directories
        if out_dir is not None:
            self.out_dir = out_dir
        elif self.args.out_dir is not None:
            self.out_dir = self.args.out_dir
        else:
            self.out_dir = type(self).default_out_dir

        self.phylo_out_dir = self.out_dir

        # Example data
        if self.args.example_text:
            self.args.data_fp = type(self).kgp_subset_fp
        elif self.args.example_vcf:
            self.args.data_fp = type(self).kgp_single_sample_vcf_fp

        if self.args.example_text or self.args.example_vcf:
            check_example_data_availability(self.args.data_fp)
            self.args.all_aux_output = True

        # All auxiliary output option
        if self.args.all_aux_output or all_aux_output:
            self.args.write_anc_der_counts = True
            self.args.write_haplogroup_paths = True
            self.args.write_haplogroup_paths_detail = True
            self.args.write_der_snps = True
            self.args.write_der_snps_detail = True
            self.args.write_anc_snps = True
            self.args.write_anc_snps_detail = True

    def set_params_based_on_input_type(self) -> None:
        """Set parameters based on input type."""

        if self.args.data_fp:
            self.parse_data_fp()
        else:
            self.run_from_ablocks = bool(self.iid_to_ablock)
            self.out_fn_label = ""

    def parse_data_fp(self) -> None:
        """Set parameters based on data file name."""

        sample_major_suffixes = {".genos.txt", ".genos.txt.gz"}
        vcf_suffixes = {".vcf.gz", ".bcf"}
        supported_suffixes = sample_major_suffixes | vcf_suffixes
        suffix = ""
        for supported_suffix in supported_suffixes:
            if self.args.data_fp.endswith(supported_suffix):
                suffix = supported_suffix

        if not suffix:
            raise ValueError(f"Unknown data type: {self.args.data_fp}\n\n")

        self.run_from_sample_major_txt = suffix in sample_major_suffixes
        self.run_from_vcf = suffix in vcf_suffixes
        self.out_fn_label = basename_no_suffix(self.args.data_fp, suffix)
        if self.args.single_sample_id:
            self.out_fn_label = f"{self.out_fn_label}.{self.args.single_sample_id}."

    def make_output_directories(self) -> None:
        """Make output directories."""

        if not self.suppress_output:
            for dir_ in [self.out_dir, self.phylo_out_dir]:
                os.makedirs(dir_, exist_ok=True)

    def set_output_file_paths_and_open_some(
        self,
        root_logger: Optional[logging.Logger] = None,
    ) -> None:
        """Set log and output file paths.

        If (and only if) `root_logger` is supplied, add a file handler.
        In general, library code should not be in the business of adding
        handlers, but we (conditionally) do so here since the log filepath
        is dynamically determined.

        Open output files to which yhaplo will write in real time.

        Parameters
        ----------
        root_logger : logging.Logger | None, optional
            If supplied, add a file handler.

        """
        self.log_fp = self.construct_out_path(type(self).log_fn_tp)
        self.root_logger = root_logger
        if self.root_logger is not None:
            self.root_logger.addHandler(logging.FileHandler(self.log_fp, "w"))

        self.haplogroup_calls_fp = self.construct_out_path(
            type(self).haplogroup_calls_fn_tp
        )

        self.aligned_primary_tree_fp = self.construct_phylo_out_path(
            type(self).aligned_primary_tree_fn
        )
        self.ycc_tree_fp = self.construct_phylo_out_path(type(self).ycc_tree_fn)
        self.hg_snp_tree_fp = self.construct_phylo_out_path(type(self).hg_snp_tree_fn)
        self.aligned_ycc_tree_fp = self.construct_phylo_out_path(
            type(self).aligned_ycc_tree_fn
        )
        self.aligned_hg_snp_tree_fp = self.construct_phylo_out_path(
            type(self).aligned_hg_snp_tree_fn
        )
        self.platform_ycc_tree_fp_tp = self.construct_phylo_out_path(
            type(self).platform_ycc_tree_fn_tp
        )
        self.platform_hg_snp_tree_fp_tp = self.construct_phylo_out_path(
            type(self).platform_hg_snp_tree_fn_tp
        )
        self.bf_tree_fp = self.construct_phylo_out_path(type(self).bf_tree_fn)
        self.df_tree_fp = self.construct_phylo_out_path(type(self).df_tree_fn)
        self.tree_table_fp = self.construct_phylo_out_path(type(self).tree_table_fn)
        self.bf_primary_tree_fp = self.construct_phylo_out_path(
            type(self).bf_primary_tree_fn
        )
        self.df_primary_tree_fp = self.construct_phylo_out_path(
            type(self).df_primary_tree_fn
        )
        self.cleaned_isogg_fp = self.construct_phylo_out_path(
            type(self).cleaned_isogg_fn
        )
        self.unique_isogg_fp = self.construct_phylo_out_path(type(self).unique_isogg_fn)
        self.dropped_isogg_fp = self.construct_phylo_out_path(
            type(self).dropped_isogg_fn
        )
        self.multiallelic_found_fp = self.construct_phylo_out_path(
            type(self).multiallelic_found_fn
        )

        if self.args:
            if self.args.write_anc_der_counts:
                self.counts_anc_der_fp = self.construct_out_path(
                    type(self).counts_anc_der_fn_tp
                )

            if (
                self.args.write_haplogroup_paths
                or self.args.write_haplogroup_paths_detail
            ):
                self.haplogroup_paths_fp = self.construct_out_path(
                    type(self).haplogroup_paths_fn_tp
                )

            if self.args.write_der_snps:
                self.der_snps_fp = self.construct_out_path(type(self).der_snps_fn_tp)

            if self.args.write_der_snps_detail:
                self.der_snps_detail_fp = self.construct_out_path(
                    type(self).der_snps_detail_fn_tp
                )

            if self.args.write_anc_snps:
                self.anc_snps_fp = self.construct_out_path(type(self).anc_snps_fn_tp)

            if self.args.write_anc_snps_detail:
                self.anc_snps_detail_fp = self.construct_out_path(
                    type(self).anc_snps_detail_fn_tp
                )

            # Files to write to in real time. Open now.
            if self.args.write_haplogroups_real_time:
                self.haplogroup_real_time_fp = self.construct_out_path(
                    type(self).haplogroup_real_time_fn_tp
                )
                self.haplogroup_real_time_file = open(  # noqa SIM115
                    self.haplogroup_real_time_fp, "w", 1
                )

            if self.args.haplogroup_to_list_genotypes_for:
                self.hg_genos_fp = self.construct_out_path(
                    type(self).hg_genos_fn_tp
                ).format(haplogroup=self.args.haplogroup_to_list_genotypes_for)
                self.hg_genos_file = open(self.hg_genos_fp, "w", 1)  # noqa SIM115

    def construct_out_path(self, fn_tp: str) -> str:
        """Return an output file path, given a filename template."""

        file_path = os.path.join(
            self.out_dir,
            fn_tp.format(fn_label=self.out_fn_label).replace("..", "."),
        )

        return file_path

    def construct_phylo_out_path(self, fn: str) -> str:
        """Return an output file path, given a filename or filename template."""

        file_path = os.path.join(self.phylo_out_dir, fn)
        return file_path

    def log_welcome_message(self) -> None:
        """Log welcome message."""

        logger.info(
            f"\n{DASHED_LINE}   yhaplo {__version__} "
            "| Y-chromosome haplogroup caller"
        )

        if self.invoked_from_command_line:
            command = os.path.basename(sys.argv[0])
            args = " ".join(sys.argv[1:])
            logger.info(f"      Command: {command} {args}")

        if self.root_logger is not None:
            logger.info(f"      Log:     {self.log_fp}")

        logger.info(DASHED_LINE)

    def override_output_generating_args(self) -> None:
        """Turn off all auxiliary output options."""

        self.args.traverse_bf = False
        self.args.traverse_df = False
        self.args.write_tree_table = False
        self.args.write_platform_trees = False

        self.args.write_anc_der_counts = False
        self.args.write_haplogroup_paths = False
        self.args.write_haplogroup_paths_detail = False
        self.args.write_der_snps = False
        self.args.write_der_snps_detail = False
        self.args.write_anc_snps = False
        self.args.write_anc_snps_detail = False

        self.args.write_haplogroups_real_time = False
        self.args.haplogroup_to_list_genotypes_for = None

    def close_real_time_output_files(self) -> None:
        """Close optional real-time output files."""

        if self.args.write_haplogroups_real_time:
            self.haplogroup_real_time_file.close()

        if self.args.haplogroup_to_list_genotypes_for:
            self.hg_genos_file.close()
            logger.info(
                "Wrote genotypes at SNPs associated with haplogroup "
                f"{self.args.haplogroup_to_list_genotypes_for}:\n"
                f"    {self.hg_genos_fp}\n"
            )


def basename_no_suffix(file_path: str, suffix: str) -> str:
    """Return the basename of a file path, with the supplied suffix removed."""

    basename_no_suffix = os.path.basename(file_path).removesuffix(suffix)
    return basename_no_suffix


def check_example_data_availability(filepath: str) -> None:
    """Check whether example data is available.

    Raises
    ------
    FileNotFoundError
        If example data is unavailable.

    """
    if not os.path.isfile(filepath):
        raise FileNotFoundError(
            f"Example input file not available: {filepath}\n\n"
            "There are two ways to run on example data:\n\n"
            "- Clone the repo and install yhaplo as editable:\n\n"
            "      cd <path/to/repo>\n"
            "      pip install --editable .\n\n"
            "- Download fixture data from tests/fixtures/input and run:\n\n"
            "      yhaplo --input <path/to/input_file> --all_aux_output\n"
        )
