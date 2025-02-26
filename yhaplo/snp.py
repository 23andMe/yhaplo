"""Define SNP and related classes.

Classes defined herein include:
* SNP
* DroppedMarker

"""

from __future__ import annotations

import itertools
import logging
import re
from collections.abc import Mapping, Sequence
from operator import attrgetter
from typing import TypeVar

import pandas as pd
import yaml

from yhaplo import node as node_module  # noqa F401
from yhaplo import tree as tree_module  # noqa F401
from yhaplo.config import Config
from yhaplo.utils.loaders import (
    DataFile,
    load_data,
    load_data_lines,
    load_dataframe,
)

logger = logging.getLogger(__name__)


class SNP:
    """Class representing a SNP.

    Attributes
    ----------
    name_list : list[str]
        SNP names.
    haplogroup : str
        Haplogroup the SNP is associated with.
    position : int
        GRCh37 physical position.
    ancestral : str
        Ancestral allele.
    derived : str
        Derived allele.

    is_representative : bool
        Indicator as to whether this SNP should be used to represent the haplogroup.
    allele_set : set[str]
        Set of ancestral and derived alleles.
    node : Node
        Node associated with the haplogroup the SNP is associated with.

    """

    tree: tree_module.Tree
    pos_to_block_indexes: dict[int, list[int]] = {}
    platform_to_pos_set: dict[str, set[int]] = {}

    def __init__(
        self,
        name: str,
        haplogroup: str,
        position: int,
        ancestral: str,
        derived: str,
    ):
        """Instantiate SNP.

        Parameters
        ----------
        name : str
            SNP name.
        haplogroup : str
            Haplogroup the SNP is associated with.
        position : int
            Physical position.
        ancestral : str
            Ancestral allele.
        derived : str
            Derived allele.

        """
        if type(self).tree is None:
            raise RuntimeError(
                "Before instantiating, call: "
                f"{self.__class__.__name__}.set_class_variables(tree)"
            )

        self.set_label(name)
        self.name_list = [name]
        self.is_representative = name in type(self).tree.representative_snp_name_set

        self.haplogroup = haplogroup
        self.position = position
        self.ancestral = ancestral
        self.derived = derived
        self.allele_set = {ancestral, derived}

        self.node = type(self).tree.find_or_create_node(haplogroup)
        self.node.add_snp(self)

    def set_label(self, label: str) -> None:
        """Set label and associated instance variables."""

        self.label = label
        (
            self.label_letters_rank,
            self.label_letters,
            self.label_number,
        ) = parse_snp_label(label, Config.snp_label_letters_rank_dict)
        self.label_cleaned = clean_snp_label(label)

    def __repr__(self) -> str:
        """Return string representation."""

        return (
            f"<{__name__}.{self.__class__.__name__}: "
            f'label="{self.label}", node.label="{self.node.label}", '
            f"position={self.position}, "
            f'mutation="{self.ancestral}->{self.derived}">'
        )

    def __str__(self) -> str:
        """Return printable string representation."""

        return (
            f"{self.label:15s} {self.node.label:25s} {self.position:8d} "
            f"{self.ancestral}->{self.derived}"
        )

    @property
    def info(self) -> str:
        """Return multiline summary of SNP."""

        names = [name for name in self.name_list if name != self.label]
        aliases_str = ", ".join(names) if names else "None"
        info = (
            f"Name: {self.label}\n"
            f"YCC haplogroup: {self.node.label}\n"
            f"GRCh37 position: {self.position:,}\n"
            f"Mutation: {self.ancestral}->{self.derived}\n"
            f"Aliases: {aliases_str}"
        )
        return info

    @property
    def str_with_all_names(self) -> str:
        """Return long string representation.

        This includes the normal string representation,
        plus a comma-separated list of names.

        """
        names_str = ",".join(self.name_list)
        str_with_all_names = f"{str(self)}     {names_str}"
        return str_with_all_names

    @property
    def str_short(self) -> str:
        """Return short string representation: Node label and SNP label."""

        str_short = f"{self.node.label}:{self.label}"
        return str_short

    @property
    def dfs_rank(self) -> int:
        """Return depth-first search rank."""

        return self.node.dfs_rank

    @property
    def hg_snp(self) -> str:
        """Return string representation with truncated haplogroup label and SNP label.

        Example: R-V88

        """
        hg_snp = f"{self.node.hg_trunc}-{self.label_cleaned}"
        return hg_snp

    def is_derived(self, geno: str) -> bool:
        """Return True if geno is the derived allele."""

        is_derived = geno == self.derived
        return is_derived

    def is_ancestral(self, geno: str) -> bool:
        """Return True if geno is the ancestral allele."""

        is_ancestral = geno == self.ancestral
        return is_ancestral

    def is_on_platform(self, platform: str) -> bool:
        """Return True if this SNP is on the supplied 23andMe platform."""

        is_on_platform = self.position in type(self).platform_to_pos_set[platform]
        return is_on_platform

    def back_trace_path(self) -> list[node_module.Node]:
        """Return the backtrace path (node list) for the corresponding node."""

        back_trace_path = self.node.back_trace_path()
        return back_trace_path

    def add_name(self, name: str) -> None:
        """Add name to list and update label if appropriate."""

        self.name_list.append(name)
        if name in type(self).tree.representative_snp_name_set:
            self.is_representative = True

        if type(self).is_a_preferred_name(self.label):
            if type(self).is_a_preferred_name(name):
                logger.warning(
                    f"WARNING. Two preferred names for one SNP: {name}, {self.label}\n"
                )
        elif type(self).is_a_preferred_name(name):
            self.set_label(name)
        else:
            label_letters_rank, label_letters, label_number = parse_snp_label(
                name,
                Config.snp_label_letters_rank_dict,
            )
            if label_letters_rank < self.label_letters_rank or (
                label_letters == self.label_letters and label_number < self.label_number
            ):
                self.set_label(name)

    @classmethod
    def set_class_variables(
        cls,
        tree: tree_module.Tree,
    ) -> None:
        """Set class variables.

        Doing so enables the SNP class to know about the tree instance,
        configuration and command-line arguments.

        """
        cls.tree = tree
        if tree.config.run_from_ablocks or tree.args.write_platform_trees:
            cls.pos_to_block_indexes = load_pos_to_block_indexes()
            cls.platform_to_pos_set = build_platform_to_pos_set()

    @classmethod
    def is_a_preferred_name(cls, name: str) -> bool:
        """Return True if a SNP name is in the set of preferred names.

        Ignore extensions like ".1".

        """
        is_a_preferred_name = (
            name in cls.tree.preferred_snp_name_set
            or name.split(".")[0] in cls.tree.preferred_snp_name_set
        )
        return is_a_preferred_name


def clean_snp_label(label: str) -> str:
    """Remove superfluous text and hyphens from a SNP label."""

    for superfluous_snp_text in Config.superfluous_snp_text_list:
        label = label.replace(superfluous_snp_text, "")

    label = label.replace("-", "_").replace("^", "").replace("≤", "<=")

    return label


def parse_snp_label(
    name: str,
    snp_label_letters_rank_dict: Mapping[str, int],
) -> tuple[int, str, int]:
    """Parse SNP label.

    Returns
    -------
    label_letters_rank
        Priority rank of SNP name.
    label_letters
        SNP-name letters.
    label_number
        SNP-name numbers.

    """
    match = re.search(r"([a-zA-Z-]*)([0-9]*)", str(name))
    if match is None:
        raise ValueError(f"SNP name unparsable: {name}")

    label_letters, label_number = match.group(1), match.group(2)
    label_number = int(label_number) if len(label_number) > 0 else 0
    if label_letters in snp_label_letters_rank_dict:
        label_letters_rank = snp_label_letters_rank_dict[label_letters]
    else:
        label_letters_rank = len(snp_label_letters_rank_dict)  # max value

    return label_letters_rank, label_letters, label_number


def load_pos_to_block_indexes() -> dict[int, list[int]]:
    """Load mapping from physical position to block indexes."""

    pos_to_block_indexes = yaml.safe_load(
        load_data(Config.pos_to_block_indexes_data_file)
    )
    num_positions = len(pos_to_block_indexes)
    num_block_indexes = len(
        list(itertools.chain.from_iterable(pos_to_block_indexes.values()))
    )
    logger.info(
        "Loaded mapping from physical position to block indexes\n"
        f"{num_positions:7d} positions: "
        f"{Config.pos_to_block_indexes_data_file.filename}\n"
        f"{num_block_indexes:7d} ablock indexes\n"
    )

    return pos_to_block_indexes


def build_platform_to_pos_set() -> dict[str, set[int]]:
    """Build mapping from platform to a set of physical positions."""

    logger.info("Loading platform positions...")

    platform_to_pos_set = {}
    for platform in Config.platforms:
        pos_set = load_platform_position_set(platform)
        exclude_pos_set = load_platform_exclude_position_set(platform)
        pos_set -= exclude_pos_set
        platform_to_pos_set[platform] = pos_set

        if exclude_pos_set:
            logger.info(f"{platform}: {len(pos_set):5d} unique positions remain")

    logger.info("")

    return platform_to_pos_set


def load_platform_position_set(
    platform: str,
    log_loader: bool = False,
) -> set[int]:
    """Load physical positions of Y-chromosome SNPs on a 23andMe platform.

    Parameters
    ----------
    platform : str
        Platform name.
    log_loader : bool, optional
        Indicator to log loader activity.

    Returns
    -------
    pos_set : set[int]
        Set of physical positions.

    """
    platform_pos_data_file = DataFile(
        Config.platform_pos_data_subdir,
        Config.platform_pos_fn_tp.format(platform=platform),
        f"Platform {platform} SNP positions",
        ttam_only=True,
    )
    pos_set = {
        int(line.strip().split()[0])
        for line in load_data_lines(platform_pos_data_file, log=log_loader)
    }
    logger.info(
        f"{platform}: {len(pos_set):5d} unique positions loaded: "
        f"{platform_pos_data_file.filename}"
    )

    return pos_set


def load_platform_exclude_position_set(platform: str) -> set[int]:
    """Load physical positions of SNPs to exclude from a 23andMe platform.

    Parameters
    ----------
    platform : str
        Platform name.

    Returns
    -------
    exclude_pos_set : set[int]
        Set of physical positions to exclude.

    """
    try:
        exclude_df = load_platform_exclude_table(platform)
        exclude_pos_set = set(exclude_df["position"])
        num_exclude_pos = len(exclude_pos_set)
        logger.info(f"{platform}: {num_exclude_pos:5d} unique positions to exclude")
    except FileNotFoundError:
        exclude_pos_set = set()

    return exclude_pos_set


def load_platform_exclude_table(platform: str) -> pd.DataFrame:
    """Load table of SNPs to exclude from a 23andMe platform.

    Parameters
    ----------
    platform : str
        Platform name.

    Returns
    -------
    exclude_df : pd.DataFrame
        Table of SNPs to exclude.

    """
    platform_pos_exclude_data_file = DataFile(
        Config.platform_pos_data_subdir,
        Config.platform_qc_exclude_fn_tp.format(platform=platform),
        f"Platform {platform} SNP QC exclusions table",
        ttam_only=True,
    )
    exclude_df = load_dataframe(platform_pos_exclude_data_file)
    logger.info(
        f"{platform}: {len(exclude_df):5d} QC exclusions loaded: "
        f"{platform_pos_exclude_data_file.filename}"
    )

    return exclude_df


class DroppedMarker:
    """Class representing a marker not used for classification.

    Such a marker may be useful for node labeling. Examples:
    - Non-SNPs
    - Multiallelic SNPs
    - SNPs not meeting ISOGG quality guidelines

    Parameters
    ----------
    name : str
        Marker name.
    haplogroup : str
        Haplogroup the marker is associated with.
    tree : Tree
        Y-chromosome tree.

    """

    def __init__(
        self,
        name: str,
        haplogroup: str,
        tree: tree_module.Tree,
    ):
        self.name = clean_snp_label(name)
        self.haplogroup = haplogroup
        self.tree = tree

    def __repr__(self) -> str:
        """Return string representation."""

        return (
            f"<{__name__}.{self.__class__.__name__}: "
            f'name="{self.name}", haplogroup="{self.haplogroup}">'
        )

    def add_to_node(self) -> bool:
        """Add this dropped marker to the corresponding node, if it exists."""

        added = False
        if self.haplogroup in self.tree.haplogroup_to_node:
            (
                self.label_letters_rank,
                self.label_letters,
                self.label_number,
            ) = parse_snp_label(self.name, Config.snp_label_letters_rank_dict)
            self.is_representative = self.name in self.tree.representative_snp_name_set
            self.tree.haplogroup_to_node[self.haplogroup].add_dropped_marker(self)
            added = True

        return added


Marker = TypeVar("Marker", SNP, DroppedMarker)


def priority_sort_marker_list(marker_list: Sequence[Marker]) -> list[Marker]:
    """Sort a list of markers by priority ranking.

    Preference is given to those deemed representative for the corresponding haplogroup.

    """
    sorted_marker_list = sorted(
        sorted(
            marker_list,
            key=attrgetter(
                "label_letters_rank",
                "label_letters",
                "label_number",
            ),
        ),
        key=attrgetter("is_representative"),
        reverse=True,
    )

    return sorted_marker_list
