"""Define Node class."""

from __future__ import annotations

import argparse
import logging
from collections import deque
from operator import attrgetter
from typing import Optional, TextIO

from yhaplo import sample as sample_module  # noqa F401
from yhaplo import snp as snp_module  # noqa F401
from yhaplo import tree as tree_module  # noqa F401
from yhaplo.config import Config

logger = logging.getLogger(__name__)


class Node:

    """Class representing one node of a haplogroup tree.

    Each node represents the branch that leads to it.

    A node knows its:
    - Parent
    - Depth
    - Children
    - Diagnostic SNPs

    """

    tree: "tree_module.Tree"
    config: Config
    args: argparse.Namespace
    hg_snp_set: set[str] = set()

    def __init__(
        self,
        parent: Optional[Node],
        tree: Optional["tree_module.Tree"] = None,
    ):
        self.parent = parent
        if parent is None:
            if tree is not None:
                type(self).set_tree_config_and_args(tree)
                self.depth = 0
            else:
                raise ValueError(
                    "A tree instance must be supplied when instantiating a root node."
                )
        else:
            parent.add_child(self)
            self.depth = parent.depth + 1
            if self.depth > type(self).tree.max_depth:
                type(self).tree.max_depth = self.depth

        self.haplogroup: str = ""  # YCC haplogroup name (e.g., "R1b1c")
        self.label: str = ""  # YCC including alt names (e.g., "P/K2b2")
        self.hg_trunc: str = ""  # Truncated haplogroup (e.g., "R")
        self.hg_snp: str = ""  # Haplogroup with representative SNP (e.g., "R-V88")
        self.child_list: list[Node] = []
        self.snp_list: list["snp_module.SNP"] = []
        self.dropped_marker_list: list["snp_module.DroppedMarker"] = []
        self.branch_length: Optional[float] = None
        self.dfs_rank: int = 0

    # String representations
    # ----------------------------------------------------------------------
    def __str__(self) -> str:
        """Return string representation."""

        return self.str_simple

    @property
    def str_simple(self) -> str:
        """Return string representation with label and representative SNP."""

        return f"{self.label:25s} {self.hg_snp}"

    @property
    def str_snp_list(self) -> str:
        """Return string representation with label and list of snps."""

        snp_string = " ".join(snp.label for snp in self.snp_list)
        str_snp_list = f"{self.label:25s} {snp_string}"

        return str_snp_list

    @property
    def str_dot_pipe_depth(self) -> str:
        """Return string representation indicating depth with dots and pipes."""

        dot_list = list("." * (self.depth))
        for i in range(0, len(dot_list), 5):
            dot_list[i] = "|"

        dots = "".join(dot_list)
        str_dot_pipe_depth = f"{dots}{self.label} {self.hg_snp}"

        return str_dot_pipe_depth

    # Other properties
    # ----------------------------------------------------------------------
    @property
    def tree_table_data(self) -> tuple[str, str, str, str, str]:
        """Return a tuple of data summarizing the node.

        Returns
        -------
        tree_table_tuple : tuple[str, str, str, str, str]
            Depth-first-search rank, YCC haplogroup label, SNP-based haplogroup,
            parent DFS rank, parent SNP-based haplogroup.

        """
        if self.parent is not None:
            parent_dfs_rank = str(self.parent.dfs_rank)
            parent_hg_snp = self.parent.hg_snp
        else:
            parent_dfs_rank = "root"
            parent_hg_snp = "root"

        tree_table_row = (
            str(self.dfs_rank),
            self.haplogroup,
            self.hg_snp,
            parent_dfs_rank,
            parent_hg_snp,
        )

        return tree_table_row

    @property
    def most_highly_ranked_snp(self) -> "snp_module.SNP":
        """Return the most highly ranked SNP."""

        return self.snp_list[0]

    @property
    def most_highly_ranked_dropped_marker(self) -> "snp_module.DroppedMarker":
        """Return the most highly ranked dropped marker."""

        return self.dropped_marker_list[0]

    # Class methods
    # ----------------------------------------------------------------------
    @classmethod
    def set_tree_config_and_args(cls, tree: "tree_module.Tree") -> None:
        """Set tree, config, and args."""

        cls.tree = tree
        cls.config = tree.config
        cls.args = tree.args

    @classmethod
    def truncate_haplogroup_label(cls, haplogroup: str) -> str:
        """Return first truncated haplogroup label.

        Truncation here means the first two to five characters of specified haplogroups
        and the first letter of others.

        """
        truncated_haplogroup_label = haplogroup[0]
        for num_chars in range(cls.config.multi_char_hg_trunc_max_len, 1, -1):
            if haplogroup[:num_chars] in cls.config.multi_char_hg_trunc_set:
                truncated_haplogroup_label = haplogroup[:num_chars]
                break

        return truncated_haplogroup_label

    # Setters and mutaters
    # ----------------------------------------------------------------------
    def set_label(self, label: str) -> None:
        """Set label, haplogroup, and hg_trunc."""

        self.label = label
        label_list = label.split("/")

        if self.is_root():
            self.haplogroup = self.hg_trunc = self.config.root_haplogroup
            type(self).tree.haplogroup_to_node[self.haplogroup] = self
        else:
            self.haplogroup = label_list[0]
            self.hg_trunc = type(self).truncate_haplogroup_label(self.haplogroup)

        for key in label_list:
            type(self).tree.haplogroup_to_node[key] = self

    def set_branch_length(self, branch_length: float) -> None:
        """Set branch length."""

        self.branch_length = branch_length

    def set_dfs_rank(self, dfs_rank: int) -> None:
        """Set depth-first search rank."""

        self.dfs_rank = dfs_rank

    def add_snp(self, snp: "snp_module.SNP") -> None:
        """Append a SNP to the SNP list."""

        self.snp_list.append(snp)

    def add_dropped_marker(
        self,
        dropped_marker: "snp_module.DroppedMarker",
    ) -> None:
        """Append a dropped marker to the list."""

        self.dropped_marker_list.append(dropped_marker)

    def priority_sort_snp_list_and_set_hg_snp(self) -> None:
        """Sort SNP list and set SNP-based haplogroup.

        First, sort SNP list (or dropped marker list) by priority ranking.
        Then, set reresentative-SNP-based label: self.hg_snp.
        The standard form incudes the truncated haplogroup label
        and the label of a representative SNP, separated by a hyphen (e.g. R-V88).

        """
        # Root: no markers
        if self.is_root():
            self.hg_snp = self.haplogroup

        # Normal case
        elif self.snp_list:
            self.snp_list = snp_module.priority_sort_marker_list(self.snp_list)
            self.hg_snp = self.most_highly_ranked_snp.hg_snp

        # Backup: use discarded marker name
        elif self.dropped_marker_list:
            self.dropped_marker_list = snp_module.priority_sort_marker_list(
                self.dropped_marker_list
            )
            marker_name = self.most_highly_ranked_dropped_marker.name
            self.hg_snp = f"{self.hg_trunc}-{marker_name}"

        # No markers to use
        else:
            if self.parent is not None and self.parent.hg_snp:
                symbol = "*" if self.is_leaf() else "+"
                self.hg_snp = self.parent.hg_snp + symbol

                # Uniquify if necessary
                if self.hg_snp in type(self).hg_snp_set:
                    i = 1
                    hg_snp_uniqe = f"{self.hg_snp}{i}"
                    while hg_snp_uniqe in type(self).hg_snp_set:
                        i += 1
                        hg_snp_uniqe = f"{self.hg_snp}{i}"

                    self.hg_snp = hg_snp_uniqe
            else:
                logger.warning(
                    "WARNING. Attempted to set star label, "
                    f"but parent.hg_snp not set yet: {self.haplogroup}\n"
                )
                self.hg_snp = self.haplogroup

        type(self).hg_snp_set.add(self.hg_snp)

    # Queries
    # ----------------------------------------------------------------------
    def is_root(self) -> bool:
        """Return a Boolean indicating whether or not the Node is root."""

        return self.parent is None

    def is_leaf(self) -> bool:
        """Return a Boolean indicating whether or not the Node is a leaf."""

        return len(self.child_list) == 0

    def get_branch_length(
        self,
        align_tips: bool = False,
        platform: Optional[str] = None,
    ) -> Optional[float]:
        """Get branch length."""

        if self.branch_length:
            branch_length = self.branch_length
        elif align_tips and self.is_leaf():
            branch_length = type(self).tree.max_depth - self.depth + 1
        elif align_tips:
            branch_length = 1
        elif platform:
            branch_length = 0
            for snp in self.snp_list:
                if snp.is_on_platform(platform):
                    branch_length += 1
        else:
            branch_length = None

        return branch_length

    def back_trace_path(self) -> list[Node]:
        """Return a list of nodes from root to self."""

        node_list = [self]
        parent = self.parent
        while parent is not None:
            node_list.append(parent)
            parent = parent.parent

        node_list.reverse()

        return node_list

    def assess_genotypes(
        self,
        sample: "sample_module.Sample",
    ) -> tuple[list["snp_module.SNP"], list["snp_module.SNP"]]:
        """Assess an individual's genotypes with respect to self.snp_list.

        Returns
        -------
        anc_snp_list : list[SNP]
            SNPs for which ancestral genotypes were observed.
        der_snp_list : list[SNP]
            SNPs for which derived genotypes were observed.

        """
        anc_snp_list, der_snp_list = [], []
        for snp in self.snp_list:
            genotype = sample.get_genotype(snp.position)
            if genotype != Config.missing_genotype:
                if snp.is_ancestral(genotype):
                    anc_snp_list.append(snp)
                elif snp.is_derived(genotype):
                    der_snp_list.append(snp)

                if type(self).args.haplogroup_to_list_genotypes_for == self.haplogroup:
                    file = type(self).config.hg_genos_file
                    derived_flag = "*" if snp.is_derived(genotype) else ""
                    file.write(
                        f"{str(sample.iid):8s} {snp} {genotype} {derived_flag}\n"
                    )

        return anc_snp_list, der_snp_list

    # Children
    # ----------------------------------------------------------------------
    def add_child(self, child: Node) -> None:
        """Append a child to the child list."""

        self.child_list.append(child)

    def serial_split(self, target_haplogroup: str) -> Node:
        """Split node serially until there is a spot for the target haplogroup."""

        current_node = self
        start_length = len(self.haplogroup)
        end_length = len(target_haplogroup)
        for str_len in range(start_length, end_length):
            next_node = None
            target_hg_substring = target_haplogroup[: (str_len + 1)]
            if current_node.num_children < 2:
                current_node.bifurcate()
            for node in current_node.child_list:
                if node.haplogroup == target_hg_substring:
                    next_node = node
            if next_node is None:
                next_node = type(self)(parent=current_node)
                next_node.set_label(target_hg_substring)
                current_node.sort_children()

            current_node = next_node

        return current_node

    @property
    def num_children(self) -> int:
        """Return number of children."""

        return len(self.child_list)

    def bifurcate(self) -> tuple[Node, Node]:
        """Split a node and return the two children."""

        left_child = type(self)(parent=self)
        right_child = type(self)(parent=self)
        if self.haplogroup[-1].isalpha():
            left_child.set_label(self.haplogroup + "1")
            right_child.set_label(self.haplogroup + "2")
        else:
            left_child.set_label(self.haplogroup + "a")
            right_child.set_label(self.haplogroup + "b")

        return left_child, right_child

    def sort_children(self) -> None:
        """Sort children by haplogroup."""

        self.child_list = sorted(self.child_list, key=attrgetter("haplogroup"))

    def reverse_children(self) -> None:
        """Reverse child list."""

        self.child_list.reverse()

    # Tree traversals
    # ----------------------------------------------------------------------
    def write_breadth_first_traversal(self, bf_tree_file: TextIO) -> None:
        """Write breadth-first traversal."""

        bf_tree_file.write(self.str_dot_pipe_depth + "\n")
        node_deque = deque(self.child_list)
        while node_deque:
            node = node_deque.popleft()
            bf_tree_file.write(node.str_dot_pipe_depth + "\n")
            node_deque.extend(node.child_list)

    def get_depth_first_node_list(self) -> list[Node]:
        """Conduct depth-first pre-order traversal."""

        depth_first_node_list = [self]
        self.traverse_depth_first_pre_order_recursive(depth_first_node_list)

        return depth_first_node_list

    def traverse_depth_first_pre_order_recursive(
        self,
        depth_first_node_list: list[Node],
    ) -> None:
        """Append each node in depth-first pre order, recursively."""

        for child in self.child_list:
            depth_first_node_list.append(child)
            child.traverse_depth_first_pre_order_recursive(depth_first_node_list)

    def mrca(self, other_node: Node) -> Node:
        """Return the most recent common ancestor of this node and another."""

        if self.depth < other_node.depth:
            higher_node, lower_node = self, other_node
        else:
            higher_node, lower_node = other_node, self

        while higher_node.depth < lower_node.depth:
            assert lower_node.parent is not None
            lower_node = lower_node.parent
        while lower_node != higher_node:
            assert lower_node.parent is not None
            assert higher_node.parent is not None
            lower_node = lower_node.parent
            higher_node = higher_node.parent

        return higher_node

    # Writing tree to file in Newick format
    # ----------------------------------------------------------------------
    def write_newick(
        self,
        newick_fp: str,
        use_hg_snp_label: bool = False,
        align_tips: bool = False,
        platform: Optional[str] = None,
    ) -> None:
        """Write Newick string for the subtree rooted at this node."""

        if not type(self).config.suppress_output:
            with open(newick_fp, "w") as out_file:
                out_file.write(
                    self.build_newick_string_recursive(
                        use_hg_snp_label,
                        align_tips,
                        platform,
                    )
                    + ";\n"
                )

            if align_tips:
                tree_descriptor = "aligned "
            elif platform:
                tree_descriptor = f"platform {platform} "
            else:
                tree_descriptor = ""

            if use_hg_snp_label:
                label_type = "representative-SNP"
            else:
                label_type = "YCC"

            logger.info(
                f"Wrote {tree_descriptor}tree with {label_type} labels:\n"
                f"    {newick_fp}\n"
            )

    def build_newick_string_recursive(
        self,
        use_hg_snplabel: bool = False,
        align_tips: bool = False,
        platform: Optional[str] = None,
    ) -> str:
        """Build Newick string recursively for the subtree rooted at this node."""

        if not self.is_leaf():
            child_string_list = []
            for child in self.child_list[::-1]:
                child_string = child.build_newick_string_recursive(
                    use_hg_snplabel,
                    align_tips,
                    platform,
                )
                child_string_list.append(child_string)

            children = ",".join(child_string_list)
            tree_string_part_1 = f"({children})"
        else:
            tree_string_part_1 = ""

        branch_label = self.hg_snp if use_hg_snplabel else self.label
        branch_length = self.get_branch_length(align_tips, platform)
        if align_tips:
            branch_string = f"{branch_label}:{branch_length}"
        elif branch_length is None or (self.is_leaf() and branch_length == 0):
            branch_string = branch_label
        elif branch_length > 0:
            branch_string = f"{branch_label}|{branch_length}:{branch_length}"
        else:
            branch_string = ":0.5"

        tree_string = f"{tree_string_part_1}{branch_string}"

        return tree_string
