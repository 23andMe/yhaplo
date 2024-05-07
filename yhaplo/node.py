"""Define Node class."""

from __future__ import annotations

import argparse
import logging
from collections import deque
from collections.abc import Iterator
from operator import attrgetter

import numpy as np

from yhaplo import sample as sample_module  # noqa F401
from yhaplo import snp as snp_module  # noqa F401
from yhaplo import tree as tree_module  # noqa F401
from yhaplo.config import Config

logger = logging.getLogger(__name__)


class Node:
    """Class representing one node of a haplogroup tree.

    Each node represents the branch that leads to it.

    Attributes
    ----------
    parent : Node | None
        Parent node. None for the root node.
    depth : int
        Number of edges between this node and the root node.
    child_list : list[Node]
        List of descendant nodes.
    snp_list : list[SNP]
        SNPs associated with the branch incident upon this node.

    haplogroup : str
        YCC haplogroup name (e.g., "R1b1c").
    label : str
        YCC haplogroup name, including alternative names (e.g., "P/K2b2").
    hg_trunc : str
        Truncated haplogroup (e.g., "R").
    hg_snp : str
        Haplogroup with representative SNP (e.g., "R-V88").

    dropped_marker_list : list[DroppedMarker]
        List of dropped markers.
    branch_length : float | None
        Branch length.
    dfs_rank : int
        Rank in depth-first listing of all nodes.

    """

    tree: tree_module.Tree
    config: Config
    args: argparse.Namespace
    hg_snp_set: set[str]

    def __init__(
        self,
        parent: Node | None,
        tree: tree_module.Tree | None = None,
    ):
        """Instantiate Node.

        Parameters
        ----------
        parent : Node | None
            Parent node. None for the root node.
        tree : Tree | None, optional
            Required for root node and ignored otherwise.
            Tree to which the root node belongs.

        """
        self.parent = parent
        if parent is None:
            if tree is not None:
                type(self).set_class_variables(tree)
            else:
                raise ValueError("Root node requires a tree instance")

            self.depth = 0
        else:
            parent.add_child(self)
            self.depth = parent.depth + 1
            if self.depth > type(self).tree.max_depth:
                type(self).tree.max_depth = self.depth

        self.haplogroup: str = ""
        self.label: str = ""
        self.hg_trunc: str = ""
        self.hg_snp: str = ""
        self.child_list: list[Node] = []
        self.snp_list: list[snp_module.SNP] = []
        self.dropped_marker_list: list[snp_module.DroppedMarker] = []
        self.branch_length: float | None = None
        self.dfs_rank: int = 0

    # String representations
    # ----------------------------------------------------------------------
    def __repr__(self) -> str:
        """Return string representation."""

        return (
            f"<{__name__}.{self.__class__.__name__}: "
            f'label="{self.label}", hg_snp="{self.hg_snp}">'
        )

    def __str__(self) -> str:
        """Return printable string representation."""

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
    def is_root(self) -> bool:
        """Whether or not the Node is the root of the tree."""

        return self.parent is None

    @property
    def is_leaf(self) -> bool:
        """Whether or not the Node is a leaf."""

        return self.num_children == 0

    @property
    def num_children(self) -> int:
        """Return number of children."""

        return len(self.child_list)

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
    def most_highly_ranked_snp(self) -> snp_module.SNP:
        """Return the most highly ranked SNP."""

        return self.snp_list[0]

    @property
    def most_highly_ranked_dropped_marker(self) -> snp_module.DroppedMarker:
        """Return the most highly ranked dropped marker."""

        return self.dropped_marker_list[0]

    # Class methods
    # ----------------------------------------------------------------------
    @classmethod
    def set_class_variables(cls, tree: tree_module.Tree) -> None:
        """Set tree, config, and args."""

        cls.tree = tree
        cls.config = tree.config
        cls.args = tree.args
        cls.hg_snp_set = set()

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
        """Set label, haplogroup, and hg_trunc.

        Also, add tree.haplogroup_to_node entry mapping haplogroup label to node.

        """
        self.label = label
        label_list = label.split("/")

        if self.is_root:
            label = self.config.root_haplogroup
            self.haplogroup = label
            self.hg_trunc = label
            label_list = [label]
        else:
            self.haplogroup = label_list[0]
            self.hg_trunc = type(self).truncate_haplogroup_label(self.haplogroup)

        for label in label_list:
            type(self).tree.haplogroup_to_node[label] = self

    def set_branch_length(self, branch_length: float) -> None:
        """Set branch length."""

        self.branch_length = branch_length

    def set_dfs_rank(self, dfs_rank: int) -> None:
        """Set depth-first search rank."""

        self.dfs_rank = dfs_rank

    def add_snp(self, snp: snp_module.SNP) -> None:
        """Append a SNP to the SNP list."""

        self.snp_list.append(snp)

    def add_dropped_marker(
        self,
        dropped_marker: snp_module.DroppedMarker,
    ) -> None:
        """Append a dropped marker to the list."""

        self.dropped_marker_list.append(dropped_marker)

    def priority_sort_snp_list_and_set_hg_snp(self) -> None:
        """Sort SNP list and set SNP-based haplogroup.

        First, sort SNP list (or dropped marker list) by priority ranking.
        Then, set reresentative-SNP-based label: self.hg_snp.
        The standard form incudes the truncated haplogroup label
        and the label of a representative SNP, separated by a hyphen (e.g. R-V88).

        Also, add tree.haplogroup_to_node entry mapping hg_snp to node.

        """
        # Root: no markers
        if self.is_root:
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
                symbol = "*" if self.is_leaf else "+"
                self.hg_snp = self.parent.hg_snp + symbol

                # Uniquify if necessary
                original_hg_snp = self.hg_snp
                i = 0
                while self.hg_snp in type(self).hg_snp_set:
                    i += 1
                    self.hg_snp = f"{original_hg_snp}{i}"

            else:
                logger.warning(
                    "WARNING. Attempted to set star label, "
                    f"but parent.hg_snp not set yet: {self.haplogroup}\n"
                )
                self.hg_snp = self.haplogroup

        type(self).hg_snp_set.add(self.hg_snp)
        type(self).tree.haplogroup_to_node[self.hg_snp] = self

    # Queries
    # ----------------------------------------------------------------------
    def get_branch_length(
        self,
        align_tips: bool = False,
        subtree_max_depth: int | None = None,
        platform: str | None = None,
    ) -> float | None:
        """Get branch length.

        Parameters
        ----------
        align_tips : bool, optional
            When True, set internal branch lengths to one and leaf branch lengths
            in such a manner as to align the tips of the tree.
        subtree_max_depth : int | None, optional
            Maximum depth of subtree.
            Used to set leaf branch lengths when aligning tips.
            Default to maximum depth of full tree.
        platform : str | None, optional
            23andMe platform to use for computing branch length.

        Returns
        -------
        branch_length : float | None
            Branch length.

        """
        if self.branch_length is not None:
            branch_length = self.branch_length
        elif align_tips and self.is_leaf:
            subtree_max_depth = subtree_max_depth or type(self).tree.max_depth
            branch_length = subtree_max_depth - self.depth + 1
        elif align_tips:
            branch_length = 1
        elif platform:
            branch_length = np.sum(
                [snp.is_on_platform(platform) for snp in self.snp_list]
            )
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
        sample: sample_module.Sample,
    ) -> tuple[list[snp_module.SNP], list[snp_module.SNP]]:
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

    # Children methods
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

    def remove_children(self) -> None:
        """Make this node a leaf by purging children."""

        self.child_list.clear()

    # Tree traversals
    # ----------------------------------------------------------------------
    def iter_depth_first(self) -> Iterator[Node]:
        """Traverse tree depth first, pre order."""

        yield self
        for child in self.child_list:
            yield from child.iter_depth_first()

    def iter_breadth_first(self) -> Iterator[Node]:
        """Traverse tree breadth first."""

        yield self
        node_deque = deque(self.child_list)
        while node_deque:
            node = node_deque.popleft()
            node_deque.extend(node.child_list)
            yield node

    def mrca(self, other_node: Node) -> Node:
        """Return the most recent common ancestor of this node and another."""

        if self.depth > other_node.depth:
            lower_node, higher_node = self, other_node
        else:
            lower_node, higher_node = other_node, self

        # Get to same level of tree
        while lower_node.depth > higher_node.depth:
            assert lower_node.parent is not None
            lower_node = lower_node.parent

        # Move up one branch at a time
        node_a, node_b = lower_node, higher_node
        while node_a is not node_b:
            assert node_a.parent is not None
            assert node_b.parent is not None
            node_a = node_a.parent
            node_b = node_b.parent

        return node_a

    # Writing tree to file in Newick format
    # ----------------------------------------------------------------------
    def write_newick(
        self,
        newick_fp: str,
        use_hg_snp_label: bool = False,
        align_tips: bool = False,
        platform: str | None = None,
        rotate: bool = False,
    ) -> None:
        """Write Newick representation of the subtree rooted at this node.

        Parameters
        ----------
        newick_fp : str
            File path to which to write Newick representation.
        use_hg_snp_label : bool, optional
            Use SNP-based haplogroup labels rather than YCC haplogroup labels.
        align_tips : bool, optional
            When True, set branch lengths to align the tips of the tree.
        platform : str | None, optional
            23andMe platform to use for computing branch lengths.
        rotate : bool, optional
            Rotate nodes. By default, branches will be ordered top to bottom.
            Rotating nodes orders branches bottom to top, which is left to right
            when an image is rotated 90 degrees to the right.

        """
        if not type(self).config.suppress_output:
            newick = self.build_newick(
                use_hg_snp_label=use_hg_snp_label,
                align_tips=align_tips,
                platform=platform,
                rotate=rotate,
            )
            with open(newick_fp, "w") as out_file:
                out_file.write(newick + "\n")

            if align_tips:
                tree_descriptor = "aligned "
            elif platform:
                tree_descriptor = f"platform {platform} "
            else:
                tree_descriptor = ""

            label_type = "representative-SNP" if use_hg_snp_label else "YCC"

            logger.info(
                f"Wrote {tree_descriptor}tree with {label_type} labels:\n"
                f"    {newick_fp}\n"
            )

    def build_newick(
        self,
        use_hg_snp_label: bool = False,
        align_tips: bool = False,
        platform: str | None = None,
        rotate: bool = False,
    ) -> str:
        """Build Newick string for the subtree rooted at this node.

        Parameters
        ----------
        use_hg_snp_label : bool, optional
            Use SNP-based haplogroup labels rather than YCC haplogroup labels.
        align_tips : bool, optional
            When True, set branch lengths to align the tips of the tree.
        platform : str | None, optional
            23andMe platform to use for computing branch lengths.
        rotate : bool, optional
            Rotate nodes. By default, branches will be ordered top to bottom.
            Rotating nodes orders branches bottom to top, which is left to right
            when an image is rotated 90 degrees to the right.

        Returns
        -------
        newick : str
            Newick representation of the tree.

        """
        subtree_max_depth = np.max([node.depth for node in self.iter_depth_first()])
        newick = (
            self.build_newick_recursive(
                use_hg_snp_label=use_hg_snp_label,
                align_tips=align_tips,
                subtree_max_depth=subtree_max_depth,
                platform=platform,
                rotate=rotate,
            )
            + ";"
        )

        return newick

    def build_newick_recursive(
        self,
        use_hg_snp_label: bool = False,
        align_tips: bool = False,
        subtree_max_depth: int | None = None,
        platform: str | None = None,
        rotate: bool = False,
    ) -> str:
        """Build Newick string recursively for the subtree rooted at this node.

        Parameters
        ----------
        use_hg_snp_label : bool, optional
            Use SNP-based haplogroup labels rather than YCC haplogroup labels.
        align_tips : bool, optional
            When True, set branch lengths to align the tips of the tree.
        subtree_max_depth : int | None, optional
            Maximum depth of subtree.
            Used to set leaf branch lengths when aligning tips.
            Default to maximum depth of full tree.
        platform : str | None, optional
            23andMe platform to use for computing branch lengths.
        rotate : bool, optional
            Rotate nodes. By default, branches will be ordered top to bottom.
            Rotating nodes orders branches bottom to top, which is left to right
            when an image is rotated 90 degrees to the right.

        Returns
        -------
        tree_string : str
            Component of Newick representation of the tree.

        """
        subtree_max_depth = subtree_max_depth or type(self).tree.max_depth
        child_list = self.child_list if not rotate else self.child_list[::-1]
        children_string = (
            (
                "("
                + ",".join(
                    [
                        child.build_newick_recursive(
                            use_hg_snp_label=use_hg_snp_label,
                            align_tips=align_tips,
                            subtree_max_depth=subtree_max_depth,
                            platform=platform,
                            rotate=rotate,
                        )
                        for child in child_list
                    ]
                )
                + ")"
            )
            if not self.is_leaf
            else ""
        )
        branch_label = self.hg_snp if use_hg_snp_label else self.label
        branch_length = self.get_branch_length(
            align_tips=align_tips,
            subtree_max_depth=subtree_max_depth,
            platform=platform,
        )
        if align_tips:
            branch_string = f"{branch_label}:{branch_length}"
        elif branch_length is None or (self.is_leaf and branch_length == 0):
            branch_string = branch_label
        elif branch_length > 0:
            branch_string = f"{branch_label}|{branch_length}:{branch_length}"
        else:
            branch_string = ":0.5"

        tree_string = f"{children_string}{branch_string}"

        return tree_string
