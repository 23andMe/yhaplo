"""Define Path class."""

from __future__ import annotations

from collections import deque
from collections.abc import Sequence
from typing import Optional

from yhaplo import node as node_module  # noqa F401
from yhaplo import snp as snp_module  # noqa F401


class Path:
    """Class representing a path through a haplogroup tree.

    Attributes
    ----------
    node : Node
        The next node to visit.
    der_snp_list : list[SNP]
        List of SNPs observed in the derived state.
    most_derived_snp : SNP | None
        The most derived SNP observed.
    num_ancestral : int
        The number of ancestral alleles encountered.

    """

    def __init__(
        self,
        node: "node_module.Node",
    ):
        """Instantiate Path.

        Parameters
        ----------
        node : Node
            The next node to visit.

        """
        self.node = node
        self.der_snp_list: list["snp_module.SNP"] = []
        self.most_derived_snp: Optional["snp_module.SNP"] = None
        self.num_ancestral = 0
        self.init_push_through_vars()

    def init_push_through_vars(self) -> None:
        """Initialize "push-through" variables.

        These track progress subsequent to "pushing through" a branch with
        one ancestral allele and no derived alleles.

        """
        self.node_when_pushed_through: Optional["node_module.Node"] = None
        self.most_derived_snp_when_pushed_through: Optional["snp_module.SNP"] = None
        self.num_anc_since_push_through = 0
        self.num_der_since_push_through = 0

    def set_push_through_vars(self) -> None:
        """Set memory of push-though state to current state."""

        self.node_when_pushed_through = self.node
        self.most_derived_snp_when_pushed_through = self.most_derived_snp

    def update_push_through_vars(
        self,
        num_ancestral: int,
        num_derived: int,
    ) -> None:
        """Update push-through state with data from most recent branch assessment."""

        self.num_anc_since_push_through += num_ancestral
        self.num_der_since_push_through += num_derived
        if self.num_der_since_push_through > self.num_anc_since_push_through:
            self.init_push_through_vars()

    def copy_all_attributes_other_than_node(self, other: Path) -> None:
        """Copy all attributes of another path, other than its node."""

        self.der_snp_list = list(other.der_snp_list)
        self.most_derived_snp = other.most_derived_snp
        self.num_ancestral = other.num_ancestral

        self.node_when_pushed_through = other.node_when_pushed_through
        self.most_derived_snp_when_pushed_through = (
            other.most_derived_snp_when_pushed_through
        )
        self.num_anc_since_push_through = other.num_anc_since_push_through
        self.num_der_since_push_through = other.num_der_since_push_through

    def __repr__(self) -> str:
        """Return string representation."""

        return (
            f"<{__name__}.{self.__class__.__name__}: "
            f"num_ancestral={self.num_ancestral}, num_derived={self.num_derived}, "
            f'node_string="{self.node_string}", snp_string="{self.snp_string}">'
        )

    def __str__(self) -> str:
        """Return printable string representation."""

        return (
            f"{self.num_ancestral} {self.num_derived}\n"
            f"{self.node_string}\n"
            f"{self.snp_string}"
        )

    # Properties
    # ----------------------------------------------------------------------
    @property
    def has_pushed_through(self) -> bool:
        """Whether or not this path has "pushed through".

        That is, whether or not a path has proceeded past a branch with
        one ancestral allele and no derived alleles.

        """
        has_pushed_through = self.node_when_pushed_through is not None
        return has_pushed_through

    @property
    def node_string(self) -> str:
        """String concatenation of nodes visited."""

        node_string = " ".join([node.label for node in self.node.back_trace_path()])
        return node_string

    @property
    def num_derived(self) -> int:
        """Number of derived SNPs in the list."""

        num_derived = len(self.der_snp_list)
        return num_derived

    @property
    def snp_string(self) -> str:
        """String concatenation of derived SNPs observed."""

        snp_string = " ".join([snp.label for snp in self.der_snp_list])
        return snp_string

    # Regular methods
    # ----------------------------------------------------------------------
    def better_than(self, other: Path) -> bool:
        """Evaluate whether this path is better than another."""

        better_than = (
            other is None
            or self.num_derived > other.num_derived
            or (
                self.num_derived == other.num_derived
                and self.num_ancestral < other.num_ancestral
            )
        )
        return better_than

    def fork(self, node_list: Sequence["node_module.Node"]) -> deque[Path]:
        """Fork path.

        Returns
        -------
        path_deque : deque[Path]
            Deque of paths, each of which is identical to self,
            but with a new current node.

        """
        path_deque: deque[Path] = deque()
        for node in node_list:
            path = type(self)(node)
            path.copy_all_attributes_other_than_node(self)
            path_deque.append(path)

        return path_deque

    def revert_if_pushed_through_too_far(self) -> None:
        """Revert path to its state prior to pushing through.

        Do so if the path has pushed through a branch with one ancestral allele
        and no derived alleles and, since doing so, it has encountered just
        one derived allele and a nonzero number of ancestral alleles.

        """
        if (
            self.has_pushed_through
            and self.num_anc_since_push_through > 0
            and self.num_der_since_push_through == 1
        ):
            assert isinstance(self.node_when_pushed_through, node_module.Node)
            self.node = self.node_when_pushed_through
            del self.der_snp_list[-1]
            self.most_derived_snp = self.most_derived_snp_when_pushed_through
            self.num_ancestral -= self.num_anc_since_push_through
            self.init_push_through_vars()

    def update_with_branch_assessment(
        self,
        anc_snp_list: Sequence["snp_module.SNP"],
        der_snp_list: Sequence["snp_module.SNP"],
    ) -> None:
        """Update with branch assessment.

        Extend derived SNP list.
        Set most derived SNP.
        Add number of ancestral alleles seen.
        Track whether or not path has pushed through an (anc, der) == (1, 0) branch.

        """
        num_ancestral, num_derived = len(anc_snp_list), len(der_snp_list)
        self.num_ancestral += num_ancestral
        self.der_snp_list.extend(der_snp_list)
        if der_snp_list:
            self.most_derived_snp = der_snp_list[0]

        if self.has_pushed_through:
            self.update_push_through_vars(num_ancestral, num_derived)
        elif (num_ancestral, num_derived) == (1, 0):
            self.set_push_through_vars()

    # Class methods
    # ----------------------------------------------------------------------
    @classmethod
    def create_path_deque(cls, node_list: Sequence["node_module.Node"]) -> deque[Path]:
        """Return a deque of paths, each corresponding to one node in node_list."""

        path_deque: deque[Path] = deque()
        for node in node_list:
            path_deque.append(cls(node))

        return path_deque


def post_process_path_list_and_select_best(path_list: Sequence[Path]) -> Path:
    """Post-processes each Path return the best."""

    for path in path_list:
        path.revert_if_pushed_through_too_far()

    best_path = path_list[0]
    for path in path_list[1:]:
        if path.better_than(best_path):
            best_path = path

    return best_path
