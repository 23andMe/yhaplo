"""Define Tree class."""

import logging
import re
from collections import defaultdict, deque
from operator import attrgetter
from typing import Optional, Union

import pandas as pd

from yhaplo import node as node_module  # noqa F401
from yhaplo import path as path_module  # noqa F401
from yhaplo import sample as sample_module  # noqa F401
from yhaplo import snp as snp_module  # noqa F401
from yhaplo.config import Config
from yhaplo.utils.loaders import load_data, load_data_lines, load_dataframe

logger = logging.getLogger(__name__)


class Tree:
    """Class representing a haplogroup tree.

    Attributes
    ----------
    root : Node
        Root node.
    config : Config
        Yhaplo configuration.
    args : argparse.Namespace
        Command-line arguments or defaults.
    haplogroup_to_node : dict[str, Node]
        Maps haplogroup labels (YCC or SNP-based) to nodes.
    depth_first_node_list : list[Node]
        List of nodes in depth-first order.
    max_depth : int
        Maximum depth of tree.

    snp_dict : dict[Union[str, tuple[str, int], int], SNP]
        Maps SNP names, positions, and (haplogroup, position) tuples to SNPs.
    snp_list : list[SNP]
        List of SNPs
    snp_pos_set : set[int]
        Set of SNP physical positions.
    snp_name_set : set[str]
        Set of SNP names.
    preferred_snp_name_set : set[str]
        Set of preferred SNP names.
    representative_snp_name_set : set[str]
        Set of names of "representative" SNPs.

    isogg_omit_name_set : set[str]
        Set of SNP names to omit from ISOGG variant metadata.
    isogg_omit_pos_str_mutation_set : set[tuple[str, str]]
        Set of (position_str, mutation) tuples representing SNPs to omit.
    isogg_correction_dict : dict[str, tuple[str, str, str]]
        Maps SNP name to a tuple of haplogroup, position, and mutation.
        Either position or mutation (or both) are corrections from the ISOGG table.
    isogg_counts_dict : dict[str, int]
        Counts of various categories of ISOGG SNPs.
    num_snps_corrected : int
        Number of SNPs corrected, as compared to ISOGG raw data.
    multiallelic_old_pos_set : set[int]
        Set of positions to exclude due to multiple alleles.
    multiallelic_new_pos_set : set[int]
        Set of positions found to have multiple alleles.

    """

    def __init__(
        self,
        config: Optional[Config] = None,
    ):
        """Instantiate Tree.

        Parameters
        ----------
        config : Config | None, optional
            Configuration. If None, use default configuration.

        """
        self.config = config if config is not None else Config()
        self.args = self.config.args

        self.max_depth = 0
        self.haplogroup_to_node: dict[str, "node_module.Node"] = {}
        self.depth_first_node_list: list["node_module.Node"] = []

        self.snp_dict: dict[
            Union[str, tuple[str, int], int],
            "snp_module.SNP",
        ] = {}  # Possible keys: name, (haplogroup, position), position
        self.snp_list: list["snp_module.SNP"] = []
        self.snp_pos_set: set[int] = set()
        self.snp_name_set: set[str] = set()
        self.preferred_snp_name_set: set[str] = set()
        self.representative_snp_name_set: set[str] = set()
        self.multiallelic_old_pos_set: set[int] = set()
        self.multiallelic_new_pos_set: set[int] = set()
        self.isogg_omit_name_set: set[str] = set()
        self.isogg_omit_pos_str_mutation_set: set[tuple[str, str]] = set()
        self.isogg_correction_dict: dict[str, tuple[str, str, str]] = {}
        self.isogg_counts_dict: dict[str, int] = defaultdict(int)
        self.num_snps_corrected = 0

        self.root = self.build_tree_from_newick()
        if self.args.primary_only:
            self.set_depth_first_node_list()
        else:
            self.import_isogg_snps()

        self.set_search_root()
        self.write_optional_traversal_output()

    def __repr__(self) -> str:
        """Return string representation."""

        return (
            f"<{__name__}.{self.__class__.__name__}: "
            f"{len(self.depth_first_node_list)} nodes, {len(self.snp_list)} SNPs, "
            f"max_depth={self.max_depth}>"
        )

    # Setters
    # ----------------------------------------------------------------------
    def set_search_root(self) -> None:
        """Set node from which to start haplogroup-calling traversals."""

        if self.args.alternative_root:
            alternative_root_hg = self.args.alternative_root
            if alternative_root_hg in self.haplogroup_to_node:
                self.search_root = self.haplogroup_to_node[alternative_root_hg]
                logger.info(
                    "Will start haplogroup assignment traversal from:\n"
                    f"    {alternative_root_hg}\n\n"
                )
            else:
                raise ValueError(
                    "Cannot start traversal "
                    f"from non-existant haplogroup: {alternative_root_hg}\n"
                )
        else:
            self.search_root = self.root

    def set_depth_first_node_list(self) -> None:
        """Build Node list from depth-first pre-order traversal."""

        self.depth_first_node_list = [node for node in self.root.iter_depth_first()]
        for dfs_rank, node in enumerate(self.depth_first_node_list):
            node.set_dfs_rank(dfs_rank)

    # Traversals
    # ----------------------------------------------------------------------
    def write_optional_traversal_output(self) -> None:
        """Write optional tree-traversal output."""

        if self.args.traverse_bf:
            self.write_breadth_first()
        if self.args.traverse_df:
            self.write_depth_first_pre_order()
        if self.args.write_tree_table:
            self.write_tree_table()

        if self.args.mrca_haplogroup_list:
            self.query_mrca(*self.args.mrca_haplogroup_list)
        if self.args.query_haplogroups:
            for haplogroup in self.args.query_haplogroups.split(","):
                self.query_haplogroup(haplogroup)
        if self.args.query_snp_names:
            for snp_name in self.args.query_snp_names.split(","):
                self.query_snp_path(snp_name)

    def write_breadth_first(self) -> None:
        """Write bread-first traversal in pipe/dot format."""

        bf_tree_fp = (
            self.config.bf_primary_tree_fp
            if self.args.primary_only
            else self.config.bf_tree_fp
        )
        with open(bf_tree_fp, "w") as bf_tree_file:
            for node in self.root.iter_breadth_first():
                bf_tree_file.write(f"{node.str_dot_pipe_depth}\n")

        logger.info(f"Wrote breadth-first tree traveral:\n    {bf_tree_fp}\n")

    def write_depth_first_pre_order(self) -> None:
        """Write depth-first pre-order traversal in pipe/dot format."""

        df_tree_fp = (
            self.config.df_primary_tree_fp
            if self.args.primary_only
            else self.config.df_tree_fp
        )
        with open(df_tree_fp, "w") as df_tree_file:
            for node in self.depth_first_node_list:
                df_tree_file.write(f"{node.str_dot_pipe_depth}\n")

        logger.info(f"Wrote depth-first tree traveral:\n    {df_tree_fp}\n")

    def write_tree_table(self) -> None:
        """Write depth-first pre-order traversal in table format."""

        tree_df = self.generate_tree_table()
        tree_df.to_csv(self.config.tree_table_fp, sep="\t", index=False)
        logger.info(f"Wrote tree table:\n    {self.config.tree_table_fp}\n")

    def generate_tree_table(self) -> pd.DataFrame:
        """Generate tree table from depth-first pre-order traversal."""

        tree_df = pd.DataFrame(
            [node.tree_table_data for node in self.depth_first_node_list],
            columns="#index ycc_label label parent_index parent_label".split(),
        ).astype("string")

        return tree_df

    def query_mrca(self, haplogroup1: str, haplogroup2: str) -> None:
        """Write MRCA of two haplogroups.

        Parameters
        ----------
        haplogroup1 : str
            Haplogroup label.
        haplogroup2 : str
            Haplogroup label.

        """
        node1 = self.haplogroup_to_node[haplogroup1]
        node2 = self.haplogroup_to_node[haplogroup2]
        mrca = node1.mrca(node2)
        logger.info(
            "\nMRCA Query\n\n"
            f"Haplogroup 1: {node1.haplogroup} ({node1.hg_snp})\n"
            f"Haplogroup 2: {node2.haplogroup} ({node2.hg_snp})\n"
            f"MRCA: {mrca.haplogroup} ({mrca.hg_snp})\n"
        )

    def query_haplogroup(self, haplogroup: str) -> None:
        """List phylogenetic path for a query haplogroup.

        Parameters
        ----------
        haplogroup : str
            Name of haplogroup to query.

        """
        logger.info(f"\nHaplogroup Query: {haplogroup}\n")
        haplogroup_node = self.haplogroup_to_node.get(haplogroup)

        if haplogroup_node:
            logger.info("Descent\n------")
            for node in haplogroup_node.back_trace_path():
                logger.info(node.str_simple)

            logger.info("\nSNPs\n----")
            for snp in haplogroup_node.snp_list:
                logger.info(snp)

        else:
            logger.info(f'Haplogroup "{haplogroup}" not found')

        logger.info("")

    def query_snp_path(self, snp_name: str) -> None:
        """List phylogenetic path for a query SNP.

        Parameters
        ----------
        snp_name : str
            Name of SNP to query.

        """
        logger.info(f"\nSNP Query: {snp_name}\n")
        snp = self.snp_dict.get(snp_name)

        if snp:
            logger.info(f"{snp.info}\n")
            for node in snp.back_trace_path():
                logger.info(node.str_simple)

            if snp_name != snp.label:
                logger.info(f"\nNote: {snp_name} is an alias of {snp.label}.")
            elif not node.hg_snp.endswith(snp_name):
                logger.info(f"\nNote: {snp_name} is on the {node.hg_snp} branch.")
        else:
            logger.info(f'No SNP found with the name "{snp_name}"')

        logger.info("")

    # Write Newick files
    # ----------------------------------------------------------------------
    def write_newick(self) -> None:
        """Write tree as-is and with aligned terminal branch lengths."""

        if not self.config.suppress_output:
            logger.info("\nTree output\n")
            self.root.write_newick(self.config.ycc_tree_fp)
            self.root.write_newick(self.config.hg_snp_tree_fp, use_hg_snp_label=True)
            self.root.write_newick(self.config.aligned_ycc_tree_fp, align_tips=True)
            self.root.write_newick(
                self.config.aligned_hg_snp_tree_fp,
                use_hg_snp_label=True,
                align_tips=True,
            )
            if self.args.write_platform_trees:
                self.write_platform_trees()

    def write_platform_trees(self) -> None:
        """Write trees whose branch lengths are numbers of platform sites."""

        for platform in Config.platforms:
            self.root.write_newick(
                self.config.platform_ycc_tree_fp_tp.format(platform=platform),
                platform=platform,
            )
            self.root.write_newick(
                self.config.platform_hg_snp_tree_fp_tp.format(platform=platform),
                use_hg_snp_label=True,
                platform=platform,
            )

    # Query
    # ----------------------------------------------------------------------
    def identify_phylogenetic_path(
        self,
        sample: "sample_module.Sample",
    ) -> tuple[
        "path_module.Path",
        list["snp_module.SNP"],
        list[tuple["node_module.Node", int, int]],
    ]:
        """Identify phylogenetic path for haplogroup call.

        Conduct a modified breadth-first search (BFS) to identify
        the phylogenetic path leading from the root to the most terminal branch
        representing a Sample's haplogroup.

        Parameters
        ----------
        sample : Sample
            Sample instance.

        Returns
        -------
        best_path : Path
            The best phylogenetic path.
        anc_snp_full_list : list[SNP]
            List of SNPs observed in the ancestral state.
        anc_der_count_tuples : list[tuple[Node, int, int]]
            List of (node, num_ancestral, num_derived) tuples.

        Notes
        -----
        The key differences from a standard BFS are:
        - Stopping condition is robust to genotype error, homoplasy, etc.
        - Collapsing condition to speed up and (marginally) improve accuracy

        When the stopping condition is met, add the current path to a list.
        At the end, post-processes this list and select the best element.

        The stopping condition is a disjunction of three atomic conditions.
        The first is trivial:

        a. node.is_leaf
           We cannot go any further.

        The following table enumerates possible cases for the other two
        atomic conditions.

        #Anc: Number of ancestral alleles observed on a branch.
        #Der: Number of derived alleles observed on the branch.
              These are only considered if #anc == 2.
        Stop: Whether or not to stop.

        | #Anc | #Der | Stop | Reason
        |------|------|------|--------------------------------------------------------
        | 0, 1 |    . |   no | Insufficient evidence to stop
        |    2 |   1+ |   no | Given evidence to continue, do so for robustness
        |    2 |    0 |  yes | Reasonable evidence to stop and no evidence to continue
        |   3+ |    . |  yes | Strong evidence to stop

        b. Row 4: num_ancestral >= 3
           num_derived == 0: Compelling evidence to stop.
           num_derived >= 1: The sample's lineage probably diverges from the known tree here.

        c. Row 3: num_ancestral == 2 and num_derived == 0
           It is safe to assume that this path will not yield fruit.

        These conditions are robust to the most challenging case:
        when just a single SNP is genotyped on a branch, and the observed genotype
        corresponds to the ancestral allele due to genotype error, homoplasy,
        or an uncorrected ISOGG error. When at least one derived allele is observed,
        the conditions are also robust to two false ancestral alleles on a branch.

        """
        path_deque = path_module.Path.create_path_deque(self.search_root.child_list)
        stopped_path_list = []
        anc_snp_full_list = []
        anc_der_count_tuples = []
        while path_deque:
            path = path_deque.popleft()
            anc_snp_list, der_snp_list = path.node.assess_genotypes(sample)
            path.update_with_branch_assessment(anc_snp_list, der_snp_list)
            anc_snp_full_list.extend(anc_snp_list)
            num_ancestral, num_derived = len(anc_snp_list), len(der_snp_list)
            anc_der_count_tuples.append((path.node, num_ancestral, num_derived))

            if (
                path.node.is_leaf
                or (num_ancestral > self.config.args.anc_stop_thresh)
                or (
                    num_ancestral == self.config.args.anc_stop_thresh
                    and num_derived == 0
                )
            ):
                stopped_path_list.append(path)
            else:
                if num_derived >= self.config.args.der_collapse_thresh:
                    path_deque = deque()

                path_deque.extend(path.fork(path.node.child_list))

        best_path = path_module.post_process_path_list_and_select_best(
            stopped_path_list
        )

        return best_path, anc_snp_full_list, anc_der_count_tuples

    def get_dfs_rank(self, haplogroup: str) -> int:
        """Return the DFS rank of a haplogroup."""

        dfs_rank = self.haplogroup_to_node[haplogroup].dfs_rank
        return dfs_rank

    # Build tree from Newick-formatted text file
    # ----------------------------------------------------------------------
    def build_tree_from_newick(self) -> "node_module.Node":
        """Read a Newick-formatted tree and build a Tree instance.

        Discard bootstrap values.

        Returns
        -------
        root : Node
            Root of the tree.

        """
        logger.info("\nPrimary tree\n")
        tree_string = load_data(self.config.primary_tree_data_file, log=True).strip()

        # Tokenization:
        # a. Strip out bootstraps: text within brackets.
        # b. Split on any semantic token.
        # c. Group to retain retain tokens themselves.
        # d. Drop empty tokens from splitting adjacent semantic tokens.
        tree_string = re.subn(r"\[.*?\]", "", tree_string)[0]
        tree_list = re.split(
            f"([{self.config.newick_semantic_token_string}])",
            tree_string,
        )
        tree_list = [token for token in tree_list if token != ""]
        tree_deque = deque(tree_list)

        has_lengths = ":" in tree_deque  # determine whether tree has lengths
        root = self.add_child_subtree_from_newick_deque(None, tree_deque, has_lengths)
        root.write_newick(self.config.aligned_primary_tree_fp, align_tips=True)

        return root

    def add_child_subtree_from_newick_deque(
        self,
        parent: Optional["node_module.Node"],
        tree_deque: deque[str],
        has_lengths: bool,
    ) -> "node_module.Node":
        """Process a deque of Newick tokens to build a tree.

        Each call constructs one subtree and returns its root.
        1. Recursive case
               An open paren indicates a compound subtree.
               The function calls itself to add the first child.
        2. Base case
               An alphanumeric label indicates a leaf.
               Return a simple leaf node.
        3. Following the first child subtree
               There will be an arbitrary number of sibling subtrees,
               each preceeded by a comma.
               The function calls itself to add each in turn.
        4. The end of a subtree
               Signaled by a close paren.
               At this point, add a label and/or length, if either are provided.

        """
        # -------------------------------------------------------------------------
        # First node of subtree
        node = node_module.Node(parent=parent, tree=self)
        token = tree_deque.popleft()
        if token == "(":  # Recursive case: compound subtree
            self.add_child_subtree_from_newick_deque(node, tree_deque, has_lengths)
        else:  # Base case: leaf tree
            node.set_label(token)
            if has_lengths:
                type(self).process_newick_length(node, tree_deque)

            return node

        # -------------------------------------------------------------------------
        # Second through n-th nodes of subtree
        token = tree_deque.popleft()
        while token == ",":
            self.add_child_subtree_from_newick_deque(node, tree_deque, has_lengths)
            token = tree_deque.popleft()

        # -------------------------------------------------------------------------
        # End of subtree
        verify_newick_token(token, ")")
        node.reverse_children()
        token = tree_deque.popleft()
        if token not in self.config.newick_semantic_token_set:
            node.set_label(token)

        if has_lengths and tree_deque[0] != ";":
            self.process_newick_length(node, tree_deque)

        return node

    @classmethod
    def process_newick_length(
        cls,
        node: "node_module.Node",
        tree_deque: deque[str],
    ) -> None:
        """Set branch length from Newick tokens."""

        verify_newick_token(tree_deque.popleft(), ":")  # Next token should be colon
        branch_length = float(tree_deque.popleft())  # Branch length
        node.set_branch_length(branch_length)

    # Import SNPs and assign to branches
    # ----------------------------------------------------------------------
    def import_isogg_snps(self) -> None:
        """Import ISOGG SNPs."""

        snp_module.SNP.set_class_variables(self)
        self.load_preferred_snp_name_set()
        self.load_representative_snp_name_set()
        self.load_isogg_multiallelic_pos_set()
        self.load_isogg_omit_name_set()
        self.load_isogg_omit_pos_str_mutation_set()
        self.load_isogg_corrections()
        self.load_and_parse_isogg_table()
        self.set_depth_first_node_list()
        self.sort_snp_lists_and_set_representatives()
        self.log_isogg_counts()
        self.write_unique_snp_table()
        self.write_newick()
        self.check_multiallelics()

    def load_preferred_snp_name_set(self) -> None:
        """Load preferred SNP names.

        Presence on this list is the primary selection criterion for SNP labels.

        Sets
        ----
        self.preferred_snp_name_set : set[str]
            Set of preferred SNP names.

        """
        for line in load_data_lines(self.config.preferred_snp_names_data_file):
            self.preferred_snp_name_set.add(line.strip())

        logger.info(
            "\nVariant names\n\n"
            "Loaded preferred SNP names\n"
            f"{len(self.preferred_snp_name_set):6d} SNP names: "
            f"{self.config.preferred_snp_names_data_file.filename}\n"
        )

    def load_representative_snp_name_set(self) -> None:
        """Load the names of SNPs deemed representative for their respective lineages.

        Sets
        ----
        self.representative_snp_name_set : set[str]
            Set of names of SNPs deemed representative for their respective lineages.

        """
        counts_dict: dict[str, int] = defaultdict(int)

        set1 = set()
        for line in load_data_lines(self.config.isogg_rep_snp_data_file):
            counts_dict["lines"] += 1
            snp_aliases_string = line.strip().split()[1]
            if snp_aliases_string != ".":
                counts_dict["haplogroups"] += 1
                for snp_aliases in snp_aliases_string.split(","):
                    counts_dict["snps"] += 1
                    for snp_name in snp_aliases.split("/"):
                        set1.add(snp_name)

        set2 = set()
        for line in load_data_lines(self.config.other_rep_snp_data_file):
            set2.add(line.strip().split()[1])

        self.representative_snp_name_set = set1 | set2
        logger.info(
            "Loaded representative SNPs\n"
            f"{counts_dict['lines']:6d} Haplogroups in: "
            f"{self.config.isogg_rep_snp_data_file.filename}\n"
            f"{counts_dict['haplogroups']:6d} "
            "Haplogroups with at least one ISOGG-designated representative SNP\n"
            f"{counts_dict['snps']:6d} "
            "SNPs, as some haplogroups have more than one representative\n"
            f"{len(set1):6d} SNP names, including aliases\n"
            f"{len(set2):6d} Additional representative SNPs loaded from: "
            f"{self.config.other_rep_snp_data_file.filename}\n"
            f"{len(self.representative_snp_name_set):6d} Total SNP names\n"
        )

    def load_isogg_multiallelic_pos_set(self) -> None:
        """Load set of positions to exclude due to multiple alleles.

        Sets
        ----
        self.multiallelic_old_pos_set : set[int]
            Set of positions to exclude due to multiple alleles.

        """
        for line in load_data_lines(self.config.isogg_multiallelic_data_file):
            position = int(line.strip())
            self.multiallelic_old_pos_set.add(position)

    def load_isogg_omit_name_set(self) -> None:
        """Load SNPs to omit from ISOGG database based on name.

        Sets
        ----
        self.isogg_omit_name_set : set[str]
            Set of SNP names to omit.

        """
        for line in load_data_lines(self.config.isogg_omit_name_data_file):
            name = line.strip().split()[0]
            self.isogg_omit_name_set.add(name)

    def load_isogg_omit_pos_str_mutation_set(self) -> None:
        """Load SNPs to omit from ISOGG database based on position and mutation.

        Sets
        ----
        self.isogg_omit_pos_str_mutation_set : set[tuple[str, str]]
            Set of (position_str, mutation) tuples representing SNPs to omit.

        """
        for data_file in self.config.isogg_omit_pos_str_mutation_data_files:
            for line in load_data_lines(data_file):
                line_list = line.strip().split()
                if len(line_list) > 0 and line_list[0] != "#":
                    position_str, mutation = line_list[2:4]
                    self.isogg_omit_pos_str_mutation_set.add((position_str, mutation))

    def load_isogg_corrections(self) -> None:
        """Load SNPs to correct from ISOGG database.

        Sets
        ----
        self.name_correction_df : pd.DataFrame
            DataFrame of SNPs whose names were incorrect in ISOGG table.
            Index: incorrect_name
            Columns: name, haplogroup, position, mutation
        self.isogg_correction_dict : dict[str, tuple[str, str, str]]
            Maps SNP name to a tuple of haplogroup, position, and mutation.
            Either position or mutation (or both) are corrections from the ISOGG table.

        """
        logger.info("\nISOGG Corrections\n")

        data_file = self.config.isogg_correct_name_data_file
        self.name_correction_df = load_dataframe(data_file).set_index("incorrect_name")
        logger.info(f"{len(self.name_correction_df):6d} {data_file.filename}")

        for data_file in self.config.isogg_corrections_data_file_dict.values():
            num_corrections = 0
            for line in load_data_lines(data_file):
                line_list = line.strip().split()
                if len(line_list) > 0 and line_list[0] != "#":
                    num_corrections += 1
                    haplogroup, position_str, mutation, aliases = line_list[1:5]
                    for alias in aliases.split(","):
                        self.isogg_correction_dict[alias] = (
                            haplogroup,
                            position_str,
                            mutation,
                        )

            logger.info(f"{num_corrections:6d} {data_file.filename}")

        total_corrections = len(self.name_correction_df) + len(
            self.isogg_correction_dict
        )
        logger.info(f"\n{total_corrections:6d} total corrections, including aliases\n")

    def load_and_parse_isogg_table(self) -> None:
        """Load and parse ISOGG table."""

        logger.info("\nISOGG variant data\n")
        if self.config.suppress_output:
            isogg_out_file = None
            isogg_drop_out_file = None
        else:
            isogg_out_file = open(self.config.cleaned_isogg_fp, "w")
            isogg_drop_out_file = open(self.config.dropped_isogg_fp, "w")

        dropped_marker_list = []
        for line in load_data_lines(self.config.isogg_data_file, log=True)[1:]:
            line_list = line.split("\t")
            self.isogg_counts_dict["read"] += 1

            # Clean up data row and extract values
            line_list = [element.strip() for element in line_list]
            if line_list[1] == "":  # When present, remove extra tab after SNP name
                del line_list[1]

            if len(line_list) != 6:
                self.isogg_counts_dict["bad_lines"] += 1
                continue

            name, haplogroup, _, _, position_str, mutation = line_list

            # Apply corrections
            try:
                snp_ser = self.name_correction_df.loc[name]
                if (
                    snp_ser["haplogroup"],
                    snp_ser["position"],
                    snp_ser["mutation"],
                ) == (haplogroup, int(position_str), mutation):
                    name = snp_ser["name"]
                    self.num_snps_corrected += 1
            except KeyError:
                pass

            if name in self.isogg_correction_dict:
                haplogroup, position_str, mutation = self.isogg_correction_dict[name]
                self.num_snps_corrected += 1

            # Identify markers to drop
            record_is_bad, marker_is_ok_to_represent_node = self.check_isogg_record(
                name,
                haplogroup,
                position_str,
                mutation,
            )
            if record_is_bad:
                self.isogg_counts_dict["dropped"] += 1
                if isogg_drop_out_file:
                    isogg_drop_out_file.write(
                        f"{name:10s} {haplogroup:25s} {position_str:>8s} {mutation}\n"
                    )

                if marker_is_ok_to_represent_node:
                    dropped_marker = snp_module.DroppedMarker(name, haplogroup, self)
                    dropped_marker_list.append(dropped_marker)

                continue

            # Process retained SNPs
            self.isogg_counts_dict["retained"] += 1
            position = int(position_str)
            if isogg_out_file:
                isogg_out_file.write(
                    f"{name:10s} {haplogroup:25s} {position:8d} {mutation}\n"
                )

            self.construct_snp(name, haplogroup, position, mutation)

        self.add_dropped_markers_to_nodes(dropped_marker_list)
        for file in [isogg_out_file, isogg_drop_out_file]:
            if file is not None:
                file.close()

    def construct_snp(
        self,
        name: str,
        haplogroup: str,
        position: int,
        mutation: str,
    ) -> None:
        """Construct SNP.

        Typically, instantiate a SNP and add it to various containers.
        When SNPs are instantiated, they are added to the tree.
        This process may entail growing the tree to include the corresponding node.
        More specialized things occur if a SNP already exists at this position.

        """
        if self.haplogroup_to_node:
            ancestral, derived = mutation[0], mutation[3]
            snp_key = (haplogroup, position)

            if snp_key in self.snp_dict:  # SNP exists under an alias
                snp = self.snp_dict[snp_key]
                if snp.is_ancestral(ancestral) and snp.is_derived(derived):
                    snp.add_name(name)
                    self.snp_dict[name] = snp
                else:
                    new_snp = snp_module.SNP(
                        name,
                        haplogroup,
                        position,
                        ancestral,
                        derived,
                    )
                    raise ValueError(f"Conflicting SNPs:\n{snp}\n{new_snp}\n")
            else:
                if position in self.snp_dict:  # Another SNP with same position
                    old_snp = self.snp_dict[position]
                    if (
                        ancestral not in old_snp.allele_set
                        or derived not in old_snp.allele_set
                    ):
                        self.multiallelic_new_pos_set.add(position)

                # Typical behavior
                snp = snp_module.SNP(name, haplogroup, position, ancestral, derived)
                self.snp_dict[(haplogroup, position)] = snp
                self.snp_dict[name] = snp
                self.snp_dict[position] = snp
                self.snp_list.append(snp)
                self.snp_pos_set.add(position)
                self.snp_name_set.add(name)
                self.isogg_counts_dict["unique"] += 1

    def add_dropped_markers_to_nodes(
        self,
        dropped_marker_list: list["snp_module.DroppedMarker"],
    ) -> None:
        """Add dropped markers to coresponding nodes."""

        for dropped_marker in dropped_marker_list:
            dropped_marker.add_to_node()

    def sort_snp_lists_and_set_representatives(self) -> None:
        """Sort SNPs by priority ranking and select the best representative.

        Repeat for each Node.

        """
        if not self.depth_first_node_list:
            self.set_depth_first_node_list()

        for node in self.depth_first_node_list:
            node.priority_sort_snp_list_and_set_hg_snp()

    def log_isogg_counts(self) -> None:
        """Log counts of ISOGG SNPs."""

        config = self.config
        counts_dict = self.isogg_counts_dict
        num_alt_names = counts_dict["retained"] - counts_dict["unique"]
        correction_fps_str = ("\n" + " " * 12).join(
            [config.isogg_correct_name_data_file.filename]
            + [
                isogg_corrections_data_file.filename
                for isogg_corrections_data_file in config.isogg_corrections_data_file_dict.values()
            ]
        )
        omit_fps_str = ("\n" + " " * 18).join(
            [config.isogg_omit_name_data_file.filename]
            + [
                isogg_omit_data_file.filename
                for isogg_omit_data_file in config.isogg_omit_pos_str_mutation_data_files
            ]
        )

        log_text_list = [
            f"      {counts_dict['read']:5d} SNPs loaded",
            f"  {self.num_snps_corrected:5d} Corrected based on:\n"
            f"            {correction_fps_str}\n",
            f"- {counts_dict['dropped']:5d} SNPs Dropped",
            f"        {counts_dict['qc']:5d} Flagged as not meeting quality guidelines",
            f"        {counts_dict['approx_loc']:5d} Tree location approximate",
            f"        {counts_dict['provisional']:5d} Removed, flagged as provisional, "
            "or otherwise problematic",
            f"        {counts_dict['non_snp']:5d} Non-SNPs",
            f"        {counts_dict['multiallelic']:5d} Excluded as multiallelic "
            f"based on: {config.isogg_multiallelic_data_file.filename}",
            f"        {counts_dict['duplicated_names']:5d} Duplicated names",
            f"        {counts_dict['omitted']:5d} Explicitly excluded based on:\n"
            f"                  {omit_fps_str}",
            f"- {counts_dict['bad_lines']:5d} Bad lines",
            f"= {counts_dict['retained']:5d} SNPs retained\n",
            f"- {num_alt_names:5d} Alternative names",
            f"= {counts_dict['unique']:5d} Unique SNPs added to the tree\n",
        ]

        if not config.suppress_output:
            log_text_list.extend(
                [
                    "Wrote summary tables:",
                    f"- Dropped:  {config.dropped_isogg_fp}",
                    f"- Retained: {config.cleaned_isogg_fp}",
                    f"- Unique:   {config.unique_isogg_fp}\n",
                ]
            )

        logger.info(("\n" + " " * 4).join(log_text_list))

    def write_unique_snp_table(self) -> None:
        """Sort unique SNP list by phylogeny and position, then write to file."""

        if not self.config.suppress_output:
            self.snp_list = sorted(
                self.snp_list,
                key=attrgetter("dfs_rank", "position"),
            )
            with open(self.config.unique_isogg_fp, "w") as unique_isogg_file:
                for snp in self.snp_list:
                    unique_isogg_file.write(f"{snp.str_with_all_names}\n")

    def check_isogg_record(
        self,
        name: str,
        haplogroup: str,
        position_str: str,
        mutation: str,
    ) -> tuple[bool, bool]:
        """Check an ISOGG record.

        Returns
        -------
        record_is_bad : bool
            When True, do not use this marker for classification.
        name_is_ok_to_represent_node : bool
            When True, if no SNPs are retained for the corresponding node,
            it is OK to use this marker name for the node's hg_snp representation.

        """
        if name.endswith("^"):
            self.isogg_counts_dict["qc"] += 1
            return True, True

        if haplogroup.find("~") >= 0:
            self.isogg_counts_dict["approx_loc"] += 1
            return True, False  # Second value irrelevant: no corresponding node

        if (
            haplogroup.find("Investigation") >= 0
            or haplogroup.find("Notes") >= 0
            or haplogroup.find("Private") >= 0
            or haplogroup.find("Removed") >= 0
            or haplogroup.find("Withdrawn") >= 0
            or haplogroup.find("Freq. Mut.") >= 0
            or len(haplogroup) < 1
        ):
            self.isogg_counts_dict["provisional"] += 1
            return True, False  # Second value irrelevant: no corresponding node

        if (
            len(mutation) != 4
            or mutation.find("?") >= 0
            or position_str.find("..") >= 0
        ):
            self.isogg_counts_dict["non_snp"] += 1
            return True, True

        try:
            position = int(position_str)
        except ValueError:
            logger.info(f"\nERROR. Invalid position: {position}\n")
            return True, False

        if position in self.multiallelic_old_pos_set:
            self.isogg_counts_dict["multiallelic"] += 1
            return True, True

        if name in self.snp_name_set:
            self.isogg_counts_dict["duplicated_names"] += 1
            return True, False

        if (
            name in self.isogg_omit_name_set
            or (position_str, mutation) in self.isogg_omit_pos_str_mutation_set
        ):
            self.isogg_counts_dict["omitted"] += 1
            return True, False

        return False, True

    def check_multiallelics(self) -> None:
        """Check for mutliallelic variants and write list to file."""

        if not self.config.suppress_output and len(self.multiallelic_new_pos_set) > 0:
            with open(self.config.multiallelic_found_fp, "w") as out_file:
                for position in sorted(list(self.multiallelic_new_pos_set)):
                    out_file.write(f"{position:8d}\n")

            num_multiallelic = len(self.multiallelic_new_pos_set)
            logger.info(
                f"\n*** Detected {num_multiallelic} multiallelic positions. ***\n\n"
                "Please do the following and then re-run:\n"
                f"    cat {self.config.multiallelic_found_fp}"
                f" >> {self.config.isogg_multiallelic_data_file.filename}\n\n"
            )

    def find_or_create_node(self, haplogroup: str) -> "node_module.Node":
        """Return Node corresponding to a haplogroup, if it exists.

        If no Node corresponds to the haplogroup, serially split the
        most recent ancestor that does exist until there is a place for a new Node.

        """
        node = None
        if haplogroup in self.haplogroup_to_node:
            node = self.haplogroup_to_node[haplogroup]
        else:
            for num_chars_to_chop in range(1, len(haplogroup)):
                ancestor_string = haplogroup[:-num_chars_to_chop]
                if ancestor_string in self.haplogroup_to_node:
                    ancestor = self.haplogroup_to_node[ancestor_string]
                    node = ancestor.serial_split(haplogroup)
                    break

        if node is None:
            raise ValueError(f"Unplaceable haplogroup: {haplogroup}")

        return node


def verify_newick_token(observed: str, expected: str) -> None:
    """Raise ValueError if observed and expected strings do not match."""

    if observed != expected:
        raise ValueError(
            "Malformed Newick file.\n"
            f"Expected this token: {expected}\n"
            f"Got this one:        {observed}\n"
        )


def get_bounded_subtree(
    root_haplogroup: str,
    target_haplogroup: str,
) -> "node_module.Node":
    """Generate a bounded subtree.

    Set the root node to correspond to `root_haplogroup`.
    Prune children of the target node and children of all nodes not on a direct path
    from the root to the target.

    Parameters
    ----------
    root_haplogroup : str
        Haplogroup corresponding to the root of the subtree.
    target_haplogroup : str
        Haplogroup corresponding to the target node.

    Returns
    -------
    root_node : Node
        Subtree root node.

    """
    tree = Tree()
    root_node = tree.haplogroup_to_node[root_haplogroup]
    target_node = tree.haplogroup_to_node[target_haplogroup]
    for node in root_node.iter_depth_first():
        if node is target_node or node is not node.mrca(target_node):
            node.remove_children()

    return root_node
