from yhaplo.config import Config
from yhaplo.tree import Tree
from yhaplo.utils.context_managers import logging_disabled


def test_hg_snp_idempotency():
    config = Config(suppress_output=True)
    with logging_disabled():
        tree_1 = Tree(config)
        tree_2 = Tree(config)

    assert extract_hg_snps(tree_1) == extract_hg_snps(tree_2)


def extract_hg_snps(tree):
    return [node.hg_snp for node in tree.depth_first_node_list]
