import numpy as np

from yhaplo.config import Config
from yhaplo.tree import Tree
from yhaplo.utils.context_managers import logging_disabled


def test_hg_snp_idempotency():
    config = Config(suppress_output=True)
    with logging_disabled():
        tree_1 = Tree(config)
        tree_2 = Tree(config)

    hg_snps_1 = np.array([node.hg_snp for node in tree_1.depth_first_node_list])
    hg_snps_2 = np.array([node.hg_snp for node in tree_2.depth_first_node_list])

    assert (hg_snps_1 == hg_snps_2).all()
