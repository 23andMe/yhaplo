"""Plot a Newick-formatted tree.

For details, run:
    yhaplo_plot_tree --help

"""

import argparse
import os

from yhaplo.config import Config
from yhaplo.utils.optional_dependencies import optional_import_error_message

try:
    from Bio import Phylo
except ImportError as error:
    error.msg = error.msg + optional_import_error_message(
        "Biopython",
        "plot haplogroup trees",
        "plot",
    )
    raise error


def main() -> None:
    """Plot a Newick-formatted tree."""

    args = get_args()
    phylo_tree = Phylo.read(args.newick_fp, "newick")
    if args.draw:
        Phylo.draw(phylo_tree)
    else:
        Phylo.draw_ascii(phylo_tree)


def get_args() -> argparse.Namespace:
    """Get command-line arguments."""

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-d",
        "--draw",
        action="store_true",
        default=False,
        help="Draw tree, rather than printing ASCII version",
    )
    parser.add_argument(
        "-n",
        "--newick_fp",
        type=str,
        default=os.path.join(
            os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
            "data",
            Config.primary_tree_data_file.data_subdir,
            Config.primary_tree_data_file.filename,
        ),
        help="Path to file containing Newick tree to plot",
    )

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    main()
