"""Command-line arguments."""

import argparse

from yhaplo import __version__

DESCRIPTION = """
Yhaplo identifies the Y-chromosome haplogroup of each male in a sample of one
to millions. Sequence data will yield the most highly resolved classifications,
but the algorithm also works well with chip-based genotype data, provided a
reasonable number of phylogenetically informative sites have been assayed.
"""

ANC_STOP_THRESH_DEFAULT = 2  # BFS stopping condition parameter default
DER_COLLAPSE_THRESH_DEFAULT = 2  # BFS collapsing parameter default


def get_command_line_args(set_defaults: bool = False) -> argparse.Namespace:
    """Process command-line arguments or set defaults.

    Parameters
    ----------
    set_defaults : bool
        When True, ignore command-line arguments and return defaults.

    Returns
    -------
    args : argparse.Namespace
        Command-line options and arguments.

    """
    parser = argparse.ArgumentParser(
        description=DESCRIPTION,
        formatter_class=RawTextWithDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"yhaplo {__version__}",
    )

    group = parser.add_argument_group("Input")
    group.add_argument(
        "-i",
        "--input",
        dest="data_fp",
        metavar="file_name",
        help="Input file. Formats:\n"
        "* Indexed BCF: .bcf, .bcf.csi\n"
        "* Indexed VCF: .vcf.gz, .vcf.gz.tbi\n"
        "* Sample-major text: .genos.txt or .genos.txt.gz\n"
        "  Row 0: Physical coordinates (GRCh37)\n"
        "  Column 0: Individual identifiers\n"
        "  Cell (i, j): Genotype for individual i at position j.\n"
        '      Values include {"A", "C", "G", "T", "."}, \n'
        '      with "." indicating an unobserved value.',
    )

    group = parser.add_argument_group("Output")
    group.add_argument(
        "-o",
        "--out_dir",
        dest="out_dir",
        metavar="dir_name",
        help="Output directory",
    )

    group = parser.add_argument_group("Example data")
    group.add_argument(
        "-ex-txt",
        "--example_text",
        action="store_true",
        help="Run yhaplo on a subset of 1000 Genomes data\n"
        "and produce all auxiliary output",
    )
    group.add_argument(
        "-ex-vcf",
        "--example_vcf",
        action="store_true",
        help="Run yhaplo on a single-sample 1000 Genomes VCF\n"
        "and produce all auxiliary output",
    )

    group = parser.add_argument_group(
        "Auxiliary output",
        "Generate files detailing haplogroup calling for each individual.",
    )
    group.add_argument(
        "-aao",
        "--all_aux_output",
        action="store_true",
        help="Generate all auxiliary output.\n"
        "Equivalent to these seven options:\n"
        "--ancDerCounts --haplogroupPaths --haplogroupPathsDetail\n"
        "--derSNPs --derSNPsDetail --ancSNPs --ancSNPsDetail",
    )
    group.add_argument(
        "-c",
        "--anc_der_counts",
        dest="write_anc_der_counts",
        action="store_true",
        help="Counts of ancestral and derived alleles encountered\n"
        "at each node visited (omits nodes with zero of each)",
    )
    group.add_argument(
        "-hp",
        "--haplogroup_paths",
        dest="write_haplogroup_paths",
        action="store_true",
        help="Sequence of branch labels from root to call,\n"
        "with counts of derived SNPs observed",
    )
    group.add_argument(
        "-hpd",
        "--haplogroup_paths_detail",
        dest="write_haplogroup_paths_detail",
        action="store_true",
        help="Sequence of branch labels from root to call,\n"
        "with counts of derived SNPs observed and lists thereof",
    )
    group.add_argument(
        "-ds",
        "--der_snps",
        dest="write_der_snps",
        action="store_true",
        help="Lists of derived SNPs on path",
    )
    group.add_argument(
        "-dsd",
        "--der_snps_detail",
        dest="write_der_snps_detail",
        action="store_true",
        help="Detailed information about each derived SNP on path",
    )
    group.add_argument(
        "-as",
        "--anc_snps",
        dest="write_anc_snps",
        action="store_true",
        help="Lists of ancestral SNPs encountered in search",
    )
    group.add_argument(
        "-asd",
        "--anc_snps_detail",
        dest="write_anc_snps_detail",
        action="store_true",
        help="Detailed information about each ancestral SNP\n" "encountered in search",
    )

    group = parser.add_argument_group(
        "Real-time auxiliary output",
        "Write haplogroup calling information as each individual is processed.",
    )
    group.add_argument(
        "-rt",
        "--write_real_time",
        dest="write_haplogroups_real_time",
        action="store_true",
        help="Write haplogroups in real time. includes DFS rank,\n"
        "to sort ex post facto: sort -nk5",
    )
    group.add_argument(
        "-hg",
        "--hg_genos",
        dest="haplogroup_to_list_genotypes_for",
        metavar="haplogroup",
        help="Write genotypes observed for SNPs associated with\n"
        "a specified node of the tree, when it is visited",
    )

    group = parser.add_argument_group("Tree traversal")
    group.add_argument(
        "-b",
        "--breadth_first",
        dest="traverse_bf",
        action="store_true",
        help="Write bread-first traversal",
    )
    group.add_argument(
        "-d",
        "--depth_first",
        dest="traverse_df",
        action="store_true",
        help="Write depth-first (pre-order) traversal",
    )
    group.add_argument(
        "-dt",
        "--depth_first_table",
        dest="write_tree_table",
        action="store_true",
        help="Write depth-first (pre-order) traversal table",
    )
    group.add_argument(
        "-m",
        "--mrca",
        nargs=2,
        dest="mrca_haplogroup_list",
        metavar=("haplogroup1", "haplogroup2"),
        help="Output MRCA of two haplogroups",
    )
    group.add_argument(
        "-sq",
        "--snp_query",
        dest="query_snp_names",
        metavar="snp_names",
        help="List phylogenetic path for each SNP in comma-separated list",
    )
    group.add_argument(
        "-pt",
        "--platform_trees",
        dest="write_platform_trees",
        action="store_true",
        help="23andMe: Write trees whose branch lengths are numbers\n"
        "of platform sites",
    )

    group = parser.add_argument_group("Search parameters")
    group.add_argument(
        "-ast",
        "--anc_stop_thresh",
        dest="anc_stop_thresh",
        metavar="anc_stop_thresh",
        type=int,
        default=ANC_STOP_THRESH_DEFAULT,
        help="BFS ancestral allele stopping condition",
    )
    group.add_argument(
        "-dct",
        "--der_collapse_thresh",
        dest="der_collapse_thresh",
        metavar="der_collapse_thresh",
        type=int,
        default=DER_COLLAPSE_THRESH_DEFAULT,
        help="BFS derived allele collapsing parameter",
    )

    group = parser.add_argument_group("Restrictions")
    group.add_argument(
        "-po",
        "--primary_only",
        action="store_true",
        help="Do NOT import ISOGG SNPs",
    )
    group.add_argument(
        "-r",
        "--root",
        dest="alternative_root",
        metavar="haplogroup",
        help="Start searching tree from this branch",
    )
    group.add_argument(
        "-s",
        "--single_sample",
        dest="single_sample_id",
        metavar="ID",
        help="Restrict to a single sample",
    )

    if set_defaults:
        # Set default values for all options and arguments
        args = parser.parse_args([])
    else:
        # Read options and arguments from from sys.argv[1:]
        args = parser.parse_args()

    return args


def get_command_line_arg_defaults() -> argparse.Namespace:
    """Get default values of command-line arguments.

    Returns
    -------
    args : argparse.Namespace
        Default values of all command-line options and arguments.

    """
    args = get_command_line_args(set_defaults=True)

    return args


class RawTextWithDefaultsHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Help message formatter that retains formatting and adds defaults.

    Combines argparse.RawTextHelpFormatter and argparse.ArgumentDefaultsHelpFormatter.

    """

    def _split_lines(self, text, _):
        return text.splitlines()

    def _get_help_string(self, action):
        help_message = action.help
        if "%(default)" not in action.help:
            if action.default is not argparse.SUPPRESS:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help_message += "\n(default: %(default)s)"

        return help_message
