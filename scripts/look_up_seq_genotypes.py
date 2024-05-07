#!/usr/bin/env python

"""Look up genotypes in sequencing data.

Example:
-------
$ look_up_seq_genotypes.py E1b1b1b 7317593

   1244 haplogroups loaded: data/haplogroups.1000Y.all.txt
     19 with a haplogroup in clade E1b1b1b

bcftools view --regions Y:7317593 --samples HG01699,HG02317,HG02798,HG01088,HG01104,\
HG01161,HG01680,HG02150,NA19676,NA19759,NA19792,HG01110,HG01250,HG01325,HG01112,\
HG01455,NA19311,NA19331,NA19334 data/1000Y.all.bcf \
| bcftools query --format '%POS %REF %ALT [%GT]'

         ref alt  num_ref  num_alt  num_miss
position
7317593    C   G       10        8         1

"""

import logging
import os
import subprocess

import click
import numpy as np
import pandas as pd

from yhaplo.utils.loaders import load_haplogroup_df

DATA_DIR = "data"
KGP_BCF_FP = os.path.join(DATA_DIR, "1000Y.all.bcf")
KGP_HAPLOGROUPS_FP = os.path.join(DATA_DIR, "haplogroups.1000Y.all.txt")

logger = logging.getLogger()


@click.command(no_args_is_help=True, context_settings={"show_default": True})
@click.help_option("-h", "--help")
@click.argument("haplogroup", type=str)
@click.argument("position", type=int)
@click.option(
    "-hg",
    "--haplogroup_fp",
    type=str,
    default=KGP_HAPLOGROUPS_FP,
    help="File path of previously computed haplogroups.",
)
@click.option(
    "-b",
    "--bcf_fp",
    type=str,
    default=KGP_BCF_FP,
    help="File path of BCF to scan for genotypes.",
)
def main(
    haplogroup: str,
    position: int,
    haplogroup_fp: str,
    bcf_fp: str,
) -> None:
    """Look up genotypes in sequencing data.

    Select individuals whose previously called haplogroups descend from HAPLOGROUP,
    and look up genotypes at POSITION.

    """
    logging.basicConfig(format="%(message)s", level=logging.INFO)
    look_up_seq_genotypes(haplogroup, position, haplogroup_fp, bcf_fp)


def look_up_seq_genotypes(
    haplogroup: str,
    position: int,
    haplogroup_fp: str = KGP_HAPLOGROUPS_FP,
    bcf_fp: str = KGP_BCF_FP,
) -> pd.DataFrame:
    """Look up genotypes in sequencing data.

    Select individuals whose previously called haplogroups descend from the specified
    haplogroup, and look up genotypes at the specified position.

    Parameters
    ----------
    haplogroup : str
        YCC haplogroup.
    position : int
        GRCh37 Y-chromosome position.
    haplogroup_fp : str, optional
        File path of previously computed haplogroups.
    bcf_fp : str, optional
        File path of BCF to scan for alleles.

    Returns
    -------
    genotypes_df : pd.DataFrame
        DataFrame of bcftools genotype lookups.
        Index: "name"
        Columns: "position", "ref", "alt", "genotypes"

    """
    logger.info("")
    haplogroup_df = load_haplogroup_df(haplogroup_fp)
    haplogroup_ser = haplogroup_df["ycc_haplogroup"]
    # Note: This mask is a hack to cover current needs. It is by no means complete.
    # A more thorough treatment would instantiate the tree and include haplogroups
    # descending from the clade of interest.
    haplogroup_mask = (
        haplogroup_ser.str.startswith(haplogroup)
        | (
            (haplogroup == "A1b")
            & (haplogroup_ser != "A1")
            & ~haplogroup_ser.str.startswith(("A0", "A1a"))
        )
        | ((haplogroup == "BT") & ~haplogroup_ser.str.startswith("A"))
        | ((haplogroup == "CT") & ~haplogroup_ser.str.startswith(("A", "B")))
    )
    iids = haplogroup_df.loc[haplogroup_mask].index.to_list()
    genotypes_df = pd.DataFrame()
    if iids:
        logger.info(f"{len(iids):7d} with a haplogroup in clade {haplogroup}\n")
        iids_str = ",".join(iids)
        bcftools_cmd = " ".join(
            [
                "bcftools view",
                f"--regions Y:{position}",
                f"--samples {iids_str}",
                bcf_fp,
                "| bcftools query --format '%POS %REF %ALT [%GT]'",
            ]
        )
        logger.info(f"{bcftools_cmd}\n")
        bcftools_df = run_bcftools(
            bcftools_cmd,
            {
                "position": "Int64",
                "ref": "string",
                "alt": "string",
                "genotypes": "string",
            },
        )
        if not bcftools_df.empty:
            genotypes_df = (
                bcftools_df.set_index("position")
                .assign(
                    num_ref=lambda df: count_char(df["genotypes"], "0"),
                    num_alt=lambda df: count_char(df["genotypes"], "1"),
                    num_miss=lambda df: count_char(df["genotypes"], "."),
                )
                .drop(["genotypes"], axis="columns")
            )
            logger.info(f"{genotypes_df}\n")
    else:
        logger.info(f"{haplogroup}: No individuals belong to this clade\n")

    return genotypes_df


count_char = np.vectorize(lambda string_, char: string_.count(char))


def run_bcftools(
    bcftools_cmd: str,
    col_name_to_dtype: dict[str, str],
) -> pd.DataFrame:
    """Run bcftools and return DataFrame of output.

    Parameters
    ----------
    bcftools_cmd : str
        Command to run.
    col_name_to_dtype : dict[str, str]
        Maps column name to dtype.

    Returns
    -------
    bcftools_df : pd.DataFrame
        DataFrame of results.

    """
    try:
        bcfools_output = subprocess.run(
            bcftools_cmd,
            capture_output=True,
            text=True,
            shell=True,
            check=True,
        ).stdout
    except subprocess.CalledProcessError as error:
        raise RuntimeError(
            "Unable to run bcftools command.\n"
            f"    {bcftools_cmd}\n"
            f"    Check bcftools installation and availability of BCF file."
        ) from error

    if bcfools_output:
        bcftools_df = pd.DataFrame(
            [record.split() for record in bcfools_output.strip().split("\n")],
            columns=list(col_name_to_dtype.keys()),
        ).astype(col_name_to_dtype)
    else:
        logger.warning("WARNING. bcftools query returned no results.\n")
        bcftools_df = pd.DataFrame()

    return bcftools_df


if __name__ == "__main__":
    main()
