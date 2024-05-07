#!/usr/bin/env python

"""Update SNP metadata corrections.

Yhaplo uses a specific freeze of ISOGG's Y-chromosome haplogroup tree and SNP Index.
This script compares the frozen SNP Index to the current SNP Index in order to identify
mutations and positions that have been corrected since the freeze.
It then patches the Yhaplo input files that correct for these errors.
The following table includes two examples of such corrections:

Haplogroup         SNP    Position_b37 Mutation_Prev  Mutation_Corrected
------------------------------------------------------------------------
R1b1a2a1a2c1f2a1a  Y4010       8606022         A->C   C->A
R1b1a2a1a2c1k      S730       14316964         C->G   G->C

The script was run in late 2023 to generate patched files for Yhaplo 2.1.
To regenerate these findings, first revert any subsequent patches:
$ git log --oneline \
      yhaplo/data/variants/isogg.correct.mutation.txt \
      yhaplo/data/variants/isogg.correct.position.txt
$ git revert <commit shas>

"""

import logging
import os
import subprocess
from typing import Literal, cast

import click
import numpy as np
import pandas as pd
from tabulate import tabulate

from yhaplo.config import Config
from yhaplo.utils.loaders import (
    SNP_TABLE_COL_NAMES,
    SNP_TABLE_DTYPE_DICT,
    load_dataframe,
    load_haplogroup_df,
    load_yhaplo_unique_snps,
)

# Input
SNP_INDEX_URL = "https://docs.google.com/spreadsheets/d/1UY26FvLE3UmEmYFiXgOy0uezJi_wOut-V5TD0a_6-bE"
INPUT_DIR = "input"
DEFAULT_INPUT_DICT = {
    key: os.path.join(INPUT_DIR, filename)
    for key, filename in [
        ("snp_index", "isogg.snp.index.2023.10.26.csv"),
        ("snp_index_removed", "isogg.snp.index.removed.items.2023.11.06.csv"),
        ("unique", "isogg.snps.unique.2016.01.04.txt"),
    ]
}

# Data
#
# Steps to generate:
# 1. Download 1000 Genomes Y-chromosome data:
#    https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/chrY/
#        ALL.chrY_10Mbp_mask.glia_freebayes_maxLikGT_siteQC.20130502.
#        60555_biallelic_snps.vcf.gz
# 2. Convert to BCF: 1000Y.all.bcf.
# 3. Call haplogroups with yhaplo.
DATA_DIR = "data"
KGP_BCF_FP = os.path.join(DATA_DIR, "1000Y.all.bcf")
KGP_HAPLOGROUPS_FP = os.path.join(DATA_DIR, "haplogroups.1000Y.all.txt")

# Output
OUTPUT_DIR = "output"
LOG_FP = os.path.join(OUTPUT_DIR, "log.snp_metadata_corrections.txt")
CORRECTION_FP_TP = os.path.join(OUTPUT_DIR, "isogg.correct.{target}.txt")
EXTRA_OUTPUT_FP_TP = os.path.join(OUTPUT_DIR, "{target}.{label}.txt")

# SNP table constants
SNP_INDEX_RENAME_DICT = {
    "Name": "name",
    "Subgroup Name": "haplogroup",
    "Build 37 Number": "position",
    "Mutation Info": "mutation",
    "Alternate Names": "aliases",
}
assert list(SNP_INDEX_RENAME_DICT.values()) == SNP_TABLE_COL_NAMES
SNP_TABLE_SORT_COLUMNS = ["haplogroup", "position"]
BASES = ["A", "C", "G", "T"]
MUTATIONS_SET = set([f"{anc}->{der}" for anc in BASES for der in BASES if anc != der])

# Manual adjustments
MANUAL_DO_NOT_CORRECT_NAMES = [
    "M278",  # Validated: G2a2b2a2 15022465 T->G
]

# Validation: We expect these updates to be in the final table.
EXPECTED_NAME_MUTATION_TUPLES = [
    ("Y4010", "C->A"),
    ("S730", "G->C"),
]

logger = logging.getLogger()


# ----------------------------------------------------------------------
# App


@click.command(context_settings={"show_default": True})
@click.help_option("-h", "--help")
@click.option(
    "-i",
    "--snp_index_fp",
    type=str,
    default=DEFAULT_INPUT_DICT["snp_index"],
    help="\b\nFilepath of current ISOGG SNP Index,\n"
    f"downloaded as CSV from:\n{SNP_INDEX_URL}\n\b",
)
@click.option(
    "-p",
    "--prev_unique_fp",
    type=str,
    default=DEFAULT_INPUT_DICT["unique"],
    help="\b\nFilepath of previous table of unique SNPs,\n"
    "generated by running `yhaplo` at the command line.\n\b",
)
def main(
    snp_index_fp: str,
    prev_unique_fp: str,
) -> None:
    """Identify ISOGG SNP Index corrections."""

    logging.basicConfig(format="%(message)s", level=logging.INFO)
    os.makedirs(os.path.dirname(LOG_FP), exist_ok=True)
    logger.addHandler(logging.FileHandler(LOG_FP, "w"))
    update_snp_metadata_corrections(snp_index_fp, prev_unique_fp)
    logger.info(f"\n\nLog: {LOG_FP}\n")


def update_snp_metadata_corrections(
    snp_index_fp: str = DEFAULT_INPUT_DICT["snp_index"],
    prev_unique_fp: str = DEFAULT_INPUT_DICT["unique"],
) -> None:
    """Update Yhaplo's ISOGG correction files.

    Parameters
    ----------
    snp_index_fp : str
        Input filepath of current ISOGG SNP Index,
        downloaded as CSV from `SNP_INDEX_URL`.
    prev_unique_fp : str
        Input filepath of previous table of unique SNPs,
        generated by running `yhaplo` at the command line.

    """
    logger.info("\nSNP metadata corrections\n")
    snp_index_df = load_snp_index(snp_index_fp)
    unique_df = load_yhaplo_unique_snps(prev_unique_fp)
    update_metadata_corrections_file(unique_df, snp_index_df, "position")
    update_metadata_corrections_file(unique_df, snp_index_df, "mutation")
    check_isogg_removals(unique_df)


# ----------------------------------------------------------------------
# Data loaders


def load_snp_index(
    snp_index_fp: str = DEFAULT_INPUT_DICT["snp_index"],
    snp_index_rename_dict: dict[str, str] = SNP_INDEX_RENAME_DICT,
) -> pd.DataFrame:
    """Load ISOGG SNP Index.

    Parameters
    ----------
    snp_index_fp : str
        Filepath of ISOGG SNP Index,
        downloaded as CSV from `SNP_INDEX_URL`.
    snp_index_rename_dict : dict[str, str]
        Mapping to rename columns of ISOGG SNP Index.

    Returns
    -------
    snp_index_df : pd.DataFrame
        ISOGG SNP Index.
        Index: name
        Columns: haplogroup, position, mutation, aliases

    """
    col_names = list(snp_index_rename_dict.values())
    snp_index_df = (
        pd.read_csv(snp_index_fp, skiprows=1)
        .rename(columns=snp_index_rename_dict)
        .dropna(subset=col_names[:-1])  # OK if aliases are missing
        .loc[
            lambda df: ~df["name"].str.startswith("[")
            & ~df["position"].str.contains(r":|;|\.")
            & df["mutation"].isin(MUTATIONS_SET)
        ]
        .filter(col_names, axis="columns")
        .astype(SNP_TABLE_DTYPE_DICT)
        .assign(aliases=lambda df: df["aliases"].str.replace("; ", ","))
        .sort_values("name")
        .reset_index(drop=True)
    )

    log_prefix = f"{len(snp_index_df):7d} current records:"
    logger.info(f"{log_prefix:32s}{snp_index_fp}")

    return snp_index_df


def load_snp_index_removed_items(
    snp_index_removed_fp: str = DEFAULT_INPUT_DICT["snp_index_removed"],
) -> pd.DataFrame:
    """Load ISOGG table of removed SNPs.

    Parameters
    ----------
    snp_index_removed_fp : str
        Filepath of ISOGG-removed SNPs,
        downloaded as CSV from `SNP_INDEX_URL`, sheet 2.

    Returns
    -------
    snp_index_removed_df : pd.DataFrame
        ISOGG table of removed SNPs, with non-specific haplogroup.
        Index: name
        Columns: haplogroup, position, mutation

    """
    rename_dict = {
        "SNP": "name",
        "Haplogroup": "haplogroup",
        "Y-position (GRCh37)": "position",
        "Mutation": "mutation",
    }
    col_name_to_dtype = {
        col_name: "Int64" if col_name == "position" else "string"
        for col_name in rename_dict.values()
    }
    snp_index_removed_df = (
        pd.read_csv(snp_index_removed_fp, skiprows=1)
        .rename(columns=rename_dict)
        .filter(list(rename_dict.values()), axis="columns")
        .dropna()
        .loc[
            lambda df: ~df["position"].str.contains(r":|;|\.")
            & df["mutation"].isin(MUTATIONS_SET)
        ]
        .astype(col_name_to_dtype)
        .assign(
            haplogroup=lambda df: df["haplogroup"].str.replace(
                "[Rr]emoved from ",
                "",
                regex=True,
            )
        )
        .sort_values("name")
    )

    log_prefix = f"{len(snp_index_removed_df):7d} ISOGG-removed records:"
    logger.info(f"{log_prefix:32s}{snp_index_removed_fp}")

    return snp_index_removed_df


def load_yhaplo_snp_corrections(
    data_file_key: Literal["mutation", "position"],
    col_names: list[str] = SNP_TABLE_COL_NAMES,
) -> pd.DataFrame:
    """Load Yhaplo data file detailing one class of ISOGG corrections.

    Parameters
    ----------
    data_file_key : Literal["mutation", "position"]
        Key for Yhaplo DataFile of ISOGG corrections, to be augmented.

    col_names : list[str]
        Column names for corrections DataFrame.

    Returns
    -------
    corrections_df : pd.DataFrame
        Dataframe of ISOGG corrections.

    """
    data_file = Config.isogg_corrections_data_file_dict[data_file_key]
    corrections_df = (
        load_dataframe(data_file, header=None)
        .set_axis(col_names, axis="columns")
        .astype(SNP_TABLE_DTYPE_DICT)
    )
    log_prefix = f"{len(corrections_df):7d} {data_file_key}s known:"
    logger.info(f"{log_prefix:32s}{data_file.package} {data_file.filename}")

    return corrections_df


# ----------------------------------------------------------------------
# Data mungers


def update_metadata_corrections_file(
    unique_df: pd.DataFrame,
    snp_index_df: pd.DataFrame,
    target: Literal["mutation", "position"],
    manual_do_not_correct_names: list[str] = MANUAL_DO_NOT_CORRECT_NAMES,
    expected_name_mutation_tuples: list[
        tuple[str, str]
    ] = EXPECTED_NAME_MUTATION_TUPLES,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Update correction file.

    Parameters
    ----------
    unique_df : pd.DataFrame
        DataFrame of Yhaplo-determined unique ISOGG SNPs.
    snp_index_df : pd.DataFrame
        ISOGG SNP Index.
    target : str
        Name of value for which to update corrections.
    manual_do_not_correct_names : list[str]
        List of SNP names to manually exclude.
    expected_name_mutation_tuples : list[tuple[str, str]]
        List of expected (name, mutation) tuples.

    Returns
    -------
    change_df : pd.DataFrame
        Table of changed values.
    patch_df : pd.DataFrame
        Patch table.
    update_df : pd.DataFrame
        Updated corrections.

    """
    logger.info(f"\n--- {target[0].upper()}{target[1:]} ---\n")
    extant_df = load_yhaplo_snp_corrections(target)
    change_df, patch_df = construct_patch(unique_df, snp_index_df, target)
    if target == "mutation":
        extant_df, patch_df = resolve_extant_vs_patch(extant_df, patch_df, change_df)

    update_df = (
        pd.concat([extant_df, patch_df])
        .sort_values(["haplogroup", "name"])
        .reset_index(drop=True)
    )
    logger.info(f"{len(update_df):7d} preliminary merge")

    concordant_df, discordant_df, missing_alleles_df = compare_alleles_to_seq(update_df)
    write_table(missing_alleles_df.reset_index(), target, "no_seq")
    logger.info(f"\n{missing_alleles_df}\n")

    if target == "mutation":
        check_seq_tuple = check_seq_genotypes_of_discordants(discordant_df, unique_df)
        do_not_correct_names, manual_name_mutation_tuples = check_seq_tuple[1:]
        update_df = drop_names(update_df, do_not_correct_names, log=True)
        update_df = drop_names(update_df, manual_do_not_correct_names, log=True)
        logger.info(f"{len(update_df):7d} revised merge\n")
        update_df = set_mutations_manually(update_df, manual_name_mutation_tuples)
        alleles_no_genos_df = check_seq_genotypes_of_concordants(concordant_df)[1]
        write_table(alleles_no_genos_df.reset_index(), target, "no_seq_genos")
        compare_alleles_to_seq(update_df)
        check_expected(update_df, expected_name_mutation_tuples)
    elif target == "position":
        check_seq_genotypes_of_concordants(concordant_df)

    logger.info("\nFinalizing...\n")
    write_snp_table(update_df, target)
    check_dups(update_df)

    return change_df, patch_df, update_df


def construct_patch(
    unique_df: pd.DataFrame,
    snp_index_df: pd.DataFrame,
    target: Literal["mutation", "position"],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Construct patch table.

    Parameters
    ----------
    unique_df : pd.DataFrame
        DataFrame of Yhaplo-determined unique ISOGG SNPs.
    snp_index_df : pd.DataFrame
        ISOGG SNP Index.
    target : str
        Name of value for which to construct patch.

    Returns
    -------
    change_df : pd.DataFrame
        Table of changed values.
    patch_df : pd.DataFrame
        Patch table.

    """
    if target == "mutation":
        _detect_changes = detect_mutation_changes
    elif target == "position":
        _detect_changes = detect_position_changes
    else:
        raise ValueError(f"Invalid target column name: {target}")

    change_df = _detect_changes(unique_df, snp_index_df)
    patch_df = change_df.rename(
        columns={
            "haplogroup_old": "haplogroup",
            f"{target}_new": target,
            "aliases_old": "aliases",
        }
    ).filter(["name", "haplogroup", "position", "mutation", "aliases"], axis="columns")

    return change_df, patch_df


def detect_mutation_changes(
    unique_df: pd.DataFrame,
    snp_index_df: pd.DataFrame,
) -> pd.DataFrame:
    """Detect mutation changes between Yhaplo table and ISOGG SNP Index.

    Parameters
    ----------
    unique_df : pd.DataFrame
        DataFrame of Yhaplo-determined unique ISOGG SNPs.
    snp_index_df : pd.DataFrame
        ISOGG SNP Index.

    Returns
    -------
    mutation_change_df : pd.DataFrame
        DataFrame of discordant mutations.

    """
    mutation_change_df = detect_changes(unique_df, snp_index_df, "mutation").assign(
        anc_old=lambda df: df["mutation_old"].str.split("->").str[0],
        der_old=lambda df: df["mutation_old"].str.split("->").str[1],
        anc_new=lambda df: df["mutation_new"].str.split("->").str[0],
        der_new=lambda df: df["mutation_new"].str.split("->").str[1],
        flip=lambda df: (df["anc_old"] == df["der_new"])
        & (df["anc_new"] == df["der_old"]),
    )
    write_table(mutation_change_df, "mutation", "update")

    return mutation_change_df


def detect_position_changes(
    unique_df: pd.DataFrame,
    snp_index_df: pd.DataFrame,
) -> pd.DataFrame:
    """Detect position changes between Yhaplo table and ISOGG SNP Index.

    Parameters
    ----------
    unique_df : pd.DataFrame
        DataFrame of Yhaplo-determined unique ISOGG SNPs.
    snp_index_df : pd.DataFrame
        ISOGG SNP Index.

    Returns
    -------
    position_change_df : pd.DataFrame
        DataFrame of discordant positions.

    """
    position_change_df = (
        detect_changes(unique_df, snp_index_df, "position")
        .loc[
            lambda df: is_roughly_same_clade(df["haplogroup_old"], df["haplogroup_new"])
        ]
        .reset_index(drop=True)
    )
    write_table(position_change_df, "position", "update")

    return position_change_df


is_roughly_same_clade = np.vectorize(
    lambda haplogroup_old, haplogroup_new: (haplogroup_old[0] == haplogroup_new[0])
    and (haplogroup_old[1].isalpha() == haplogroup_new[1].isalpha())
)


def detect_changes(
    unique_df: pd.DataFrame,
    snp_index_df: pd.DataFrame,
    target: Literal["mutation", "position"],
) -> pd.DataFrame:
    """Detect value changes between Yhaplo table and ISOGG SNP Index.

    Parameters
    ----------
    unique_df : pd.DataFrame
        DataFrame of Yhaplo-determined unique ISOGG SNPs.
    snp_index_df : pd.DataFrame
        ISOGG SNP Index.
    target : str
        Name of value to assess for changes.

    Returns
    -------
    change_df : pd.DataFrame
        Table of changed values.

    """
    if target == "mutation":
        merge_on_col_name = "position"
    elif target == "position":
        merge_on_col_name = "mutation"
    else:
        raise ValueError(f"Invalid target column name: {target}")

    change_df = (
        unique_df.merge(
            snp_index_df,
            how="left",
            on=["name", merge_on_col_name],
            suffixes=("_old", "_new"),
        )
        .loc[lambda df: df[f"{target}_old"] != df[f"{target}_new"]]
        .filter(
            [
                "name",
                merge_on_col_name,
                "haplogroup_old",
                "haplogroup_new",
                f"{target}_old",
                f"{target}_new",
                "aliases_old",
                "aliases_new",
            ],
            axis="columns",
        )
        .sort_values(["haplogroup_old", "name"])
        .reset_index(drop=True)
    )

    return change_df


def resolve_extant_vs_patch(
    extant_df: pd.DataFrame,
    patch_df: pd.DataFrame,
    change_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Resolve conflicts between extant corrections and preliminary patch.

    Parameters
    ----------
    extant_df : pd.DataFrame
        Dataframe of extant ISOGG corrections.
    patch_df : pd.DataFrame
        Preliminary patch table.
    change_df : pd.DataFrame
        Table of ISOGG changes.

    Returns
    -------
    extant_df : pd.DataFrame
        Pruned extant ISOGG corrections.
    patch_df : pd.DataFrame
        Pruned patch table.

    """
    overlap_names = (
        extant_df.set_index("name")
        .filter(patch_df["name"].to_list(), axis="index")
        .index.to_list()
    )
    logger.info(f"{len(overlap_names):7d} overlap names")
    if overlap_names:
        overlap_change_df = (
            change_df.set_index("name")
            .loc[overlap_names]
            .assign(use_update=lambda df: df["der_new"] != df["anc_old"])
        )
        use_update_names = overlap_change_df.query("use_update").index.to_list()
        use_update_names_str = ", ".join(use_update_names)
        extant_df = drop_names(extant_df, use_update_names)
        retain_names = overlap_change_df.query("not use_update").index.to_list()
        retain_names_str = ", ".join(retain_names)
        patch_df = drop_names(patch_df, retain_names)
        logger.info(
            f"\n{overlap_change_df}\n\n"
            f"{len(use_update_names):7d} use update: {use_update_names_str}\n"
            f"{len(retain_names):7d} drop update: {retain_names_str}\n\n"
            f"{len(extant_df):7d} revised extant\n"
            f"{len(patch_df):7d} revised patch"
        )

    return extant_df, patch_df


def drop_names(
    df: pd.DataFrame,
    names_to_drop: list[str],
    log: bool = False,
) -> pd.DataFrame:
    """Drop DataFrame rows with specified names.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to prune.
    names_to_drop : list[str]
        Name values indicating which rows to drop.
    log : bool
        Whether to log actions.

    Returns
    -------
    pd.DataFrame
        Pruned DataFrame.

    """
    df = df.set_index("name").drop(names_to_drop).reset_index()
    if log:
        names_str = ", ".join(names_to_drop)
        logger.info(f"{len(names_to_drop):7d} SNPs excluded: {names_str}")

    return df


def set_mutations_manually(
    df: pd.DataFrame,
    manual_name_mutation_tuples: list[tuple[str, str]],
) -> pd.DataFrame:
    """Set mutation values manually.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with "name" and "mutation" columns.
    manual_name_mutation_tuples : list[tuple[str, str]]
        List of (name, mutation) tuples to use as manual overrides.

    Returns
    -------
    pd.DataFrame
        Updated DataFrame.

    """
    name_mutation_tuples = []
    for name, mutation in manual_name_mutation_tuples:
        name_mask = df["name"] == name
        if name_mask.sum():
            name_mutation_tuples.append((name, mutation))
            df.loc[name_mask, "mutation"] = mutation

    if name_mutation_tuples:
        name_mutation_str = ", ".join(
            [f"{name} {mutation}" for name, mutation in name_mutation_tuples]
        )
        logger.info(
            f"{len(name_mutation_tuples):7d} SNP mutations set manually: "
            f"{name_mutation_str}"
        )

    return df


# ----------------------------------------------------------------------
# Data checkers


def check_expected(
    snp_df: pd.DataFrame,
    expected_name_mutation_tuples: list[
        tuple[str, str]
    ] = EXPECTED_NAME_MUTATION_TUPLES,
) -> None:
    """Check SNP table for presence of expected values.

    Parameters
    ----------
    snp_df : pd.DataFrame
        SNP table.
    expected_name_mutation_tuples : list[tuple[str, str]]
        List of expected (name, mutation) tuples.

    """
    num_expected = len(expected_name_mutation_tuples)
    logger.info(f"\nValidating {num_expected} expected updates...\n")
    snp_df = snp_df.set_index("name")
    for name, expected_mutation in expected_name_mutation_tuples:
        observed_mutation = snp_df.loc[name, "mutation"]
        prefix = (
            " " * 4
            if observed_mutation == expected_mutation
            else f"*** WARNING *** Observed: {observed_mutation}"
        )
        logger.info(f"{prefix}{name:8s} {expected_mutation}")


def check_dups(df: pd.DataFrame) -> None:
    """Check DataFrame for duplicates in "name" column."""

    num_dup_names = df["name"].duplicated().sum()
    if num_dup_names:
        logger.info("*** WARNING ***")

    logger.info(f"{num_dup_names:7d} duplicate names")


def check_isogg_removals(unique_df: pd.DataFrame) -> pd.DataFrame:
    """Check ISOGG-removed SNPs."""

    logger.info("\n--- ISOGG Removals ---\n")
    snp_index_removed_df = load_snp_index_removed_items()
    remove_df = unique_df.merge(snp_index_removed_df, on="name")
    logger.info(
        f"{len(remove_df):7d} of {len(unique_df)} currently used SNPs\n\n"
        "No action taken."
    )

    return remove_df


# ----------------------------------------------------------------------
# Validation vs. sequencing data


def compare_alleles_to_seq(
    snp_df: pd.DataFrame,
    bcf_fp: str = KGP_BCF_FP,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compare alleles in SNP table to those in sequencing data.

    Parameters
    ----------
    snp_df : pd.DataFrame
        SNP table, with "position" and "mutation" columns.
    bcf_fp : str, optional
        File path of BCF to scan for alleles.

    Returns
    -------
    concordant_alleles_df : pd.DataFrame
        Table of allele data and genotypes, restricted to rows concordant
        between SNP table and sequencing data.
    discordant_alleles_df : pd.DataFrame
        Table of allele data and genotypes, restricted to rows with discordances
        between SNP table and sequencing data.
    missing_alleles_df : pd.DataFrame
        Subset of snp_df for which positions are not present in BCF.

    """
    logger.info("\nChecking alleles versus sequencing data...\n")
    positions = snp_df["position"].to_list()
    snp_df = snp_df.assign(
        anc=lambda df: df["mutation"].str.split("->").str[0],
        der=lambda df: df["mutation"].str.split("->").str[1],
    )
    bcf_alleles_df = get_seq_alleles(positions, bcf_fp)
    snp_seq_df = (
        snp_df.assign(alleles=lambda df: make_sorted_tuple(df["anc"], df["der"]))
        .drop(["aliases"], axis="columns")
        .merge(
            bcf_alleles_df,
            on="position",
            suffixes=("_isogg", "_bcf"),
        )
        .set_index("name")
    )
    concordant_mask = snp_seq_df["alleles_isogg"] == snp_seq_df["alleles_bcf"]
    concordant_alleles_df = snp_seq_df.loc[concordant_mask]
    discordant_alleles_df = snp_seq_df.loc[~concordant_mask]
    missing_alleles_df = snp_df.loc[
        lambda df: ~df["position"].isin(bcf_alleles_df["position"])
    ].set_index("name")
    logger.info(
        f"        {len(concordant_alleles_df):3d} SNPs with concordant alleles\n"
        f"        {len(discordant_alleles_df):3d} SNPs with discordant alleles\n"
        f"{len(missing_alleles_df):7d} positions not present in BCF"
    )

    return concordant_alleles_df, discordant_alleles_df, missing_alleles_df


make_sorted_tuple = np.vectorize(lambda x, y: tuple(sorted([x, y])), otypes=[tuple])


def get_seq_alleles(
    positions: list[int],
    bcf_fp: str = KGP_BCF_FP,
) -> pd.DataFrame:
    """Look up alleles for specified positions in BCF file.

    Parameters
    ----------
    positions : list[int]
        Physical positions of interest.
    bcf_fp : str, optional
        File path of BCF to scan for alleles.

    Returns
    -------
    bcf_alleles_df : pd.DataFrame
        Alleles observed in BCF file.
        Columns: "position", "ref", "alt", "alleles".

    """
    regions_str = ",".join([f"Y:{position}" for position in positions])
    bcftools_cmd = " ".join(
        [
            "bcftools view",
            "--drop-genotypes",
            f"--regions {regions_str}",
            bcf_fp,
            "| bcftools query -f '%POS %REF %ALT'",
        ]
    )
    bcf_alleles_df = run_bcftools(
        bcftools_cmd,
        {"position": "Int64", "ref": "string", "alt": "string"},
    ).assign(alleles=lambda df: make_sorted_tuple(df["ref"], df["alt"]))
    logger.info(f"{len(bcf_alleles_df):7d} allele sets extracted: {bcf_fp}")

    return bcf_alleles_df


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


def check_seq_genotypes_of_discordants(
    discordant_alleles_df: pd.DataFrame,
    unique_df: pd.DataFrame,
    haplogroup_fp: str = KGP_HAPLOGROUPS_FP,
    bcf_fp: str = KGP_BCF_FP,
) -> tuple[pd.DataFrame, list[str], list[tuple[str, str]]]:
    """Check sequencing genotypes for discordant-allele case.

    Parameters
    ----------
    discordant_alleles_df : pd.DataFrame
        Table of allele data, from SNP table and BCF.
        The {anc, der} alleles are assumed to be discordant with {ref, alt}.
    unique_df : pd.DataFrame
        DataFrame of Yhaplo-determined unique ISOGG SNPs.
    haplogroup_fp : str, optional
        File path of previously computed haplogroups.
    bcf_fp : str, optional
        File path of BCF in which to look up genotypes.

    Returns
    -------
    alleles_genos_unique_df : pd.DataFrame
        Table of allele data, sequencing genotypes, and previous metadata.
    do_not_correct_names : list[str]
        Names of SNPs to exclude from corrections files.
    manual_name_mutation_tuples : list[tuple[str, str]]
        List of (name, mutation) tuples to use as manual overrides.

    """
    num_snps = len(discordant_alleles_df)
    logger.info(f"\nChecking sequencing genotypes of {num_snps} discordants...\n")

    genotypes_df = look_up_seq_genotypes(discordant_alleles_df, haplogroup_fp, bcf_fp)
    if genotypes_df.empty:
        return pd.DataFrame(), [], []

    alleles_genos_unique_df = (
        pd.merge(
            discordant_alleles_df.reset_index(),
            genotypes_df.drop(["ref", "alt"], axis="columns"),
            on="position",
            how="outer",
        )
        .set_index("name")
        .join(
            unique_df.set_index("name")
            .rename(columns={"mutation": "mutation_prev"})
            .assign(
                anc_prev=lambda df: df["mutation_prev"].str.split("->").str[0],
                der_prev=lambda df: df["mutation_prev"].str.split("->").str[1],
                alleles_prev=lambda df: make_sorted_tuple(
                    df["anc_prev"],
                    df["der_prev"],
                ),
            )
            .filter(["mutation_prev", "alleles_prev"], axis="columns")
        )
        .assign(
            # For discordants, the observed alt allele is, presumably, derived.
            no_ref_obs=lambda df: ~df["genotypes"].str.contains("0"),
            prev_eq_bcf=lambda df: df["alleles_bcf"] == df["alleles_prev"],
            do_not_correct=lambda df: df["prev_eq_bcf"]
            & df["genotypes"].notna()
            & df["no_ref_obs"],
            set_manually=lambda df: ~df["prev_eq_bcf"],
            investigate=lambda df: ~(df["do_not_correct"] | df["set_manually"]),
        )
        .drop(
            [
                "alleles_isogg",
                "alleles_bcf",
                "alleles_prev",
                "no_ref_obs",
                "prev_eq_bcf",
            ],
            axis="columns",
        )
    )
    logger.info(f"\n{alleles_genos_unique_df}\n")
    do_not_correct_names = alleles_genos_unique_df.query(
        "do_not_correct"
    ).index.to_list()
    manual_name_mutation_tuples = [
        (cast(str, name), f'{row["ref"]}->{row["alt"]}')
        for name, row in alleles_genos_unique_df.query("set_manually").iterrows()
    ]

    return alleles_genos_unique_df, do_not_correct_names, manual_name_mutation_tuples


def check_seq_genotypes_of_concordants(
    concordant_alleles_df: pd.DataFrame,
    haplogroup_fp: str = KGP_HAPLOGROUPS_FP,
    bcf_fp: str = KGP_BCF_FP,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Check sequencing genotypes for concordant-allele case.

    Parameters
    ----------
    concordant_alleles_df : pd.DataFrame
        Table of allele data, from SNP table and BCF.
        The {anc, der} alleles are assumed to be concordant with {ref, alt}.
    haplogroup_fp : str, optional
        File path of previously computed haplogroups.
    bcf_fp : str, optional
        File path of BCF in which to look up genotypes.

    Returns
    -------
    alleles_genos_df : pd.DataFrame
        Table of allele data and sequencing genotypes.
    alleles_no_genos_df : pd.DataFrame
        Table of allele data of SNPs for which the sequencing data includes
        no representatives of the corresponding haplogroup.

    """
    num_snps = len(concordant_alleles_df)
    logger.info(f"\nChecking sequencing genotypes of {num_snps} concordants...\n")

    genotypes_df = look_up_seq_genotypes(
        concordant_alleles_df,
        haplogroup_fp,
        bcf_fp,
        log_bcftools_cmd=False,
    )
    alleles_genos_df = (
        pd.merge(
            concordant_alleles_df.reset_index(),
            genotypes_df.drop(["ref", "alt"], axis="columns"),
            on="position",
            how="inner",
        )
        .set_index("name")
        .drop(["alleles_isogg", "alleles_bcf"], axis="columns")
        .assign(
            der_code=lambda df: (df["der"] == df["alt"]).astype("int").astype("string"),
            num_der=lambda df: count_char(df["genotypes"], df["der_code"]),
            num_obs=lambda df: df["genotypes"].str.len(),
            der_pct=lambda df: (df["num_der"] / df["num_obs"]).round(2),
            validated=lambda df: df["der_pct"] > 0.95,
        )
        .drop(["genotypes"], axis="columns")
    )
    num_validated = alleles_genos_df["validated"].sum()
    num_invalidated = len(alleles_genos_df) - num_validated
    logger.info(
        f"        {num_validated:3d} validated\n"
        f"        {num_invalidated:3d} invalidated\n\n"
        f"{alleles_genos_df}\n"
    )
    if num_invalidated:
        logger.info("*** WARNING ***")

    alleles_no_genos_df = concordant_alleles_df.loc[
        lambda df: ~df["position"].isin(alleles_genos_df["position"])
    ].drop(["alleles_isogg", "alleles_bcf"], axis="columns")
    if alleles_no_genos_df.size:
        logger.info(f"{alleles_no_genos_df}")

    return alleles_genos_df, alleles_no_genos_df


count_char = np.vectorize(lambda string_, char: string_.count(char))


def look_up_seq_genotypes(
    snp_df: pd.DataFrame,
    haplogroup_fp: str = KGP_HAPLOGROUPS_FP,
    bcf_fp: str = KGP_BCF_FP,
    log_bcftools_cmd: bool = True,
) -> pd.DataFrame:
    """Look up genotypes in sequencing data.

    For each row of the DataFrame, look up genotypes for the position of interest,
    in individuals whose previously called haplogroups descend from the clade
    indicated in the row.

    Parameters
    ----------
    snp_df : pd.DataFrame
        SNP table with "haplogroup" and "position" columns.
    haplogroup_fp : str, optional
        File path of previously computed haplogroups.
    bcf_fp : str, optional
        File path of BCF to scan for genotypes.
    log_bcftools_cmd : bool, optional
        Log bcftools command.

    Returns
    -------
    genotypes_df : pd.DataFrame
        DataFrame of bcftools genotype lookups.
        Index: "name"
        Columns: "position", "ref", "alt", "genotypes"

    """
    haplogroup_df = load_haplogroup_df(haplogroup_fp)
    logger.info("")
    genotypes_df_list = []
    num_haplogroups_not_present = 0
    for name, row in snp_df.iterrows():
        haplogroup = row["haplogroup"]
        position = row["position"]
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
        if iids:
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
            if log_bcftools_cmd:
                logger.info(bcftools_cmd)

            genotypes_df_list.append(
                run_bcftools(
                    bcftools_cmd,
                    {
                        "position": "Int64",
                        "ref": "string",
                        "alt": "string",
                        "genotypes": "string",
                    },
                )
                .assign(name=name)
                .set_index("name")
            )
        else:
            num_haplogroups_not_present += 1
            logger.info(f"{haplogroup}: No individuals belong to this clade")

    genotypes_df = pd.concat(genotypes_df_list).assign(
        genotypes=lambda df: df["genotypes"].str.replace(".", "", regex=False)
    )

    if num_haplogroups_not_present:
        logger.info("")
    logger.info(
        f"{num_haplogroups_not_present:7d} with no data for these haplogroups\n"
        f"{len(genotypes_df):7d} with genotype data for these haplogroups"
    )

    return genotypes_df


# ----------------------------------------------------------------------
# Data writers


def write_table(
    df: pd.DataFrame,
    target: Literal["mutation", "position"],
    label: str,
) -> None:
    """Write a table.

    Parameters
    ----------
    df : pd.DataFrame
        Table to write.
    target : Literal["mutation", "position"],
        Corrections file target.
    label : str
        Label, for filename and logging.

    """
    output_fp = EXTRA_OUTPUT_FP_TP.format(target=target, label=label)
    with open(output_fp, "w") as output_file:
        table = tabulate(
            df,  # type: ignore
            tablefmt="plain",
            headers="keys",
            showindex=False,
        )
        output_file.write(f"{table}\n")

    log_prefix = f"{len(df):7d} {target}s {label}:"
    logger.info(f"{log_prefix:32s}{output_fp}")


def write_snp_table(
    df: pd.DataFrame,
    target: Literal["mutation", "position"],
) -> None:
    """Write a SNP table to file.

    Parameters
    ----------
    df : pd.DataFrame
        SNP table to write.
        Required columns: "name", "haplogroup", "position", "mutation", "aliases"
        Other columns will be ignored.
    target : Literal["mutation", "position"],
        Corrections file target.

    """
    correction_fp = CORRECTION_FP_TP.format(target=target)
    os.makedirs(os.path.dirname(correction_fp), exist_ok=True)
    with open(correction_fp, "w") as correction_file:
        for _, row in df.iterrows():
            correction_file.write(
                f'{row["name"]:15} {row["haplogroup"]:25} {row["position"]:>8} '
                f'{row["mutation"]}     {row["aliases"]}\n'
            )

    log_prefix = f"{len(df):7d} {target}s merged:"
    logger.info(f"{log_prefix:32s}{correction_fp}")


# ----------------------------------------------------------------------
if __name__ == "__main__":
    main()
