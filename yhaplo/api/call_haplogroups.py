"""Call haplogroups."""

import argparse
import logging
from collections.abc import Iterable, Mapping
from typing import Optional, Union

import pandas as pd

from yhaplo.config import ABLOCK_TYPE, IID_TYPE, Config
from yhaplo.sample import call_haplogroups_from_config


def call_haplogroups(
    command_line_args: Optional[argparse.Namespace] = None,
    iid_to_ablock: Optional[Mapping[IID_TYPE, Union[bytes, ABLOCK_TYPE]]] = None,
    iid_to_platforms: Optional[Mapping[IID_TYPE, Union[str, Iterable[str]]]] = None,
    suppress_output: bool = False,
    out_dir: Optional[str] = None,
    all_aux_output: bool = False,
    root_logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """Configure run, build tree, and call haplogroups.

    Parameters
    ----------
    command_line_args : argparse.Namespace | None, optional
        Command-line arguments.
    iid_to_ablock : Mapping[IID_TYPE, bytes | ABLOCK_TYPE] | None, optional
        Mapping of individual identifiers to 23andMe ablocks.
    iid_to_platforms : Mapping[IID_TYPE, str | Iterable[str]] | None, optional
        Mapping of individual identifier to 23andMe genotyping platforms,
        each starting with "v". Values can be a comma-separated string
        or an iterable of strings.
    suppress_output : bool, optional
        When True, do not generate output files.
    out_dir : str | None, optional
        Output directory.
        When not None, override command_line_args value.
    all_aux_output : bool = False, optional
        Generate all auxiliary output.
        When True, override command_line_args value.
    root_logger : logging.Logger | None, optional
        If supplied, add a file handler.
        To populate the log file, set the logging level to INFO or lower.

    Returns
    -------
    haplogroup_df : pd.DataFrame
        DataFrame of haplogroup calling results.
        Index: Individual identifier.
        Columns:
        - "hg_snp_obs": Haplogroup using a variant of representative-SNP form.
               Rather than using one representative SNP per haplogroup,
               use the most highly ranked SNP this individual was observed
               to carry in the derived state.
        - "hg_snp": Haplogroup in representative-SNP form (e.g., "Q-M3").
        - "ycc_haplogroup": Haplogroup using YCC nomenclature (e.g., "Q1a2a1a1").

    """
    config = Config(
        command_line_args=command_line_args,
        iid_to_ablock=iid_to_ablock,
        iid_to_platforms=iid_to_platforms,
        suppress_output=suppress_output,
        out_dir=out_dir,
        all_aux_output=all_aux_output,
        root_logger=root_logger,
    )
    haplogroup_df = call_haplogroups_from_config(config)

    return haplogroup_df
