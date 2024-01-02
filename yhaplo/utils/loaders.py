"""Data loaders."""

import logging
from dataclasses import dataclass
from importlib.resources import as_file, files

import pandas as pd

DATA_SUBPACKAGE = __package__.replace("utils", "data")

SNP_TABLE_COL_NAMES = ["name", "haplogroup", "position", "mutation", "aliases"]
SNP_TABLE_DTYPE_DICT = {
    col_name: "Int64" if col_name == "position" else "string"
    for col_name in SNP_TABLE_COL_NAMES
}

logger = logging.getLogger(__name__)


# ----------------------------------------------------------------------
# General utilities for loading Yhaplo package data.


@dataclass
class DataFile:

    """Attributes of a yhaplo data file.

    data_subdir : str
        Subpackage of DATA_SUBPACKAGE.
    filename : str
        Data filename.
    description : str, optional
        Description of data.
    ttam_only : bool, optional
        Indicates whether data file is a 23andMe-internal data file.

    """

    data_subdir: str
    filename: str
    description: str = "Data"
    ttam_only: bool = False

    @property
    def package(self):
        """Return package name."""

        return f"{DATA_SUBPACKAGE}.{self.data_subdir}"


class TtamFileNotFoundError(FileNotFoundError):

    """Exception indicating that an unfound data file is not publicly available."""

    def __init__(
        self,
        data_file: DataFile,
        package: str,
    ):
        super(TtamFileNotFoundError, self).__init__(
            f'Failed to load "{data_file.filename}" from {package}.\n'
            f"{data_file.description} file only available internally at 23andMe.\n"
        )


def load_data_lines(
    data_file: DataFile,
    log: bool = False,
) -> list[str]:
    """Load yhaplo data file and split into lines."""

    lines = load_data(data_file, log=log).strip().split("\n")
    return lines


def load_data(
    data_file: DataFile,
    log: bool = False,
) -> str:
    """Load yhaplo data file as text.

    Raises
    ------
    TtamFileNotFoundError
        When failing to load a non-public data file.

    """
    try:
        data = files(data_file.package).joinpath(data_file.filename).read_text()
    except (FileNotFoundError, ModuleNotFoundError):
        if data_file.ttam_only:
            raise TtamFileNotFoundError(data_file, data_file.package)
        else:
            raise

    if log:
        logger.info(
            f"Loaded {data_file.description}:\n"
            f"    {data_file.package}: {data_file.filename}\n"
        )

    return data


def load_dataframe(
    data_file: DataFile,
    header="infer",
    log: bool = False,
) -> pd.DataFrame:
    """Load yhaplo data file as pandas DataFrame.

    Raises
    ------
    TtamFileNotFoundError
        When failing to load a non-public data file.

    """
    try:
        with as_file(files(data_file.package).joinpath(data_file.filename)) as path:
            df = pd.read_csv(path, delim_whitespace=True, header=header)

    except (FileNotFoundError, ModuleNotFoundError):
        if data_file.ttam_only:
            raise TtamFileNotFoundError(data_file, data_file.package)
        else:
            raise

    if log:
        logger.info(
            f"Loaded {data_file.description}:\n"
            f"    {data_file.package}: {data_file.filename}\n"
        )

    return df


# ----------------------------------------------------------------------
# Utilities for loading Yhaplo output files.


def load_haplogroup_df(haplogroup_fp: str) -> pd.DataFrame:
    """Load haplogroup table output by Yhaplo.

    Parameters
    ----------
    haplogroup_fp : str
        Path to haplogroup table output by Yhaplo.

    Returns
    -------
    haplogroup_df : pd.DataFrame
        DataFrame of haplogroups calls.
        Index: "iid"
        Columns: "hg_snp_obs", "hg_snp", "ycc_haplogroup"

    """
    haplogroup_df = pd.read_csv(
        haplogroup_fp,
        delim_whitespace=True,
        names=["iid", "hg_snp_obs", "hg_snp", "ycc_haplogroup"],
    ).set_index("iid")

    logger.info(f"{len(haplogroup_df):7d} haplogroups loaded: {haplogroup_fp}")

    return haplogroup_df


def load_yhaplo_unique_snps(unique_fp: str) -> pd.DataFrame:
    """Load Yhaplo-determined unique ISOGG SNPs.

    Parameters
    ----------
    unique_fp : str
        Filepath of unique ISOGG SNPs,
        generated by running `yhaplo` at the command line.

    Returns
    -------
    unique_df : pd.DataFrame
        DataFrame of Yhaplo-determined unique ISOGG SNPs.

    """
    unique_df = pd.read_csv(
        unique_fp,
        delim_whitespace=True,
        header=None,
        names=SNP_TABLE_COL_NAMES,
    ).astype(SNP_TABLE_DTYPE_DICT)

    log_prefix = f"{len(unique_df):7d} SNPs used by yhaplo:"
    logger.info(f"{log_prefix:32s}{unique_fp}")

    return unique_df
