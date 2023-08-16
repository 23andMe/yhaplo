"""Data loaders."""

import importlib.resources
import logging
from typing import NamedTuple

DATA_SUBPACKAGE = __package__.replace("utils", "data")

logger = logging.getLogger(__name__)


class DataFile(NamedTuple):

    """Attributes of a yhaplo data file."""

    data_subdir: str
    filename: str
    description: str = "Data"
    ttam_only: bool = False


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
    """Load yhaplo data file.

    Raises
    ------
    TtamFileNotFoundError
        When failing to load a non-public data file.

    """
    package = f"{DATA_SUBPACKAGE}.{data_file.data_subdir}"

    try:
        data = (
            importlib.resources.files(package).joinpath(data_file.filename).read_text()
        )
    except (FileNotFoundError, ModuleNotFoundError):
        if data_file.ttam_only:
            raise TtamFileNotFoundError(data_file, package)
        else:
            raise

    if log:
        logger.info(
            f"Loaded {data_file.description}:\n    {package}: {data_file.filename}\n"
        )

    return data
