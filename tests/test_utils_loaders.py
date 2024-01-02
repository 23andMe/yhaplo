import pytest

from yhaplo.config import Config
from yhaplo.utils.loaders import (
    DataFile,
    TtamFileNotFoundError,
    load_data,
    load_dataframe,
)


def test_load_data():
    tree_str = load_data(Config.primary_tree_data_file).strip()
    assert tree_str.startswith("(") and tree_str.endswith(";")


def test_load_dataframe():
    df = load_dataframe(Config.isogg_correct_name_data_file)
    assert df.columns[0] == "incorrect_name"


def test_load_data_missing():
    with pytest.raises(ModuleNotFoundError):
        load_data(DataFile("foo", "bar"))

    with pytest.raises(FileNotFoundError):
        load_data(
            DataFile(
                Config.primary_tree_data_file.data_subdir,
                "foo",
            )
        )


def test_load_data_missing_ttam():
    with pytest.raises(TtamFileNotFoundError):
        load_data(DataFile("foo", "bar", ttam_only=True))

    with pytest.raises(TtamFileNotFoundError):
        load_data(
            DataFile(
                Config.primary_tree_data_file.data_subdir,
                "foo",
                ttam_only=True,
            )
        )
