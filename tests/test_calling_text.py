from tests.common import (
    GENOTYPES_1000Y_SUBSET_TEXT_FP,
    HAPLOGROUPS_1000Y_SUBSET_FP,
    load_haplogroup_df,
)
from yhaplo.api.call_haplogroups import call_haplogroups
from yhaplo.api.command_line_args import get_command_line_arg_defaults


def test_text_input_1000y_subset():
    command_line_args = get_command_line_arg_defaults()
    command_line_args.data_fp = GENOTYPES_1000Y_SUBSET_TEXT_FP
    haplogroup_df = call_haplogroups(command_line_args, suppress_output=True)
    expected_haplogroups_df = load_haplogroup_df(HAPLOGROUPS_1000Y_SUBSET_FP)
    assert haplogroup_df.equals(expected_haplogroups_df)
