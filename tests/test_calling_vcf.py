import pytest

from tests.common import (
    GENOTYPES_1000Y_ALL_BCF_FP,
    GENOTYPES_1000Y_ONE_VCF_FP,
    GENOTYPES_1000Y_SUBSET_BCF_FP,
    HAPLOGROUPS_1000Y_ALL_FP,
    HAPLOGROUPS_1000Y_ONE_FP,
    HAPLOGROUPS_1000Y_SUBSET_FP,
    load_haplogroup_df,
)
from yhaplo.api.call_haplogroups import call_haplogroups
from yhaplo.api.command_line_args import get_command_line_arg_defaults


def test_vcf_input_1000y_single_sample():
    command_line_args = get_command_line_arg_defaults()
    command_line_args.data_fp = GENOTYPES_1000Y_ONE_VCF_FP
    haplogroup_df = call_haplogroups(command_line_args, suppress_output=True)
    expected_haplogroups_df = load_haplogroup_df(HAPLOGROUPS_1000Y_ONE_FP)
    assert haplogroup_df.equals(expected_haplogroups_df)


@pytest.mark.skip("Large fixture. See: ./tests/fixtures/generate_bcf_fixtures.sh")
def test_bcf_input_1000y_subset():
    command_line_args = get_command_line_arg_defaults()
    command_line_args.data_fp = GENOTYPES_1000Y_SUBSET_BCF_FP
    haplogroup_df = call_haplogroups(command_line_args, suppress_output=True)
    expected_haplogroups_df = load_haplogroup_df(HAPLOGROUPS_1000Y_SUBSET_FP)
    assert haplogroup_df.equals(expected_haplogroups_df)


@pytest.mark.skip("Large fixture. See: ./tests/fixtures/generate_bcf_fixtures.sh")
def test_vcf_input_1000y_all():
    command_line_args = get_command_line_arg_defaults()
    command_line_args.data_fp = GENOTYPES_1000Y_ALL_BCF_FP
    haplogroup_df = call_haplogroups(command_line_args, suppress_output=True)
    expected_haplogroups_df = load_haplogroup_df(HAPLOGROUPS_1000Y_ALL_FP)
    assert haplogroup_df.equals(expected_haplogroups_df)
