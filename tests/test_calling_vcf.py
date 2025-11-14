import os
import tempfile

import pytest
from pysam import VariantFile, tabix_index

from tests.common import (
    GENOTYPES_1000Y_ALL_BCF_FP,
    GENOTYPES_1000Y_ONE_VCF_FP,
    GENOTYPES_1000Y_SUBSET_BCF_FP,
    HAPLOGROUPS_1000Y_ALL_FP,
    HAPLOGROUPS_1000Y_ONE_FP,
    HAPLOGROUPS_1000Y_SUBSET_FP,
    HAPLOGROUPS_NO_DATA_FP,
)
from yhaplo.api.call_haplogroups import call_haplogroups
from yhaplo.api.command_line_args import get_command_line_arg_defaults
from yhaplo.utils.loaders import load_haplogroup_df


def _test_vcf_yields_expected_haplogroups(vcf_fp, expected_haplogroups_fp):
    command_line_args = get_command_line_arg_defaults()
    command_line_args.data_fp = vcf_fp
    haplogroup_df = call_haplogroups(command_line_args, suppress_output=True)
    expected_haplogroups_df = load_haplogroup_df(expected_haplogroups_fp)
    assert haplogroup_df.equals(expected_haplogroups_df)


def _create_modified_vcf_and_test(genotype_transformer, expected_haplogroups_fp):
    with tempfile.NamedTemporaryFile(
        mode="wb",
        suffix=".vcf.gz",
        delete=False,
    ) as vcf_temp_file:
        temp_vcf_fp = vcf_temp_file.name
        try:
            with VariantFile(GENOTYPES_1000Y_ONE_VCF_FP) as vcf_in:
                header = vcf_in.header.copy()
                with VariantFile(temp_vcf_fp, "w", header=header) as vcf_out:
                    for record in vcf_in:
                        new_record = vcf_out.new_record(
                            contig=record.contig,
                            start=record.start,
                            stop=record.stop,
                            alleles=record.alleles,
                        )
                        for iid in record.samples:
                            genotype = genotype_transformer(record, iid)
                            new_record.samples[iid]["GT"] = genotype

                        vcf_out.write(new_record)

            tabix_index(temp_vcf_fp, preset="vcf", force=True)
            _test_vcf_yields_expected_haplogroups(temp_vcf_fp, expected_haplogroups_fp)
        finally:
            if os.path.exists(temp_vcf_fp):
                os.unlink(temp_vcf_fp)
            if os.path.exists(temp_vcf_fp + ".tbi"):
                os.unlink(temp_vcf_fp + ".tbi")


def test_vcf_input_1000y_single_sample():
    _test_vcf_yields_expected_haplogroups(
        GENOTYPES_1000Y_ONE_VCF_FP,
        HAPLOGROUPS_1000Y_ONE_FP,
    )


def test_vcf_input_1000y_single_sample_diploid_hom():
    def make_homozygous(record, iid):
        gt = record.samples[iid]["GT"]
        genotype = (gt[0], gt[0]) if gt and gt[0] is not None else (None, None)
        return genotype

    _create_modified_vcf_and_test(make_homozygous, HAPLOGROUPS_1000Y_ONE_FP)


def test_vcf_input_single_sample_diploid_het():
    """Confirm that heterozygous genotypes are treated as missing data.

    With no data, the default haplogroups call is "A", the root haplogroup.

    """

    def make_heterozygous(record, iid):
        genotype = record.alleles
        return genotype

    _create_modified_vcf_and_test(make_heterozygous, HAPLOGROUPS_NO_DATA_FP)


@pytest.mark.skip("Large fixture. See: ./tests/fixtures/generate_bcf_fixtures.sh")
def test_bcf_input_1000y_subset():
    _test_vcf_yields_expected_haplogroups(
        GENOTYPES_1000Y_SUBSET_BCF_FP,
        HAPLOGROUPS_1000Y_SUBSET_FP,
    )


@pytest.mark.skip("Large fixture. See: ./tests/fixtures/generate_bcf_fixtures.sh")
def test_vcf_input_1000y_all():
    _test_vcf_yields_expected_haplogroups(
        GENOTYPES_1000Y_ALL_BCF_FP,
        HAPLOGROUPS_1000Y_ALL_FP,
    )
