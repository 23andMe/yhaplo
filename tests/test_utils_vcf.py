import pytest

from yhaplo.utils.vcf import check_vcf_index


def test_check_vcf_index_vcf(tmpdir):
    vcf_fn = "foo.vcf.gz"
    vcf_path = tmpdir.join(vcf_fn)
    vcf_path.write("")
    with pytest.raises(FileNotFoundError, match="VCF index"):
        check_vcf_index(str(vcf_path))

    tmpdir.join(f"{vcf_fn}.tbi").write("")
    check_vcf_index(str(vcf_path))


def test_check_vcf_index_bcf(tmpdir):
    bcf_fn = "foo.bcf"
    bcf_path = tmpdir.join(bcf_fn)
    bcf_path.write("")
    with pytest.raises(FileNotFoundError, match="BCF index"):
        check_vcf_index(str(bcf_path))

    tmpdir.join(f"{bcf_fn}.csi").write("")
    check_vcf_index(str(bcf_path))


def test_check_vcf_index_vcf_no_vcf(tmpdir):
    vcf_fn = "foo.vcf"
    vcf_path = tmpdir.join(vcf_fn)
    with pytest.raises(FileNotFoundError, match="VCF/BCF file not found"):
        check_vcf_index(str(vcf_path))


def test_check_vcf_index_vcf_bad_extension(tmpdir):
    vcf_fn = "foo.vcf"
    vcf_path = tmpdir.join(vcf_fn)
    vcf_path.write("")
    with pytest.raises(ValueError, match="extension"):
        check_vcf_index(str(vcf_path))
