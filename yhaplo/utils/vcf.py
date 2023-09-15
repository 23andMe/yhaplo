"""VCF utilities."""

import os


def check_vcf_index(vcf_fp: str) -> None:
    """Check that VCF or BCF index file is present.

    Raises
    ------
    FileNotFoundError
        If the VCF/BCF file or its corresponding index file is not present.
    ValueError
        If the file path does not have a .vcf.gz or .bcf extension.

    """
    if not os.path.isfile(vcf_fp):
        raise FileNotFoundError(f"VCF/BCF file not found: {vcf_fp}")

    if vcf_fp.endswith(".vcf.gz"):
        index_fp = vcf_fp.replace(".vcf.gz", ".vcf.gz.tbi")
        if not os.path.isfile(index_fp):
            raise FileNotFoundError(f"VCF index file not found: {index_fp}")

    elif vcf_fp.endswith(".bcf"):
        index_fp = vcf_fp.replace(".bcf", ".bcf.csi")
        if not os.path.isfile(index_fp):
            raise FileNotFoundError(f"BCF index file not found: {index_fp}")

    else:
        raise ValueError(
            f"VCF/BCF file path does not have .vcf.gz or .bcf extension: {vcf_fp}"
        )
