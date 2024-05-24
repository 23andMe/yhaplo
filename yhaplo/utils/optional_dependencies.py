"""Functions for checking whether optional dependencies are available."""

ROOT_PACKAGE = __package__.removesuffix(".utils")  # noqa


def check_vcf_dependencies():
    """Check that "vcf" dependencies are available."""

    try:
        from pysam import VariantFile  # noqa F401
    except ImportError as error:
        error.msg = error.msg + optional_import_error_message(
            "Pysam",
            "process VCF/BCF input",
            "vcf",
        )
        raise error


def optional_import_error_message(
    package_name: str,
    package_use: str,
    optional_dep_category: str,
    root_package: str = ROOT_PACKAGE,
) -> str:
    """Construct message indicating that an optional dependency is required."""

    message = (
        f"\n\n{package_name} is required to {package_use}.\n"
        f'Please re-install Yhaplo with the "{optional_dep_category}" '
        "optional dependencies.\n"
        "For example:\n"
        f"    pip install {root_package}[{optional_dep_category}]\n"
    )
    return message
