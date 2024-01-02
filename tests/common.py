import os

FIXTURES_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "fixtures")
FIXTURES_INPUT_DIR = os.path.join(FIXTURES_DIR, "input")
FIXTURES_OUTPUT_DIR = os.path.join(FIXTURES_DIR, "output")

GENOTYPES_1000Y_ALL_BCF_FP = os.path.join(
    FIXTURES_INPUT_DIR,
    "1000Y.all.bcf",
)
GENOTYPES_1000Y_SUBSET_BCF_FP = os.path.join(
    FIXTURES_INPUT_DIR,
    "1000Y.subset.bcf",
)
GENOTYPES_1000Y_SUBSET_TEXT_FP = os.path.join(
    FIXTURES_INPUT_DIR,
    "1000Y.subset.genos.txt",
)
GENOTYPES_1000Y_ONE_VCF_FP = os.path.join(
    FIXTURES_INPUT_DIR,
    "HG01938.vcf.gz",
)
GENOTYPES_1000Y_ONE_ABLOCK_FP = os.path.join(
    FIXTURES_INPUT_DIR,
    "HG01938.ablock",
)

HAPLOGROUPS_1000Y_ALL_FP = os.path.join(
    FIXTURES_OUTPUT_DIR,
    "haplogroups.1000Y.all.txt",
)
HAPLOGROUPS_1000Y_SUBSET_FP = os.path.join(
    FIXTURES_OUTPUT_DIR,
    "haplogroups.1000Y.subset.txt",
)
HAPLOGROUPS_1000Y_ONE_FP = os.path.join(
    FIXTURES_OUTPUT_DIR,
    "haplogroups.HG01938.txt",
)
