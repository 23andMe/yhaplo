"""Convert data to .genos.txt format for Yhaplo.

Input format options:

1. .ped and .map

2. .23andMe.txt
    Column 1: SNP identifier (ignored)
    Column 2: Chromosome (row retained only if chromosome in CHROMOSOME_SET)
    Column 3: Physical coordinate (GRCh37 assumed)
    Column 4: Allele 1 (row retained only if allele 1 in ALLELE_SET)
    Column 5: Allele 2 (if present)

For details, run:
    yhaplo_convert_to_genos --help

"""

import argparse
import os
import sys
from typing import TextIO

CHROMOSOME_SET = {"24", "Y"}
ALLELE_SET = set("ACGTDI")
FN_ENDING_TO_FN_TYPE = {
    ".ped": "ped",
    ".23andMe.txt": "ttam",
    ".acom.txt": "ttam",
}
OUT_DIR = "converted"


# ----------------------------------------------------------------------
# Ped and Map


def convert_ped(
    ped_fp: str,
    fn_root: str,
    fn_ending: str,
) -> None:
    """Convert a .ped and .map to a .genos.txt."""

    map_fp = ped_fp.replace(fn_ending, ".map")
    out_fp = os.path.join(OUT_DIR, f"{fn_root}.genos.txt")

    with open(out_fp, "w") as out_file:
        index_list = read_map(map_fp, out_file)
        process_ped(ped_fp, index_list, out_file)

    print(f"Output: {out_fp}\n")


def read_map(
    map_fp: str,
    out_file: TextIO,
) -> list[int]:
    """Read a .map file and write positions."""

    if not os.path.exists(map_fp):
        raise FileNotFoundError(f"Map file not found: {map_fp}")

    print(f"Map: {map_fp}\n")

    position_list = []
    index_list = []
    index = 0
    with open(map_fp) as map_file:
        for line in map_file:
            chromosome, _, _, position = line.strip().split()
            if chromosome in CHROMOSOME_SET:
                position_list.append(position)
                index_list.append(index)
            index += 1

    positions_str = "\t".join(position_list)
    out_file.write(f"ID\t{positions_str}\n")

    return index_list


def process_ped(
    ped_fp: str,
    index_list: list[int],
    out_file: TextIO,
) -> None:
    """Process a .ped file."""

    diploid_index_list = [2 * i for i in index_list]
    num_individuals, num_female = 0, 0
    with open(ped_fp) as in_file:
        for line in in_file:
            line_list = line.strip().split()
            sex = line_list[4]
            if sex == "2":
                num_female += 1
                continue

            diploid_geno_list = line_list[6:]
            haploid_geno_list = []
            for i in diploid_index_list:
                allele1, allele2 = diploid_geno_list[i], diploid_geno_list[i + 1]
                if allele1 in ALLELE_SET and allele1 == allele2:
                    haploid_geno_list.append(allele1)
                else:
                    haploid_geno_list.append(".")

            num_individuals += 1
            iid = "-".join(line_list[:2])
            genotypes_str = "\t".join(haploid_geno_list)
            out_file.write(f"{iid}\t{genotypes_str}\n")

    print(f"{num_female:5d} females ignored")
    print(f"{num_individuals:5d} individuals written")
    print(f"{len(index_list):5d} markers\n")


# ----------------------------------------------------------------------
# 23andMe


def convert_ttam(
    in_fp: str,
    ID: str,
) -> None:
    """Read single-sample flat format and converts to .genos.txt."""

    out_fp = os.path.join(OUT_DIR, f"{ID}.genos.txt")
    geno_tuple_list = []
    num_non_y, num_het_or_no_call = 0, 0
    with open(in_fp) as in_file:
        for line in in_file:
            if line[0] == "#" or line[:4] == "rsid":
                continue

            line_list = line.strip().split()
            num_fields = len(line_list)
            chromosome, position, allele1 = line_list[1:4]
            if num_fields == 5:
                allele2 = line_list[4]
            elif num_fields != 4:
                raise ValueError(
                    f"Encountered line with {num_fields} elements:\n" + line
                )

            if chromosome in CHROMOSOME_SET:
                if allele1 in ALLELE_SET and (num_fields == 4 or allele1 == allele2):
                    geno_tuple_list.append((position, allele1))
                else:
                    num_het_or_no_call += 1
            else:
                num_non_y += 1

    os.makedirs(OUT_DIR, exist_ok=True)
    with open(out_fp, "w") as out_file:
        write_line_from_tuple_list(0, geno_tuple_list, out_file, "ID")
        write_line_from_tuple_list(1, geno_tuple_list, out_file, ID)

    print(f"{num_non_y:6d} non-Y genotypes ignored")
    print(f"{num_het_or_no_call:6d} Y-chromosome genotypes ignored (het or no-call)")
    print(f"{len(geno_tuple_list):6d} written\n")
    print(f"Output: {out_fp}\n")


def write_line_from_tuple_list(
    index: int,
    tuple_list: list[tuple[str, str]],
    out_file: TextIO,
    row_header: str = "",
) -> None:
    """Write one line with the i-th element of each tuple."""

    out_file.write(row_header)
    for my_tuple in tuple_list:
        out_file.write(f"\t{my_tuple[index]}")

    out_file.write("\n")


def main() -> None:
    """Run script."""

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("in_fp", type=str, help="input file name")
    args = parser.parse_args()
    in_fp = args.in_fp

    if not os.path.exists(in_fp):
        sys.exit(f"ERROR. Input file does not exist: {in_fp}")

    print(f"Input: {in_fp}\n")

    fn_type = None
    for fn_ending in FN_ENDING_TO_FN_TYPE:
        if in_fp.endswith(fn_ending):
            fn_type = FN_ENDING_TO_FN_TYPE[fn_ending]
            fn_root = os.path.basename(in_fp).replace(fn_ending, "")
            break

    if fn_type == "ped":
        convert_ped(in_fp, fn_root, fn_ending)
    elif fn_type == "ttam":
        convert_ttam(in_fp, ID=fn_root)
    else:
        sys.exit(
            "ERROR. Input file must be a .ped or a .23andMe.txt "
            + "in the corresponding format"
        )


if __name__ == "__main__":
    main()
