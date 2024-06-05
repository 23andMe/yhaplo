#!/usr/bin/env bash
#
# This script will:
# 1. Run Yhaplo with various input types and output options.
# 2. Run other project scripts.
# 3. Compare output to expected.
#
# Options:
#     -m           Include multi-sample VCF validation
#
#----------------------------------------------------------------------
set -o nounset -o pipefail

# Command-line arguments and options
if [ $# == 1 ] && [ "$1" == "-m" ]; then INCLUDE_BIG_VCF=1; fi

# Input
ttam_data_fp=data/example.23andMe.txt
multi_sample_bcf_fp=data/1000Y.all.bcf  # See: tests/fixtures/generate_bcf_fixtures.sh
expected_output_dir=output.expected  # See: README.23andMe.txt

# Output
output_dir=output
nwk_fp=${output_dir}/y.tree.primary.aligned.ycc.2016.01.04.nwk
tree_drawing_fp=${nwk_fp%.nwk}.drawing.txt

# Colors
BOLD_CYAN="\033[1;36m"
GREEN="\033[0;32m"
RESET_COLOR="\033[0m"


rm -fr ${output_dir}
echo -e "\n${BOLD_CYAN}Removed${RESET_COLOR}: ${GREEN}${output_dir}\n\n${RESET_COLOR}"


echo -e "${BOLD_CYAN}Text Input${RESET_COLOR}\n"
yhaplo --example_text \
    --hg_genos Q \
    --breadth_first \
    --depth_first \
    --depth_first_table \
    --mrca Q-M3 R-V88 \
    --haplogroup_query E1b1,Q-M3,foo \
    --snp_query L1335,S730,S530,foo
echo -e "\n"


echo -e "${BOLD_CYAN}Single-Sample VCF Input\n${RESET_COLOR}"
yhaplo --example_vcf --hg_genos Q
echo -e "\n"


echo -e "${BOLD_CYAN}Multi-Sample VCF Input\n${RESET_COLOR}"
if [ ${INCLUDE_BIG_VCF:-""} ]; then
    if [ -e ${multi_sample_bcf_fp} ]; then
        yhaplo -i ${multi_sample_bcf_fp} --hg_genos Q
    else
        echo "File not found: ${multi_sample_bcf_fp}"
        echo "See: tests/fixtures/generate_bcf_fixtures.sh"
    fi
else
    echo "Skipping. To test multi-sample VCF, use -m option."
fi
echo -e "\n"


echo -e "${BOLD_CYAN}Tree Plotter\n${RESET_COLOR}"
yhaplo_plot_tree -n ${nwk_fp} | tee ${tree_drawing_fp}
echo -e "\n"


echo -e "${BOLD_CYAN}Format Converter\n${RESET_COLOR}"
if [ -e ${ttam_data_fp} ]; then
    yhaplo_convert_to_genos ${ttam_data_fp}
    mkdir -p ${output_dir}
    mv converted/* ${output_dir}/
    rmdir converted
else
    echo "Skipping. File not found: ${ttam_data_fp}"
fi
echo -e "\n"


echo -e "${BOLD_CYAN}Expected versus Observed\n${RESET_COLOR}"

if [ -d ${expected_output_dir} ]; then
    for fn in $(comm -13 <(ls ${output_dir}/) <(ls ${expected_output_dir}/)); do
        echo "* Not found: ${fn}"
    done
    echo

    for fn in $(ls ${expected_output_dir}); do
        if [ -e ${output_dir}/${fn} ]; then
            echo "Checking: ${fn}"
            diff \
            <(cat ${expected_output_dir}/${fn} \
                | grep -v " | Y-chromosome haplogroup caller") \
            <(cat ${output_dir}/${fn} \
                | grep -v " | Y-chromosome haplogroup caller"  \
                | sed 's| yhaplo\.| ttam.yhaplo.|g' \
                | sed 's|/yhaplo/|/ttam.yhaplo/|g')
        fi
    done
    echo

    for fn in $(comm -23 <(ls ${output_dir}/) <(ls ${expected_output_dir}/)); do
        echo "* Unexpected: ${fn}"
    done
    echo
else
    echo "Directory not found: ${expected_output_dir}"
    echo "See: README.23andMe.txt"
fi
