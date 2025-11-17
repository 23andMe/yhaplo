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
ttam_data_fp=./data/example.23andMe.txt
multi_sample_bcf_fp=./data/1000Y.all.bcf  # See: tests/fixtures/generate_bcf_fixtures.sh
expected_output_dir=output.expected  # See: README.23andMe.txt

# Output
output_dir=./output
nwk_fp=${output_dir}/y.tree.primary.aligned.ycc.2016.01.04.nwk
tree_drawing_fp=${nwk_fp%.nwk}.drawing.txt

# Colors
CYAN="\033[0;36m"
GREEN="\033[0;32m"
BOLD_CYAN="\033[1;36m"
BOLD_GREEN="\033[1;32m"
RESET_COLOR="\033[0m"

# Helper functions
function echo-run {
    # Echo and run a command.
    # Usage: echo-run "executable --options arguments"

    local command="$1"

    echo -e "\n${CYAN}$ ${RESET_COLOR}${command}"
    eval "${command}"
}
function echo-split-run {
    # Echo split version of command, then run.
    # Usage: echo-split-run "VAR_NAME=value command --options <br> arguments"

    local command="$1"
    echo-split "${command}"
    eval $(echo ${command} | sed 's|<br> ||g')
    echo
}
function echo-split {
    # Split a long command string across lines.
    #
    # 1. Compress spaces.
    # 2. Postpend VAR_NAME=value instances with "\\n" (a backslash and a line break).
    # 3. Prepend instances of {"--", "&&"} with "\\n    ".
    # 4. Replace instances of "<br> " with "\\n    ".

    local command="$1"
    echo
    echo "${command}" \
    | tr -s ' ' \
    | sed 's|[A-Z0-9_]=\S* |&\\\n|g' \
    | sed 's/\(--\|&&\)/\\\n    &/g' \
    | sed 's|<br> |\\\n    |g'
    echo
}



echo -e "\n${BOLD_GREEN}Validating Yhaplo${RESET_COLOR}\n"

echo -e "${BOLD_CYAN}Clearing previous output${RESET_COLOR}..."
echo-run "rm -fr ${output_dir}"


echo -e "\n\n${BOLD_CYAN}Text Input${RESET_COLOR}"
echo-split-run "yhaplo --example_text \
    --hg_genos Q \
    --breadth_first \
    --depth_first \
    --depth_first_table \
    --mrca Q-M3 R-V88 \
    --haplogroup_query E1b1,Q-M3,foo \
    --snp_query L1335,S730,S530,foo"


echo -e "${BOLD_CYAN}Single-Sample VCF Input${RESET_COLOR}"
echo-run "yhaplo --example_vcf --hg_genos Q"
echo -e ""


echo -e "${BOLD_CYAN}Multi-Sample VCF Input${RESET_COLOR}"
if [ ${INCLUDE_BIG_VCF:-""} ]; then
    if [ -e ${multi_sample_bcf_fp} ]; then
        echo-run "yhaplo --input ${multi_sample_bcf_fp} --hg_genos Q"
    else
        echo -e "File not found: ${GREEN}${multi_sample_bcf_fp}${RESET_COLOR}"
        echo "See: tests/fixtures/generate_bcf_fixtures.sh"
    fi
else
    echo -e "\nSkipping. To test multi-sample VCF, use -m option."
fi
echo -e "\n"


echo -e "${BOLD_CYAN}Tree Plotter${RESET_COLOR}"
echo-run "yhaplo-plot-tree --newick_fp ${nwk_fp} | tee ${tree_drawing_fp}"
echo -e "\n"


echo -e "${BOLD_CYAN}Format Converter\n${RESET_COLOR}"
if [ -e ${ttam_data_fp} ]; then
    echo-run "yhaplo-convert-to-genos ${ttam_data_fp}"
    echo-run "mkdir -p ${output_dir}"
    echo-run "mv converted/* ${output_dir}/"
    echo-run "rmdir converted"
else
    echo -e "Skipping. File not found: ${GREEN}${ttam_data_fp}${RESET_COLOR}"
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
