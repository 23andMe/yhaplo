#!/bin/bash
#
# David Poznik
# 2020.01.16
# test_github_yhaplo.sh
#
# Runs yhaplo with various configuations:
#     {py2, py3} x {txt, single-sample VCF, multi-sample VCF}
#     By default, multi-sample VCF is skipped. To include, supply argument: m
#
# Runs other scripts from yhaplo package
#
# Compares:
# - expected output to py2 output
# - py2 output to py3 output
#
# To generate expected output, run before making changes.
#----------------------------------------------------------------------
set -o nounset -o pipefail

# input
genos_fn=yhaplo/data/1000Y.subset.genos.txt
single_sample_vcf_fn=yhaplo/data/HG01938.vcf.gz
multi_sample_vcf_fn=data/ALL.chrY_10Mbp_mask.glia_freebayes_maxLikGT_siteQC.20130502.60555_biallelic_snps.vcf.gz
multi_sample_vcf_url=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/chrY/ALL.chrY_10Mbp_mask.glia_freebayes_maxLikGT_siteQC.20130502.60555_biallelic_snps.vcf.gz
ttam_data_fn=data/example.23andMe.txt
nwk_fn=output/y.tree.primary.aligned.ycc.2016.01.04.nwk
expected_output_dir=output.expected

# output
tree_drawing_fn=${nwk_fn%.nwk}.drawing.txt
py2_output_dir=output.yhaplo.py2
py3_output_dir=output.yhaplo.py3

# parameters
python_path=/Users/${USER}/anacondaPY_VERSION/bin/python

rm -fr output
for py_version in 2 3; do
    echo -e "Python ${py_version}\n\n"
    python_ex=${python_path/PY_VERSION/${py_version}}
    run_yhaplo="${python_ex} -m yhaplo.callHaplogroups"
    
    echo -e "Text Input\n"
    ${run_yhaplo} --input ${genos_fn} --ancDerCounts --haplogroupPathsDetail --derSNPsDetail --ancSNPsDetail --hgGenos Q

    echo -e "\n\nSingle-Sample VCF Input\n"
    ${run_yhaplo} -i ${single_sample_vcf_fn} --hgGenos Q

    echo -e "\n\nMulti-Sample VCF Input\n"
    if [ $# == 1 ] && [ "$1" == "m" ]; then
        if [ -e ${multi_sample_vcf_fn} ]; then
            ${run_yhaplo} -i ${multi_sample_vcf_fn} --hgGenos Q
        else
            echo "File not found: ${multi_sample_vcf_fn}"
            echo -e "Download from:\n${multi_sample_vcf_url}"
        fi
    else
        echo "Skipping. To test multi-sample VCF supply argument: m"
    fi

    echo -e "\n\nFormat Converter\n"
    if [ -e ${ttam_data_fn} ]; then
        ${python_ex} -m yhaplo.convert2genos ${ttam_data_fn}
        mv converted/* output/
        rmdir converted
    else
        echo "File not found: ${ttam_data_fn}"
    fi

    echo -e "\n\nTree Plotter\n"
    ${python_ex} -m yhaplo.plotTree -n ${nwk_fn} | tee ${tree_drawing_fn}

    echo
    out_dir=output.yhaplo.py${py_version}
    rm -fr ${out_dir}
    mv output ${out_dir}
done;

echo -e "\n\n\n----------------------------------------------------------------------"
echo -e "Expected versus py2\n"
if [ -d ${expected_output_dir} ]; then
    ls ${expected_output_dir} | xargs -I {} echo diff ${expected_output_dir}/{} ${py2_output_dir}/{} | sh -v
else
    echo "Directory not found: ${expected_output_dir}"
fi


echo -e "\n\n\n----------------------------------------------------------------------"
echo -e "Py2 versus Py3\n"
ls ${py2_output_dir} | xargs -I {} echo diff ${py2_output_dir}/{} ${py3_output_dir}/{} | sh -v
