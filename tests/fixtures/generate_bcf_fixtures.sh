#!/usr/bin/env bash

# Generate two BCF fixtures ("subset", "all") based on the 1000Y VCF.
# The input file is available for download here:
# https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/chrY/

source utils.sh

data_dir=input

# Input
all_vcf_fp=${data_dir}/ALL.chrY_10Mbp_mask.glia_freebayes_maxLikGT_siteQC.20130502.60555_biallelic_snps.vcf.gz
subset_sample_fp=${data_dir}/1000Y.subset.genos.txt

# Output
subset_bcf_fp=${data_dir}/1000Y.subset.bcf
all_bcf_fp=${data_dir}/1000Y.all.bcf

echo -e "\nGenerating ${subset_bcf_fp}..."
samples=$(awk '(NR == 2){ printf $1 }(NR > 2){ printf ","$1 }' ${subset_sample_fp})
echo_run "bcftools convert --samples ${samples} --output-type b --output ${subset_bcf_fp} ${all_vcf_fp}"
echo_run "bcftools index ${subset_bcf_fp}"

echo -e "\nGenerating ${all_bcf_fp}..."
echo_run "bcftools convert --output-type b --output ${all_bcf_fp} ${all_vcf_fp}"
echo_run "bcftools index ${all_bcf_fp}"

echo -e "\nDone!\n"
