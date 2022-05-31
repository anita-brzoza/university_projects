#!/bin/bash

#================================================================================
#
# The script contains commands to filter VCF file contains SNP.
#
#================================================================================

set -e

echo -e "\nFiltering VCF file..."

vcftools --vcf SNP.vcf \
    --max-missing 0.8 \
    --remove-indels \
    --maf 0.01 \
    --minQ 30 \
    --min-meanDP 30 \
    --remove-filtered-geno-all \
    --recode \
    --recode-INFO-all \
    --out filtered_SNP

echo -e "\nRemove mitochondrial genome..."
vcftools --vcf filtered_SNP.recode.vcf \
    --chr NC_001144.5 \
    --out NC_001144.5 \
    --recode
