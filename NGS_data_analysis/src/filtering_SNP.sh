#!/bin/bash

#================================================================================
#
# The script contains commands to filter VCF file contains SNP.
#
#================================================================================

set -e

echo -e "\nFiltering VCF file..."

cat SNP_anot.vcf | java -jar /home/abrzoza/snpEff/SnpSift.jar filter " ( QUAL >= 30 ) & ( DP > 10 ) " > filtered_SNP.vcf

