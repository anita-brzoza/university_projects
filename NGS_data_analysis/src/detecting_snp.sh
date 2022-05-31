#!/bin/bash

#================================================================================
#
# The script contains commands to searching polymorphisms, merge files, 
# genotyping and selecting SNPs.
#
#================================================================================

set -e

for file in 5 6 7
    do
        echo -e "\nSearching polymorphisms for SRR06454${file}..."
        gatk HaplotypeCaller \
            -R ../reference_genome/reference_genome.fna \
            -I ../reference_genome/SRR06454$file.sorted.rg.md.bam \
            -O SRR06454$file.vcf.gz \
            -ERC GVCF
    done

echo -e "\nMerging files..."
gatk CombineGVCFs \
    -R ../reference_genome/reference_genome.fna \
    -V SRR064545.vcf.gz \
    -V SRR064546.vcf.gz \
    -V SRR064547.vcf.gz \
    -O merged.vcf.gz

echo -e "\nGenotyping..."
gatk GenotypeGVCFs \
    -R ../reference_genome/reference_genome.fna \
    -V merged.vcf.gz \
    -O merged_genotyped.vcf.gz

echo -e "\nSelecting SNPs..."
gatk SelectVariants \
    -R ../reference_genome/reference_genome.fna \
    -V merged_genotyped.vcf.gz \
    --select-type-to-include SNP \
    -O SNP.vcf.gz