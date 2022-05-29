#!/bin/bash

#================================================================================
#
# The script contains commands to improve bam files - to sort bam files, to assign
# all the reads in a file to a single new read-group, to mark duplicates, 
# to indexing bam files and calculate read depth.
#
#================================================================================

set -e

for file in 5 6 7
    do
        echo "Sorting bam files for SRR06454${file}"
        samtools sort SRR06454$file.bam -o SRR06454$file.sorted.bam

        echo -e "\nAssigning reads to a single read-group for SRR06454${file}..."
        gatk AddOrReplaceReadGroups \
            -I SRR06454$file.sorted.bam \
            -O SRR06454$file.sorted.rg.bam \
            -ID SRR06454$file \
            -LB SRR06454$file \
            -PL ILLUMINA \
            -PU SRR06454$file \
            -SM SRR06454$file

        echo -e "\nMarking duplicates for SRR06454${file}..."
        gatk MarkDuplicates \
            -I SRR06454$file.sorted.rg.bam \
            -O SRR06454$file.sorted.rg.md.bam \
            -M md_metrics_SRR06454$file.txt
        
        echo -e "\nIndexing bam file for SRR06454${file}..."
        samtools index SRR06454$file.sorted.rg.md.bam
        
        echo -e "\nCalculate read depth for SRR06454${file}"
        samtools depth SRR06454$file.sorted.rg.md.bam -o SRR06454$file.txt
    done
