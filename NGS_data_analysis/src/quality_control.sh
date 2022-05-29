#!/bin/bash

#================================================================================
#
# The script contains commands to generate:
# 1) separate quality control reports on sample data using FastQC
# 2) a quality control summary report using MultiQC
#
#================================================================================

set -e

echo -e "\nGenerating separate quality control reports on samples ==========================\n"
fastqc ../raw_data/*.fastq -o .

echo -e "\nGenerating a quality control summary report =====================================\n"
multiqc . -o .