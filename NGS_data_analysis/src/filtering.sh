#!/bin/bash

#================================================================================
#
# The script contains commands to:
# 1) delete adapters from sequences.
# 2) generate quality control reports.
#
#================================================================================

set -e 

echo "Deleting adapters from sequences ================================================"

srr_number="SRR06454"
forward="_1"
reverse="_2"
fastq=".fastq"
paired="paired_"
unpaired="unpaired_"

for file in 5 6 7
    do
        java -jar /home/abrzoza/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/trimmomatic.jar \
            PE \
            -phred33 \
            ../raw_data/$srr_number$file$forward$fastq ../raw_data/$srr_number$file$reverse$fastq \
            ./$paired$srr_number$file$forward$fastq ./$unpaired$srr_number$file$forward$fastq \
            ./$paired$srr_number$file$reverse$fastq ./$unpaired$srr_number$file$reverse$fastq \
            ILLUMINACLIP:/home/abrzoza/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:10:3
        
    done

echo 

echo -e "\nGenerating separate quality control reports on samples ==========================\n"
for file in 5 6 7
    do
        fastqc ./$paired$srr_number$file$forward$fastq -o ../fastqc_II
        fastqc ./$paired$srr_number$file$reverse$fastq -o ../fastqc_II
    done

echo -e "\nGenerating a quality control summary report =====================================\n"
multiqc ../fastqc_II -o ../fastqc_II