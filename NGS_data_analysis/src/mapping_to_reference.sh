#!/bin/bash

#================================================================================
#
# The script contains commands to:
# 1) index reference genome,
# 2) mapping sequences to reference genome and convert .sam files to .bam files
#
#================================================================================

echo "Indexing reference genome..."
bowtie2-build ./reference_genome.fna ./index/index


forward="_1"
reverse="_2"

for file in 5 6 7 
do
    echo "Mapping sequences to reference genome and converting .sam files to .bam files for SRR06454${file}..."
    bowtie2 -x ./index/index -q -1 ../trimmomatic/paired_SRR06454$file$forward.fastq -2 ../trimmomatic/paired_SRR06454$file$reverse.fastq \
        | samtools view -S -b > SRR06454$file.bam
done