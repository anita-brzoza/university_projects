#!/bin/bash

#================================================================================
#
# The script contains commands to download:
# 1. the Accession List containing the identification numbers of the data
# files which contain the raw sequence data for SRP003355 SRA study;
# 2. the first 1,000,000 sequences found in each of the files containing 
# raw sequence data.
#
#================================================================================


echo "Downloading Accession List ===================================================="
accession_list=`esearch -db sra -query SRP003355 \
    | efetch -format runinfo -mode xml \
    | xtract -pattern SraRunInfo -element Run`

echo -e "Accession List:\n $accession_list\n"

echo "Downloading the first 1,000,000 raw sequences ================================="
for srr_number in $accession_list;
    do
        echo "Downloading 1,000,000 raw sequences for ${srr_number}";
        fastq-dump -X 1000000 ${srr_number} --split-files
    done
