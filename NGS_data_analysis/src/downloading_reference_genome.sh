#!/bin/bash

#================================================================================
#
# The script contains commands to download reference genome with annotation from 
# NCBI Genome.
#
#================================================================================

echo "Downloading reference genome..."
wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz -O ./reference_genome.fna.gz || { handle ; error ; }
if [ "$?" != 0 ]; then
    echo -e "The file was not downloaded.\nError code: ${?}.\nCheck it out at https://www.gnu.org/software/wget/manual/html_node/Exit-Status.html"
    exit 
else
    echo -e "File has been downloaded successfully."
fi

echo "Downloading annotation..."
wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz -O ./reference_genome.gff.gz || { handle ; error ; }
if [ "$?" != 0 ]; then
    echo -e "The file was not downloaded.\nError code: ${?}.\nCheck it out at https://www.gnu.org/software/wget/manual/html_node/Exit-Status.html"
    exit 
else
    echo -e "File has been downloaded successfully."
fi