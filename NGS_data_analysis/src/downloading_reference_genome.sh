#!/bin/bash

#================================================================================
#
# The script contains commands to download reference genome from Saccharomyces
# Genome Database.
#
#================================================================================

echo "Downloading reference genome..."
wget -q http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_Current_Release.tgz -O ./reference_genome.tgz || { handle ; error ; }
if [ "$?" != 0 ]; then
    echo -e "The file was not downloaded.\nError code: ${?}.\nCheck it out at https://www.gnu.org/software/wget/manual/html_node/Exit-Status.html"
    exit 
else
    echo -e "File has been downloaded successfully."
    echo "Extracting files..."
    tar zxvf reference_genome.tgz
fi

