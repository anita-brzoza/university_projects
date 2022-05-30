#!/bin/bash

#================================================================================
#
# The script contains commands to generate the dictionary and index files for 
# reference genome fasta file.
#
#================================================================================


set -e
echo -e "\nCreating the FASTA sequence dictionary file for reference genome..."
gatk CreateSequenceDictionary -R reference_genome.fna

echo -e "\nCreating fasta index file for reference genome..."
samtools faidx reference_genome.fna