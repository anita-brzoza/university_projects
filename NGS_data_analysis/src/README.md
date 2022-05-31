This folder contains scripts used during the analysis:

* [data_downloading.sh](data_downloading.sh) - with commands to download: 1) the Accession List containing the identification numbers of the data files which contain the raw sequence data for SRP003355 SRA study; 2) the first 1,000,000 sequences found in each of the files containing raw sequence data.
* [quality_control.sh](quality_control.sh) - with commands to generate: 1) separate quality control reports on sample data using FastQC; 2) a quality control summary report using MultiQC.
* [filtering.sh](filtering.sh) - with commands to: 1) delete adapters from sequences, 2) generate quality control reports.
* [downloading_reference_genome.sh](downloading_reference_genome.sh) - with commands to download reference genome with annotation from NCBI Genome.
* [mapping_to_reference.sh](mapping_to_reference.sh) - with commands to: 1) index reference genome, 2) mapping sequences to reference genome and convert .sam files to .bam files
* [improving_bam_files.sh](improving_bam_files.sh) - with commands to improve bam files - to sort bam files, to assign all the reads in a file to a single new read-group, to mark duplicates, to indexing bam files and calculate read depth.
* [preparing_reference_genome.sh](preparing_reference_genome.sh) - with commands to generate the dictionary and index files for reference genome fasta file.
* [detecting_snp.sh](detecting_snp.sh) - with commands to searching polymorphisms, merge files, genotyping and selecting SNPs.
* [filtering_SNP.sh](filtering_SNP.sh) - with commands to filter VCF file contains SNP.
* [annotation_and_comparative_analysis.R](annotation_and_comparative_analysis.R) - with commands to annotate SNP and comparative analysis of samples.