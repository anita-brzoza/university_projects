
#===============================================================================
#
# This script contains commands to annotate SNP and comparative analysis 
# of samples.
# 
#===============================================================================

## Libraries
library(VariantAnnotation)
library(GenomicFeatures)

## annotating SNP ==============================================================
## Geting data
vcf <- readVcf("./NC_001144.5.recode.vcf")
annotation_info <- makeTxDbFromGFF("./annotation.gff")

# Annotating
annotated_variants <- locateVariants(vcf, annotation_info, AllVariants())

## Comparative analysis of samples =============================================
# What genotypes are present
genotypes <- geno(vcf)$GT

# Finding in which regions and genes the gene variants occur
unique_location_of_variants <- unique(annotated_variants@elementMetadata$LOCATION)

unique(annotated_variants@elementMetadata$GENEID)

variants_location <- data.frame(
  Chromosome=annotated_variants@seqnames@values,
  Position=annotated_variants@ranges@start,
  Genome_region=annotated_variants@elementMetadata$LOCATION,
  Gene_ID=annotated_variants@elementMetadata$GENEID
)





