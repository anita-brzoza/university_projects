
#===============================================================================
#
# This script contains commands to annotate SNP.
# 
#===============================================================================

## Libraries
library(VariantAnnotation)
library(GenomicFeatures)


## Geting data
vcf <- readVcf("C:/Users/pourp/OneDrive/Pulpit/NC_001144.5.recode.vcf")
annotation_info <- makeTxDbFromGFF("C:/Users/pourp/OneDrive/Pulpit/annotation.gff")

# Annotating
intersect(seqlevels(vcf), seqlevels(annotation_info))
loc_all <- locateVariants(vcf, annotation_info, AllVariants())
