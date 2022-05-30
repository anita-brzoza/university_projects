
#===============================================================================
#
# This script contains commands to filter SNP, after using SnpSift.
# 
#===============================================================================
## Libraries
lib
## Geting data form a VCF file -------------------------------------------------

# geting columns names
SNP_col <- readLines("C:/Users/pourp/OneDrive/Pulpit/filtered_SNP.vcf")
SNP_col<- SNP_col[-(grep("#CHROM",SNP_col)+1):-(length(SNP_col))]
SNP_col_names <- unlist(strsplit(SNP_col[length(SNP_col)],"\t"))

# geting data 
SNP <- read.table("C:/Users/pourp/OneDrive/Pulpit/filtered_SNP.vcf", stringsAsFactors = FALSE)

# merging columns names with data
names(SNP) <- SNP_col_names

## Deleting variants -----------------------------------------------------------
# with "None" in genomic region column
unique(SNP$`genomic region`)
SNP <- SNP[ !(SNP$`genomic region` == "None"), ]

# with "None" in transcript ID
unique(SNP$`transcript ID`)
SNP <- SNP[ !(SNP$`transcript ID` == "None"), ]

## 

