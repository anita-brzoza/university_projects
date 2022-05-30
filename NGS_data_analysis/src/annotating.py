#================================================================================
#
# This script contains command to annotate VCF file using GFF file.
#
#================================================================================

from bioinfokit.analys import marker
marker.vcf_anot(file='./SNP.vcf', gff_file='./annotation.gff')