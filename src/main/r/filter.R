#source("http://bioconductor.org/biocLite.R")
#biocLite("VariantAnnotation")
library(VariantAnnotation)
#setwd("C:\\dev\\idsv\\src\\main\\r")
source("libgridss.R")

args<-commandArgs(trailingOnly=TRUE)
file <- args[q]
print(file)
#vcf <- readVcf("C:\\dev\\idsv\\778.vcf", "hg19_random")
#vcf <- vcf[fixed(vcf)$FILTER == "."] # nrow(vcf[fixed(vcf)$FILTER != "."])
#vcf <- readVcf("C:\\dev\\CPCG0100.vcf", "hg19_random")

#file <- "CPCG0100-no-normal-assembly-support.vcf"
# file <- "W:\\dream\\p23\\gridss.bt2-sorted.bam.bam\\bt2-sorted.bam.bam.vcf"
# file <- "C:\\dev\\dream\\PCSI0072-novo.vcf"
vcf <- readVcf(file, "hg19_random")
df <- gridss.truthdetails.processvcf.vcftodf(vcf)
df[is.na(df)] <- 0
df$roiString <- paste0(seqnames(rowData(vcf)), ":", start(ranges(rowData(vcf))) - 2000, "-", start(ranges(rowData(vcf))) + 2000)

som <- df[
  (df$A_SCNormal + df$A_RPNormal) / (df$A_SC + df$A_RP) < 0.02 & # less than 2% of support is from the normal # grep "A_SC=0," CPCG0100.vcf | grep "A_RP=0," 
  df$RCNormal >= 6 & # min coverage in the normal supporting the reference allele
  df$A_SCTumour / (df$A_SCTumour + df$RCTumour) >= 0.20 & # at least 20% support for variant in the tumour data
#  df$A_RM >=2 & # assembly support from both sides
  df$A_SC >= 2 & # found an exact breakpoint on both sides
  df$A_SC + df$A_RP/2 >= 8 & # at least this many read pairs support the variant
  df$A_MQRT >= 40 & # reasonable mapping on both sides
  1==1,]
# ensure both side pass the filters
som <- som[som$EVENT %in% som$EVENT[duplicated(som$EVENT)],]
pp <- som
roi <- paste(som$roiString, collapse=" ")
somvcf <- vcf[info(vcf)$EVENT %in% som$EVENT,]
info(somvcf)$SOMATIC
#writeVcf(somvcf, paste0(file, "-som.vcf"), nchunk = NA)
# hack workaround for  VariantAnnotation bug
writeLines(c(VariantAnnotation:::.makeVcfHeader(somvcf), VariantAnnotation:::.makeVcfMatrix(NULL, somvcf)),  paste0(file, "-som.vcf"))

