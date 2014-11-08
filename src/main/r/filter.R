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
# file <- "C:\\dev\\dream\\PCSI0048-novo-2.vcf"
vcf <- readVcf(file, "hg19_random")
df <- gridss.truthdetails.processvcf.vcftodf(vcf)
df[is.na(df)] <- 0
df$roiString <- paste0(seqnames(rowData(vcf)), ":", start(ranges(rowData(vcf))) - 2000, "-", start(ranges(rowData(vcf))) + 2000)

som <- df
# less than 2% of support is from the normal # grep "A_SC=0," CPCG0100.vcf | grep "A_RP=0," 
som <- som[(som$A_SCNormal + som$A_RPNormal) / (som$A_SC + som$A_RP) < 0.02,]
# min coverage in the normal supporting the reference allele
som <- som[som$RCNormal >= 6,]
# at least 20% support for variant in the tumour data
som <- som[som$A_SCTumour / (som$A_SCTumour + som$RCTumour) >= 0.20,]
# assembly support from both sides
som <- som[som$A_RM >= 2,]
# found an exact breakpoint on both sides
som <- som[som$A_SC >= 2,]
# at least this many read pairs support the variant
som <- som[som$A_SC + som$A_RP/2 >= 8,]
# reasonably good mapping
som <- som[som$A_MQRT > 30,]
# ensure both side pass the filters
som <- som[som$EVENT %in% som$EVENT[duplicated(som$EVENT)],]
somvcf <- vcf[info(vcf)$EVENT %in% som$EVENT,]
paste(som$roiString, collapse=" ")
nrow(somvcf)
#writeVcf(somvcf, paste0(file, "-som.vcf"), nchunk = NA)
# hack workaround for  VariantAnnotation bug
writeLines(c(VariantAnnotation:::.makeVcfHeader(somvcf), VariantAnnotation:::.makeVcfMatrix(NULL, somvcf)),  paste0(file, "-som.vcf"))


