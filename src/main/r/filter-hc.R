library(VariantAnnotation)
minQual <- 1000
vcf <- readVcf("W:/778/idsv/778.filtered.vcf", "hg19_random")
hcvcf <- vcf[fixed(vcf)$QUAL > minQual & !is.na(info(vcf)$AS) & !is.na(info(vcf)$RAS),]
writeVcf(hcvcf, "W:/778/idsv/778.hc.vcf")
