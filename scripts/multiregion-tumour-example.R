#install.packages("plyr")
#install.packages("stringr")
#install.packages("reshape")
#install.packages("data.table")
#install.packages("testthat")
#biocLite("VariantAnnotation")
source("libgridss.R") # found in src/test/r

vcf <- readVcf("C:/dev/data/FR07935989_dups.gridss.filt.vcf", "hg19")
matevcf <- vcf[as.character(info(vcf)$MATEID),]
df <- gridss.vcftodf(vcf, allColumns=TRUE)
# this example assume that the normal was assigned to category 0
# and two regions exist
library(ggplot2)
ggplot(df) + aes(x=Q0 / df$QUAL) + geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(0,1))
ggplot(df) + aes(x=Q1 / df$QUAL) + geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(0,1))
ggplot(df) + aes(x=Q2 / df$QUAL) + geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(0,1))

# these tresholds depend on:
# * the amount of tumour contamination in germline 
# * total number of samples (more samples = lower threshold)
# * sequencing depth of each sample (lower normal depth w.r.t tumours = lower threshold)
isSomatic <- df$Q0 / df$QUAL == 0 & (df$Q1 + df$Q2) / df$QUAL >= 0.10
info(vcf)$SOMATIC <- isSomatic
info(vcf)$SOMATIC1 <- isSomatic & df$Q1 / df$QUAL >= 0.10
info(vcf)$SOMATIC2 <- isSomatic & df$Q2 / df$QUAL >= 0.10

# note: this just checks that breakends are close to each other, not the orientation
# of the breakends and will filter large-scale events that just happen to have
# nearby breakend
isSmallEvent <- seqnames(vcf) == seqnames(matevcf) & abs(start(vcf) - start(matevcf)) <= 50

writeVcf(vcf[isSomatic & !isSmallEvent], "FR07935989_dups.gridss.filt.somatic.vcf")
