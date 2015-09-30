#source("http://bioconductor.org/biocLite.R")
#biocLite("rtracklayer")
library(ggplot2)
library(rtracklayer)
source("libgridss.R")
source("libneochromosome.R")
source("libvcf.R")


#setwd("W:/dev/gridss/src/test/R/")

theme_set(theme_bw())

millsgr <- bedpe2grmate("~/na12878/lumpy-Mills2012-call-set.bedpe")
pacbiogr <- bedpe2grmate("~/na12878/lumpy-PacBioMoleculo-call-set.bedpe")
vcf <- readVcf("~/na12878/garvan-na12878.vcf", "hg19")

vcf <- vcf[rowRanges(vcf)$QUAL >= 100,]
vcf <- vcf[rowRanges(vcf)$FILTER %in% c(".", "PASS"),]
vcf <- vcf[gridss.vcftodf(vcf)$assembly %in% c("Both"), ] # high qual
vcf <- gridss.removeUnpartnerededBreakend(vcf)
vcf <- vcf[seqnames(rowRanges(vcf)) == seqnames(rowRanges(vcf)[as.character(info(vcf)$MATEID),]),] # intra-chromsomal events

seqlevels(millsgr) <- sub("chr", "", seqlevels(millsgr))
mills <- gridss.annotateBreakpoints(millsgr, millsgr[millsgr$mate,], vcf, maxgap=100, ignore.strand=TRUE)

seqlevels(pacbiogr) <- sub("chr", "", seqlevels(pacbiogr))
pb <- gridss.annotateBreakpoints(pacbiogr, pacbiogr[pacbiogr$mate,], vcf, maxgap=100)
pb$gridss <- pb$gridss[!is.na(pb$gridss$bedid) | pb$gridss$size >= 50,]
