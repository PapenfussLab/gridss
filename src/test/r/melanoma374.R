#setwd("W:/dev/gridss/src/test/r")
source("common.R")
source("libvcf.R")
source("libgridss.R")
library(VariantAnnotation)
library(stringr)

vcf <- readVcf(paste0(rootdir, "374/374.gridss.vcf"), "hg19")
vcf <- gridss.removeUnpartnerededBreakend(vcf)
df <- gridss.vcftodf(vcf, allColumns=TRUE)
bpgr <- vcftobpgr(vcf)

vdf <- read.table(paste0(rootdir, "374/socrates_validated.tsv"))
vleftgr <- GRanges(
  seqnames=str_extract(vdf$V1, "^[^:]+"),
  ranges=IRanges(start=as.numeric(str_extract(vdf$V1, "[^:]+$")), width=1),
  strand=vdf$V2)
vrightgr <- GRanges(
  seqnames=str_extract(vdf$V4, "^[^:]+"),
  ranges=IRanges(start=as.numeric(str_extract(vdf$V4, "[^:]+$")), width=1),
  strand=vdf$V5)
vleftgr$mateIndex <- length(vleftgr) + seq_along(vleftgr)
vrightgr$mateIndex <- seq_along(vleftgr)
vgr <- c(vleftgr, vrightgr)
vgr$event <- c(seq(1, length(vleftgr)), seq(1, length(vleftgr)))

hits <- breakpointHits(bpgr, vgr, maxgap=4)
hits <- hits[as.logical(strand(bpgr[hits$queryHits]) == strand(vgr[hits$subjectHits])),]
hitgr <- bpgr[hits$queryHits]
hitgr$event <- vgr[hits$subjectHits]$event

validatedgr <- hitgr[hitgr$event %in% c(1,2,5,11,12,13,18),]
validateddf <- df[names(validatedgr),]


unique(breakpointHits(bpgr, vgr, maxgap=100)$subjectHits)
hitgr <- bpgr[,]
# Not called by gridss
vgr[c(5,17,23,35)]



sdf <- df[df$Q0==0,]
hqdf <- sdf[sdf$FILTER==".",]

writeVcf(vcf[df$Q0==0,], "W:/374/somatic.vcf")
ggplot(sdf) + aes(x=QUAL) + geom_histogram(binwidth=25) + scale_x_continuous(lim=c(200, 1000))
table(sdf$FILTER==".")
#ggplot(sdf) + aes(x=(BUMQ0+BSCQ0) / QUAL, y=QUAL) + geom_point() + scale_x_log10() + scale_y_log10()
#ggplot(sdf[sdf$RP+sdf$SR+sdf$RSR>0,]) + aes(x=(BUM0+BSC0) / (RP+SR+RSR)) + geom_density() + scale_x_log10() + scale_y_log10()
#ggplot(sdf[sdf$RP+sdf$SR+sdf$RSR>0,]) + aes(x=(BUM+BSC) / (RP+SR+RSR)) + geom_histogram(binwidth=1) + scale_x_continuous(lim=c(0,20))
#max((hqdf$BUMQ0 + hqdf$BSCQ0) / hqdf$QUAL)
