#source("http://bioconductor.org/biocLite.R")
#biocLite("rtracklayer")
library(rtracklayer)
source("../../main/r/libgridss.R")
source("libneochromosome.R")
source("libvcf.R")

#setwd("C:/dev/idsv/src/test/R/")

theme_set(theme_bw())

sample <- "778"
rp <- getrpcalls("C:/dev/neochromosome/mmc3.xlsx", "778 (DR)")
rp <- rp[!(seqnames(rp) == seqnames(rp[rp$mate,]) & seqnames(rp) %in% c("chr7", "chr22")),] # remove 7 and 22 intrachromosomal arrangements that had really low limits
cgr <- getcgr("C:/dev/neochromosome/mmc4.xlsx", "778_CGRs")
cn <- getcn("C:/dev/neochromosome/mmc4.xlsx", "778_CN")
vcf <- readVcf("W:/Papenfuss_lab/projects/liposarcoma/data/gridss/778/778.vcf", "hg19")
out <- go(sample, vcf, rp, cgr, minimumEventSize=500)

sample <- "T1000"
rp <- getrpcalls("C:/dev/neochromosome/mmc3.xlsx", "T1000 (DR)")
cgr <- getcgr("C:/dev/neochromosome/mmc4.xlsx", "T1000_CGRs")
cn <- getcn("C:/dev/neochromosome/mmc4.xlsx", "T1000_CN")
vcf <- readVcf("Z:/projects/liposarcoma/data/gridss/T1000/T1000.vcf", "hg19")
go(sample, vcf, rp, cgr, minimumEventSize=500)

sample <- "GOT3"
rp <- getrpcalls("C:/dev/neochromosome/mmc3.xlsx", "GOT3 (DR)")
cgr <- getcgr("C:/dev/neochromosome/mmc4.xlsx", "GOT3_CGRs")
cn <- getcn("C:/dev/neochromosome/mmc4.xlsx", "GOT3_CN")
vcf <- readVcf("W:/Papenfuss_lab/projects/liposarcoma/data/gridss/GOT3/GOT3.vcf", "hg19")
go(sample, vcf, rp, cgr, minimumEventSize=500)





# Working set
df <- out$gridss
df$size <- vcftobpgr(out$vcf)$size
ggplot(df[df$confidence=="High" & is.na(df$bedid), ]) +
  aes(x=size) +
  geom_histogram() +
  scale_x_log10()

table(df[is.na(df$bedid), ]$confidence)
table(df[df$confidence=="High" & is.na(df$bedid), ]$size <= 500, useNA="ifany")

