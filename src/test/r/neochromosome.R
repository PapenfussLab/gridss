#source("http://bioconductor.org/biocLite.R")
#biocLite("rtracklayer")
library(rtracklayer)
source("../../main/r/libgridss.R")
source("libneochromosome.R")
source("libvcf.R")

#setwd("W:/dev/gridss/src/test/R/")

theme_set(theme_bw())

sample <- "778"
rp <- getrpcalls("C:/dev/neochromosome/mmc3.xlsx", "778 (DR)")
# no need to explicitly remove these calls as they're not in CGRs
#rp <- rp[!(seqnames(rp) == seqnames(rp[rp$mate,]) & seqnames(rp) %in% c("chr7", "chr22")),] # remove 7 and 22 intrachromosomal arrangements that had really low limits
cgr <- getcgr("C:/dev/neochromosome/mmc4.xlsx", "778_CGRs")
cn <- getcn("C:/dev/neochromosome/mmc4.xlsx", "778_CN")
vcf <- readVcf("W:/Papenfuss_lab/projects/liposarcoma/data/gridss/778/778.vcf", "hg19")
out <- go(sample, vcf, rp, cgr, minimumEventSize=500)

sample <- "T1000"
rp <- getrpcalls("C:/dev/neochromosome/mmc3.xlsx", "T1000 (DR)")
cgr <- getcgr("C:/dev/neochromosome/mmc4.xlsx", "T1000_CGRs")
cn <- getcn("C:/dev/neochromosome/mmc4.xlsx", "T1000_CN")
vcf <- readVcf("Z:/projects/liposarcoma/data/gridss/T1000/T1000.vcf", "hg19")
out <- go(sample, vcf, rp, cgr, minimumEventSize=500)

sample <- "GOT3"
rp <- getrpcalls("C:/dev/neochromosome/mmc3.xlsx", "GOT3 (DR)")
cgr <- getcgr("C:/dev/neochromosome/mmc4.xlsx", "GOT3_CGRs")
cn <- getcn("C:/dev/neochromosome/mmc4.xlsx", "GOT3_CN")
vcf <- readVcf("W:/Papenfuss_lab/projects/liposarcoma/data/gridss/GOT3/bt2/GOT3-bt2.vcf", "hg19")
out <- go(sample, vcf, rp, cgr, minimumEventSize=500)

  

#####################
# realignment rates


#####################
# Working set
#
df <- out$gridss
df$size <- vcftobpgr(out$vcf)$size
ggplot(df) +
  aes(x=distanceToCall, y=QUAL, color=confidence) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap( ~ found)
table(out$bed$hits)
ggplot(as.data.frame(out$bed)) +
  aes(x=distanceToMedHigh, color=is.na(gridssid)) +
  geom_histogram() +
  scale_x_log10() +
  facet_wrap( ~ hits)

out$bed[out$bed$hits=="None",]

table(df[is.na(df$bedid), ]$confidence)
table(df[df$confidence=="High" & is.na(df$bedid), ]$size <= 500, useNA="ifany")

distanceToClosest(rp, rowRanges(vcf))
distanceToClosest(rowRanges(vcf), rp)
rp$distanceToCall
rp[distanceHits]

# Breakdown of small fragment spanning calls by call type
table(out$bed$spanning, out$bed$hits)
print("gridss found ", FindFragments(vcftobpgr(vcf)), "fragments")


ggplot(out$gridss[out$gridss$confidence %in% c("High", "Medium") & is.na(out$gridss$bedid),]) +
  aes(x=RP, color=confidence) + 
  geom_histogram() +
  labs(title="Read Pair support for unmatched gridss calls")
ggsave("rpsupport.png")
  


# Spanning
spanning <- FindFragmentSpanningEvents(vcftobpgr(out$vcf), vcftobpgr(out$vcf))
spanning <- spanning[substr(as.character(spanning$query), 1, nchar(as.character(spanning$query)) - 1) != substr(as.character(spanning$subjectAlt1), 1, nchar(as.character(spanning$subjectAlt1)) - 1) &
                       substr(as.character(spanning$query), 1, nchar(as.character(spanning$query)) - 1) != substr(as.character(spanning$subjectAlt2), 1, nchar(as.character(spanning$subjectAlt2)) - 1),]
table(out$gridss[as.character(spanning$query),]$confidence)
