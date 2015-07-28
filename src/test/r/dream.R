library(ggplot2)
library(scales)
source("../../main/r/libgridss.R")
source("libvcf.R")

#rawvcf <- readVcf("W:/dream/synthetic4/synthetic4.vcf", "hg19")
rawvcf <- readVcf("W:/dream/subset5/subset5.vcf", "hg19")

vcf <- rawvcf
vcf <- gridss.removeUnpartnerededBreakend(vcf)

df <- gridss.vcftodf(vcf, allColumns=TRUE)
###################
# Aggregate counts
mate <- df[df$mate,]
df$VariantBreakpointCount1 <- df$RP1 + df$SR1 + df$RSR1 + df$ASCRP1 + df$ASCSR1
df$VariantBreakpointCount2 <- df$RP2 + df$SR2 + df$RSR2 + df$ASCRP2 + df$ASCSR2
df$VariantReadCount1 <- df$BSC1 + mate$BSC1 + df$SR1 + mate$SR1 + df$ASCSR1
df$VariantPairCount1 <- df$BUM1 + mate$BUM1 + df$RP + df$ASCRP1
df$VariantCount1 <- df$VariantReadCount1 + df$VariantPairCount1
df$ReferenceReadCount1 <- df$REF1 + mate$REF1
df$ReferencePairCount1 <- df$REFPAIR1 + mate$REFPAIR1
df$ReferenceCount1 <- df$ReferenceReadCount1 + df$ReferencePairCount1
df$VariantReadCount2 <- df$BSC2 + mate$BSC2 + df$SR2 + mate$SR2 + df$ASCSR2
df$VariantPairCount2 <- df$BUM2 + mate$BUM2 + df$RP + df$ASCRP2
df$VariantCount2 <- df$VariantReadCount2 + df$VariantPairCount2
df$ReferenceReadCount2 <- df$REF2 + mate$REF2
df$ReferencePairCount2 <- df$REFPAIR2 + mate$REFPAIR2
df$ReferenceCount2 <- df$ReferenceReadCount2 + df$ReferencePairCount2

###############
# Challenge 5
# DREAM challenge breakpoint filters
vcf <- vcf[
  df$VariantBreakpointCount1==0 # no support in normal
  & (df$VariantBreakpointCount2 > 0 | df$AS == 0 | (df$BUM1 == 0 & df$BSC1==0)) # only breakend assembly support -> breakend support must not contain normal
  & !is.na(df$size) # must be intrachromosomal
  & df$size > 200 # not small INDEL
  & df$size < 1000000 # not large event
  & df$HOMLEN < 16 # alignment homologies
  ,]
vcf <- gridss.removeUnpartnerededBreakend(vcf)
df <- gridss.vcftodf(vcf, allColumns=TRUE)

writeVcf(vcf, "W:/dream/out.vcf")
# Now run BreakendToSimpleCall on out.vcf to convert to DREAM-friendly format
simplevcf <- readVcf("W:/dream/out.vcf.simple.vcf", "hg19")
sdf <- df[rownames(simplevcf),]
overlaps <- findOverlaps(rowRanges(simplevcf), rowRanges(simplevcf), maxgap=10000)
overlaps <- overlaps[queryHits(overlaps) != subjectHits(overlaps),]
overlaps <- overlaps[sdf$QUAL[queryHits(overlaps)] < sdf$QUAL[subjectHits(overlaps)]]
# DREAM challenge overlap filter
# Remove to lower quality call of overlapping variants
simplevcf <- simplevcf[-queryHits(overlaps),]
sdf <- df[rownames(simplevcf),]
# make dream submission script happy with format
rowRanges(simplevcf)$FILTER <- "."
info(simplevcf)$SOMATIC <- TRUE
info(simplevcf)$CIEND <- info(simplevcf)$CIRPOS
info(simplevcf)$IMPRECISE[sapply(info(simplevcf)$CIPOS, function(x) is.na(x)[1])] <- FALSE
writeVcf(simplevcf, "W:/dream/somatic5_all.vcf")


###################
# Truth evaluation
#
#truthvcf <- readVcf("W:/dream/synthetic.challenge.set1.tumor.all.truth.vcf", "hg19")
truthvcf <- readVcf("W:/dream/truthsv4.vcf", "hg19")
truthvcf <- truthvcf[!is.na(info(truthvcf)$SVTYPE) & unlist(rowRanges(truthvcf)$ALT) %in% c("<DEL>", "<DUP>", "<INV>", "<INS>"),]
truth <- rowRanges(truthvcf)
end(truth) <- info(truthvcf)$END
# hacked evaluator.py to output classification per call
#./evaluator.py ../synthetic1/synthetic1-rp600.vcf ../synthetic.challenge.set1.tumor.all.truth.vcf.gz SV | grep gridss > synth1.truth.csv
#evalmatch <- read.csv("W:/dream/tools/synth1.truth.csv")

df$truth <- NA
df$truth[queryHits(findOverlaps(rowRanges(vcf), truth, maxgap=100))] <- names(truth[subjectHits(findOverlaps(rowRanges(vcf), truth, maxgap=100))])
df$distance <- pmin(
  distanceToClosest(rowRanges(vcf), GRanges(seqnames=seqnames(truth), ranges=IRanges(start=start(truth), end=start(truth)))),
  distanceToClosest(rowRanges(vcf), GRanges(seqnames=seqnames(truth), ranges=IRanges(start=end(truth), end=end(truth)))))
#df$isExpectedDistance = df$distance >= 590 & df$distance <= 601
mate <- df[df$mate,]
df$match <- !is.na(df$distance) & !is.na(mate$distance) & df$distance < 100 & mate$distance < 100
dfmatch <- df[df$match,]
dfmatch <- dfmatch[order(dfmatch$QUAL),]
truth$besthit <- NA_character_
truth[dfmatch$truth]$besthit <- rownames(dfmatch)





