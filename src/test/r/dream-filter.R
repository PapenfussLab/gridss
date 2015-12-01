#source("http://bioconductor.org/biocLite.R")
#biocLite("VariantAnnotation")
library(VariantAnnotation)
setwd("W:/dev/gridss/src/test/r/")
source("libgridss.R")
source("libvcf.R")

primaryContigs <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")

file <- "W:/dream/synthetic4/synthetic4.vcf"
vcf <- readVcf(file, "hg19_random")
vcf <- gridss.removeUnpartnerededBreakend(vcf)
df <- gridss.vcftodf(vcf, allColumns=TRUE)
df$refAF <- (df$REF1+df$REFPAIR1) / (df$SR1+df$RP1+df$RSR1)
matevcf <- vcf[df$mate,]

isSimpleEvent <- as.logical(seqnames(rowRanges(vcf)) == seqnames(rowRanges(matevcf)) & 
  pmax(start(rowRanges(vcf)), start(rowRanges(matevcf))) - pmin(start(rowRanges(vcf)), start(rowRanges(matevcf))) >= 200 &
  pmax(start(rowRanges(vcf)), start(rowRanges(matevcf))) - pmin(start(rowRanges(vcf)), start(rowRanges(matevcf))) <= 1002000)
isPrimary <- as.logical(seqnames(rowRanges(vcf)) %in% primaryContigs & seqnames(rowRanges(matevcf)) %in% primaryContigs)
hasNormalSupport <- df$Q0 > 0
hasLongHomology <- df$HOMLEN >= 10
hasNormalBreakendSupport <- df$BUM0 + df$BSC0 >= 2
hasNoDirectSupport <- is.nan(df$refAF) | is.infinite(df$refAF)
hasSuspiciousFilter <- df$FILTER %in% c("ASSEMBLY_ONLY;LOW_QUAL", "SINGLE_ASSEMBLY")
isUnexpectedSize <- df$size < 400 | df$size > 40000
hasHighRefCoverage <- df$REF1 + df$REFPAIR1 >= 200
isNotPrimaryCall <- df$QUAL < 0.8 * df$CQ

for (qual in c(0, 100, 150, 175, 200, 225, 250, 300, 400, 500, 600, 700, 800, 900, 1000)) {
  vcf2 <- vcf[isPrimary &
                isSimpleEvent &
                !hasNormalSupport &
                !hasLongHomology &
                !hasNormalBreakendSupport &
                !hasNoDirectSupport &
                !hasSuspiciousFilter &
                !isUnexpectedSize &
                !hasHighRefCoverage &
                !isNotPrimaryCall &
                df$QUAL>=qual,]
  rowRanges(vcf2)$FILTER <- "."
  vcf2 <- gridss.removeUnpartnerededBreakend(vcf2)
  writeVcf(vcf2, paste0("W:/dream/synth4-gridss0.9.2-q", qual, ".vcf"))
}


vcf <- vcf[isSimpleEvent & isPrimary,]
library(ggplot2)
mvcf <- readVcf("W:/dream/synthetic.challenge.set4.tumour.25pctmasked.truth.vcf", "hg19_random")
mvcf <- mvcf[!is.na(info(mvcf)$SVTYPE) & info(mvcf)$SVTYPE == "MSK",]
vcf <- vcf[!overlapsAny(vcf, mvcf, type="any")] # use the unmasked subset so we're not cheating by tuning to the full truth set
vcf <- gridss.removeUnpartnerededBreakend(vcf)
df <- gridss.vcftodf(vcf, allColumns=TRUE)
matevcf <- vcf[df$mate,]
tvcf <- readVcf("W:/dream/truthsv4.vcf", "hg19_random")
tgr <- GRanges(seqnames=seqnames(tvcf), ranges=IRanges(start=start(tvcf), end=info(tvcf)$END))
gr <- GRanges(seqnames=seqnames(vcf), ranges=IRanges(start=pmin(start(vcf), start(matevcf)), end=pmax(start(vcf), start(matevcf))))
shits <- findOverlaps(gr, tgr, type="start", maxgap=75)
ehits <- findOverlaps(gr, tgr, type="end", maxgap=75)
hits <- data.table(queryHits=c(queryHits(shits), queryHits(ehits)), subjectHits=c(subjectHits(shits), subjectHits(ehits)))
hits <- hits[duplicated(hits)]
gr$truth <- FALSE
gr$truth[hits$queryHits] <- TRUE
df$truth <- gr$truth
df <- data.table(df)
ggplot(df) + aes(x=Q0 + 1, fill=truth) + geom_histogram() + scale_x_log10()
df <- df[df$Q0==0,]
ggplot(df[,list(prec=sum(truth)/.N), by="HOMLEN"]) + aes(x=HOMLEN, y=prec) + geom_line() + scale_x_continuous(limits=c(1, 50))
df <- df[df$HOMLEN<10,]
ggplot(df) + aes(x=QUAL, color=truth) + geom_histogram() + scale_x_log10() + facet_wrap( ~ FILTER)
df <- df[!(df$FILTER %in% c("ASSEMBLY_ONLY;LOW_QUAL", "SINGLE_ASSEMBLY")),]
ggplot(df) + aes(x=BUM0+BSC0, fill=truth) + geom_histogram(binwidth=1) + scale_x_continuous(limits=c(1, 10))
df <- df[df$BUM0+df$BSC0 < 2,]
table(is.nan(df$refAF) | is.infinite(df$refAF), ifelse(df$truth, "tp", "fp"))
df <- df[!(is.nan(df$refAF) | is.infinite(df$refAF)),]
ggplot(df) + aes(x=size, color=truth) + geom_histogram(binwidth=200) + scale_x_continuous(limits=c(100, 100000))
df <- df[df$size >= 400 & df$size <= 40000,]
ggplot(df) + aes(x=REF1+REFPAIR1, fill=truth) + geom_histogram() + scale_x_log10()
table(df$REF1+df$REFPAIR1 > 200, ifelse(df$truth, "tp", "fp"))
df <- df[df$REF1+df$REFPAIR1 < 200,]



ggplot(df) + aes(x=QUAL, fill=truth) + geom_histogram() + scale_x_log10()



df$refreadAF <- (df$REF1) / (df$SR1+df$RSR1)
df$refpairAF <- (df$REFPAIR1) / (df$RP1)
ggplot(df) + aes(x=refAF, fill=truth) + geom_histogram(binwidth=1)+ scale_x_continuous(limits=c(0, 50)) 
ggplot(df) + aes(x=refreadAF, fill=truth) + geom_histogram(binwidth=1) + scale_x_continuous(limits=c(0, 100))
ggplot(df) + aes(x=refpairAF, fill=truth) + geom_histogram(binwidth=1) + scale_x_continuous(limits=c(0, 50))



library(ROCR) # install.packages("ROCR")
library(arm)
modeldf <- df[paste0(df$EVENT, "o")==df$variantid,]
modeldf$logsize <- log(modeldf$size)
modeldf$nlogsize <- abs((modeldf$logsize - mean(modeldf$logsize)) / sd(modeldf$logsize))
modeldf$logQUAL <- log(modeldf$QUAL)
modeldf$nlogQUAL <- (modeldf$QUAL - mean(modeldf$QUAL)) / sd(modeldf$QUAL)
traindf <- modeldf[seq.int(1, nrow(modeldf), 2),]
testdf <- modeldf[seq.int(2, nrow(modeldf), 2),]
baseline <- NULL
for (modelpara in c(
  truth ~ nlogQUAL,
  truth ~ nlogQUAL + HOMLEN,
  truth ~ nlogQUAL + nlogsize,
  truth ~ nlogQUAL + refAF,
  truth ~ nlogQUAL + nlogsize + refAF)) {
  model <- glm(formula=modelpara, family=binomial(link='logit'), data=traindf)
  testdf$response <- predict(model, newdata=testdf, type='response')
  plot(performance(prediction(testdf$response, testdf$truth), measure = "tpr", x.measure = "fpr"))
  auc <- performance(prediction(testdf$response, testdf$truth), measure = "auc")@y.values[[1]]
  if (is.null(baseline)) { baseline <- auc }
  print(paste(paste(as.character(modelpara), collapse=" "), auc, (auc - baseline)*100))
}
  




