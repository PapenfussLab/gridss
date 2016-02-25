##################
# Processing steps (approximate runtime: 400h wall, 4500h CPU)
##################
# - Follow common.R processing steps
# - Download ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/integrated_sv_map/supporting/NA12878 to ~/projects/reference_datasets/human_sequencing/NA12878/ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/integrated_sv_map/supporting/NA12878
# - Download ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ to ~/projects/reference_datasets/human_sequencing/NA12878/ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/
# - symbolic link ~/reference_genomes -> ~/projects/reference_genomes/
# - Run src/test/r/sim/longread/align.sh
# - Run src/test/r/sim/longread/tobed.sh
# - Run src/test/r/sim/platinum/align.sh
# - symbolic link ERA172924.sc.bam to data.na12878/00000000000000000000000000000000.sc.bam
# - symbolic link ERR194147.sorted.bam to data.na12878/00000000000000000000000000000000.sc.bam.bt2.bam
# - symbolic link ERR194147.sc.bam to data.na12878/00000000000000000000000000000000.sc.sr.bam
# - ./call_gridss.sh na12878
# - ./call_breakdancer.sh na12878
# - ./call_delly.sh na12878
# - ./call_pindel.sh na12878
# - ./call_socrates.sh na12878
# - ./call_lumpy.sh na12878
# - Run this script (R requires at least 6GB of memory for this script)

#source("http://bioconductor.org/biocLite.R")
#biocLite("rtracklayer")
source("common.R")
source("libna12878.R")

vcfs <- NULL
pwd <- getwd()
setwd(paste0(rootdir, "i/data.na12878"))
metadata <- LoadMetadata()
metadata <- main_callers_subset(metadata)
vcfs <- LoadVcfs(metadata, existingVcfs=vcfs)
setwd(pwd)

# truth <- "lumpyPacBioMoleculo"
# vcfs <- lapply(vcfs, function(vcf) {
#   attr(vcf, "sourceMetadata")$CX_REFERENCE_VCF <- "00000000000000000000000000000001.reference.vcf"
#   return(vcf)
# })
# truth <- "merged"
# vcfs <- lapply(vcfs, function(vcf) {
#   attr(vcf, "sourceMetadata")$CX_REFERENCE_VCF <- "00000000000000000000000000000000.reference.vcf"
#   return(vcf)
# })
truth <- "Mills"
vcfs <- lapply(vcfs, function(vcf) {
  attr(vcf, "sourceMetadata")$CX_REFERENCE_VCF <- "00000000000000000000000000000002.reference.vcf"
  return(vcf)
})
vcfs <- lapply(vcfs, function(vcf) {
  #vcf <- vcf[!isInterChromosmal(vcf),]
  write(paste0("Processing ", attr(vcf, "sourceMetadata")$Id), stderr())
  vcf <- vcf[isDeletionLike(vcf, minsize), ]
  caller <- str_extract(attr(vcf, "sourceMetadata")$CX_CALLER, "^[^/]+")
  if (!is.na(caller) && !is.null(caller) && caller %in% c("breakdancer")) {
    # strip all BND events since they are non-deletion events such as translocations
    vcf <- vcf[info(vcf)$SVTYPE == "DEL",]
  }
  if (!is.na(caller) && !is.null(caller)) {
    vcf <- withqual(vcf, caller)
  }
  return(vcf)
})

truthlist_filtered <- CalculateTruthSummary(vcfs, blacklist=blacklist, maxerrorbp=maxPositionDifference_LengthOrEndpoint, maxerrorpercent=maxLengthRelativeDifference, ignoreFilters=FALSE, ignore.strand=TRUE) # breakdancer does not specify strand
truthlist_all <- CalculateTruthSummary(vcfs, blacklist=blacklist, maxerrorbp=maxPositionDifference_LengthOrEndpoint, maxerrorpercent=maxLengthRelativeDifference, ignoreFilters=TRUE, ignore.strand=TRUE) # breakdancer does not specify strand

truthlist_filtered$calls <- truthlist_filtered$calls
truthlist_filtered$truth <- truthlist_filtered$truth
truthlist_all$calls <- truthlist_all$calls
truthlist_all$truth <- truthlist_all$truth

dtroc_filtered <- TruthSummaryToROC(truthlist_filtered, bylist=c("CX_CALLER"))
dtroc_filtered$Filter <- prettyFilter(text_default_calls)
dtroc_all <- TruthSummaryToROC(truthlist_all, bylist=c("CX_CALLER"))
dtroc_all$Filter <-prettyFilter(text_all_calls)
dtroc <- rbind(dtroc_all, dtroc_filtered)

ggplot(dtroc) + 
  aes(y=tp/2, x=fp/2, color=CX_CALLER, linetype=Filter) +
  geom_line() + 
  geom_line(data=dtroc[dtroc$CX_CALLER=="gridss",], size=1.2) +
#  scale_color_brewer(palette="Set2") + 
  coord_cartesian(xlim=c(0, 1000), ylim=c(0, 3000)) + 
  labs(y="tp", x="fp", title=paste("ROC curve NA12878 deletions", truth))
saveplot(paste0("na12878_tp_fp_", truth, "_error_", maxLengthRelativeDifference, "_", maxPositionDifference_LengthOrEndpoint, ""), width=7, height=5)
ggplot(dtroc) + 
  aes(y=prec, x=tp/2, color=CX_CALLER, linetype=Filter) +
#  scale_color_brewer(palette="Set2") + 
  coord_cartesian(xlim=c(0, 4000)) +
  geom_line() + 
  geom_line(data=dtroc[dtroc$CX_CALLER=="gridss",], size=1.2) +
  labs(y="Precision", x="True positives", title=paste("Precision-Recall curve NA12878 deletions", truth))
saveplot(paste0("na12878_prec_", truth, "_error_", maxLengthRelativeDifference, "_", maxPositionDifference_LengthOrEndpoint, ""), width=7, height=5)

dtroc[dtroc$prec > 0.75,][!duplicated(paste(dtroc[dtroc$prec > 0.75,]$CX_CALLER, dtroc[dtroc$prec > 0.75,]$Filter), fromLast=TRUE),]

###########
# Moleculo/PacBio truth set
lrbedraw <- longReadBed(vcfs)
lrbed <- lrbedraw
lrbed <- lrbed[order(-lrbed$tp, -lrbed$QUAL),] # remove the FP / lower scoring duplicate
lrbed <- lrbed[!isDuplicateCall(lrbed, maxgap=50),]
delroc <- bedToROC(lrbed)
delroc$caller <- as.factor.gridss.first(delroc$caller)
ggplot(delroc) +
  aes(x=fp, y=tp, color=caller, linetype=Filter) +
  geom_line() +
  geom_line(delroc[delroc$caller=="gridss",], size=1.2) +
  coord_cartesian(xlim=c(0, 1000), ylim=c(0, 3000)) + 
  labs(x="False Positives", y="True Positives", title="NA12878 deletions PacBio/Moleculo validated")
#  scale_color_brewer(palette="Set2") +
saveplot(paste0("na12878_tp_fp_pacbiomoleculo", "_error_", maxLengthRelativeDifference, "_", maxPositionDifference_LengthOrEndpoint, ""), width=7, height=5)
ggplot(delroc) +
  aes(x=tp, y=precision, color=caller, linetype=Filter) +
  geom_line() +
  geom_line(delroc[delroc$caller=="gridss",], size=1.2) +
  coord_cartesian(xlim=c(0, 4000)) + 
  labs(x="True Positives", y="Precision", title="NA12878 deletions PacBio/Moleculo validated")
#  scale_color_brewer(palette="Set2") +  
saveplot(paste0("na12878_prec_pacbiomoleculo", "_error_", maxLengthRelativeDifference, "_", maxPositionDifference_LengthOrEndpoint, ""), width=7, height=5)

###########
# precision thresholds
delroc[delroc$precision > 0.95,][!duplicated(paste(delroc[delroc$precision > 0.95,]$caller, delroc[delroc$precision > 0.95,]$Filter), fromLast=TRUE),]
ggplot(as.data.frame(lrbed[lrbed$caller=="gridss"])) +
  aes(x=QUAL, fill=ifelse(tp, "_tp", "fp")) +
  facet_wrap(~ filtered) +
  geom_histogram(binwidth=100) +
  coord_cartesian(xlim=c(0, 3000), ylim=c(0, 200))
saveplot("na12878_gridss_hq_treshold")

# gridss FPR by QUAL
gridssvcf_binned <- lapply(vcfs, function(vcf) {
  if (is.na(attr(vcf, "sourceMetadata")$CX_CALLER)) return(NULL)
  if (str_detect(attr(vcf, "sourceMetadata")$CX_CALLER, "gridss")) {
    rowRanges(vcf)$QUAL <- floor(rowRanges(vcf)$QUAL / 100) * 100
    return (vcf)
  }
  return(NULL)
})
gridssvcf_binned <- gridssvcf_binned[!sapply(gridssvcf_binned, is.null)]
expect_that(length(gridssvcf_binned), equals(1))
bindelbed <- longReadBed(gridssvcf_binned)
bindelroc <- bedToROC(bindelbed)
bindelroc <- bindelroc[bindelroc$Filter == text_default_calls,]
bindelroclookup <- bindelroc
bindelroclookup[nrow(bindelroclookup)]$QUAL <- 1000000
bindelroclookup <- bindelroclookup[order(-bindelroclookup$QUAL),]
bindelroc$dtp <- bindelroc$tp - bindelroclookup$tp
bindelroc$dfp <- bindelroc$fp - bindelroclookup$fp
bindelroc$tpr <- bindelroc$dtp / (bindelroc$dtp  + bindelroc$dfp )
ggplot(bindelroc) + aes(y=tpr, x=QUAL) + geom_line() + scale_x_continuous(limits=c(500, 5000)) +
  labs(title="Gridss TPR by QUAL score", y="True Positive Rate", x="QUAL score")
saveplot("na12878_gridss_tpr_by_qual")

###########
# Logicistic regression of gridss calls
# http://www.r-bloggers.com/how-to-perform-a-logistic-regression-in-r/
library(ROCR) # install.packages("ROCR")
library(arm)
gridssvcf <- lapply(vcfs, function(vcf) {
  if (is.na(attr(vcf, "sourceMetadata")$CX_CALLER)) return(NULL)
  if (attr(vcf, "sourceMetadata")$CX_CALLER == "gridss/0.10.0-SNAPSHOT") {
    return (vcf)
  }
  return(NULL)
})
gridssvcf <- gridssvcf[!sapply(gridssvcf, is.null)]
expect_that(length(gridssvcf), equals(1))
gridsslrbed <- longReadBed(gridssvcf)
gridssdf <- gridss.vcftodf(gridssvcf[[1]])
gridssdf <- gridssdf[as.character(gridsslrbed$id),]
gridssdf$blacklisted <- gridsslrbed$blacklisted
gridssdf$LongReadSupport <- gridsslrbed$tp
gridssdf$rp <- gridssdf$RP > 0
gridssdf$sr <- gridssdf$SR + gridssdf$RSR > 0
gridssdf$AQ <- gridssdf$ASQ + gridssdf$RASQ
gridssdf$SQ <- gridssdf$SRQ + gridssdf$RSRQ
gridssdf$LOW_QUAL <- gridssdf$QUAL < 1000
gridssdf$SPV <- NULL
modeldf <- gridssdf[,!(names(gridssdf) %in% c("variantid","POS","FILTER","EVENT","mate","SOMATIC","SVTYPE","HOMSEQ", "INSSEQ", "confidence"))]
#modeldf <- modeldf[modeldf$QUAL > 1000,]
for (col in c("QUAL","HOMLEN","REF","REFPAIR","SPV","CQ","BQ","AS","RP","SR","RAS","RSR","ASRP","ASSR","ASCRP","ASCSR","ASQ","RPQ","SRQ","RASQ","RSRQ","BA","BUM","BSC","BAQ","BUMQ","BSCQ","size")) {
  modeldf[col] <- (modeldf[col] - mean(unlist(modeldf[col]))) / sd(unlist(modeldf[col]))
}
traindf <- modeldf[seq.int(1, nrow(modeldf), 2),]
testdf <- modeldf[seq.int(2, nrow(modeldf), 2),]
# Try a bunch of models
baseline <- NULL
for (modelpara in c(
  LongReadSupport ~ QUAL, # 0.968486742605028 baseline area under curve # 50% QUAL threshold is 1147
  LongReadSupport ~ ASQ + RASQ,  # worse since no low quality tail
  #LongReadSupport ~ QUAL + blacklisted, # not meaningful as the blacklisted filter was already applied at runtime
  LongReadSupport ~ CQ,  # strangely better than QUAL
  LongReadSupport ~ assembly, # (only two points)
  LongReadSupport ~ QUAL + assembly,
  LongReadSupport ~ AQ + RPQ + SQ,
  LongReadSupport ~ AQ + RPQ + SQ + BUMQ + BSCQ,
  LongReadSupport ~ RP,
  LongReadSupport ~ SR + RSR,
  LongReadSupport ~ RP + SR + RSR,
  LongReadSupport ~ ASSR + ASRP,
  LongReadSupport ~ ASSR + ASRP + ASCRP + ASCSR,
  LongReadSupport ~ ASSR + ASRP + RP + SR + RSR + BUM + BSC + ASCRP + ASCSR,
  LongReadSupport ~ QUAL + BUMQ + BSCQ + rp + sr + rp + assembly,
  LongReadSupport ~ LOW_QUAL * .,
  LongReadSupport ~ ., # prediction from a rank-deficient fit may be misleading
  LongReadSupport ~ QUAL + HOMLEN + REF + REFPAIR + AS + RP + SR + ASRP + ASSR + ASCRP + ASCSR + ASQ + RPQ + SQ + BUM + BUMQ + BSC + BSCQ + rp + sr + assembly,
  LongReadSupport ~ sr * QUAL,
  LongReadSupport ~ sr * (AQ + RPQ + SQ),
  LongReadSupport ~ (assembly + sr) * (AQ + RPQ + SQ),
  LongReadSupport ~ QUAL + assembly + sr + (QUAL + assembly + sr)^2
  )) { 
  model <- glm(formula=modelpara, family=binomial(link='logit'), data=traindf)
  testdf$response <- predict(model, newdata=testdf, type='response')
  plot(performance(prediction(testdf$response, testdf$LongReadSupport), measure = "tpr", x.measure = "fpr"))
  auc <- performance(prediction(testdf$response, testdf$LongReadSupport), measure = "auc")@y.values[[1]]
  if (is.null(baseline)) { baseline <- auc }
  print(paste(paste(as.character(modelpara), collapse=" "), auc, auc - baseline))
}
coefplot(model)

#save.image("~/na12878.RData")
#load(file="~/na12878.RData")

