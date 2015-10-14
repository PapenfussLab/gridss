#source("http://bioconductor.org/biocLite.R")
#biocLite("rtracklayer")
library(ggplot2)
library(rtracklayer)
library(data.table)
library(stringr)
library(scales)
library(parallel)
library(foreach)
library(doParallel)
source("libgridss.R")
source("libneochromosome.R")
source("libvcf.R")

theme_set(theme_bw())

minsize <- 51 # minimum Mills 2012 event size
encodeblacklist <- import("consensusBlacklist.bed")
altassembly <- GRanges(seqnames=c("chr11_gl000202_random","chr17_gl000203_random","chr17_gl000204_random","chr17_gl000205_random","chr17_gl000206_random","chr18_gl000207_random","chr19_gl000208_random","chr19_gl000209_random","chr1_gl000191_random","chr1_gl000192_random","chr21_gl000210_random","chr4_gl000193_random","chr4_gl000194_random","chr7_gl000195_random","chr8_gl000196_random","chr8_gl000197_random","chr9_gl000198_random","chr9_gl000199_random","chr9_gl000200_random","chr9_gl000201_random","chrM","chrUn_gl000211","chrUn_gl000212","chrUn_gl000213","chrUn_gl000214","chrUn_gl000215","chrUn_gl000216","chrUn_gl000217","chrUn_gl000218","chrUn_gl000219","chrUn_gl000220","chrUn_gl000221","chrUn_gl000222","chrUn_gl000223","chrUn_gl000224","chrUn_gl000225","chrUn_gl000226","chrUn_gl000227","chrUn_gl000228","chrUn_gl000229","chrUn_gl000230","chrUn_gl000231","chrUn_gl000232","chrUn_gl000233","chrUn_gl000234","chrUn_gl000235","chrUn_gl000236","chrUn_gl000237","chrUn_gl000238","chrUn_gl000239","chrUn_gl000240","chrUn_gl000241","chrUn_gl000242","chrUn_gl000243","chrUn_gl000244","chrUn_gl000245","chrUn_gl000246","chrUn_gl000247","chrUn_gl000248","chrUn_gl000249"),
                       ranges=IRanges(start=0, end=1000000000))
altassembly$name <- as.character(seqnames(altassembly))
altassembly$score <- 1000
blacklist <- as(rbind(as(encodeblacklist, "RangedData"), as(altassembly, "RangedData")), "GRanges")


pwd <- getwd()
setwd(paste0(ifelse(as.character(Sys.info())[1] == "Windows", "W:/", "~/"), "i/data.na12878"))
metadata <- LoadMetadata()
vcfs <- LoadVcfs(metadata, existingVcfs=vcfs)
setwd(pwd)

truth <- "lumpyPacBioMoleculo"
vcfs <- lapply(vcfs, function(vcf) {
  attr(vcf, "metadata")$CX_REFERENCE_VCF <- "00000000000000000000000000000001.reference.vcf"
  return(vcf)
})
truth <- "Mills"
vcfs <- lapply(vcfs, function(vcf) {
  attr(vcf, "metadata")$CX_REFERENCE_VCF <- "00000000000000000000000000000002.reference.vcf"
  return(vcf)
})
truth <- "merged"
vcfs <- lapply(vcfs, function(vcf) {
  attr(vcf, "metadata")$CX_REFERENCE_VCF <- "00000000000000000000000000000000.reference.vcf"
  return(vcf)
})

# garvan <- readVcf("W:/na12878/garvan/garvan-na12878.vcf", "unknown")
# seqlevels(garvan) <- paste0("chr", seqlevels(garvan))
# attr(garvan, "metadata") <- attr(vcfs[["a136cfee2d40e62b4fd4be194366f291"]], "metadata")
# attr(garvan, "metadata")$CX_READ_LENGTH <- 150
# attr(garvan, "metadata")$CX_CALLER <- "gridss_garvan"
# attr(garvan, "metadata")$Id <- "10000000000000000000000000000000"
# attr(garvan, "metadata")$File <- "10000000000000000000000000000000"
# rownames(attr(garvan, "metadata")) <- "10000000000000000000000000000000"
# vcfs[["10000000000000000000000000000000"]] <- garvan

vcfs <- lapply(vcfs, function(vcf) {
  #vcf <- vcf[!isInterChromosmal(vcf),]
  vcf <- vcf[isDeletionLike(vcf, minsize), ]
  return(vcf)
})
vcfs <- lapply(vcfs, function(vcf) {
  withqual(vcf, attr(vcf, "metadata")$CX_CALLER)
})

truthlist_filtered <- CalculateTruthSummary(vcfs, blacklist=blacklist, maxerrorbp=100, maxerrorpercent=NULL, ignoreFilters=FALSE, ignore.strand=TRUE) # breakdancer does not specify strand
truthlist_all <- CalculateTruthSummary(vcfs, blacklist=blacklist, maxerrorbp=100, maxerrorpercent=NULL, ignoreFilters=TRUE, ignore.strand=TRUE) # breakdancer does not specify strand

dtroc_filtered <- TruthSummaryToROC(truthlist_filtered, bylist=c("CX_CALLER"))
dtroc_filtered$Filter <- "Default"
dtroc_all <- TruthSummaryToROC(truthlist_all, bylist=c("CX_CALLER"))
dtroc_all$Filter <- "Including Filtered"
dtroc <- rbind(dtroc_all, dtroc_filtered)

ggplot(dtroc) + 
  aes(y=tp/2, x=fp/2, color=CX_CALLER, linetype=Filter) +
  geom_line() + 
  scale_color_brewer(palette="Set1") + 
  scale_x_log10() + 
  scale_x_continuous(limits=c(0, 1000)) + 
  labs(y="tp", x="fp", title=paste("ROC curve NA12878 deletions", truth))
ggsave(paste0("na12878_roc_tp_fp_", truth, ".pdf"))

ggplot(dtroc) + 
  aes(y=prec, x=sens, color=CX_CALLER, linetype=Filter) +
  scale_color_brewer(palette="Set2") + 
  geom_line() + 
  labs(title=paste("Precision-Recall curve NA12878 deletions", truth))
ggsave(paste0("na12878_prec_recall_", truth, ".pdf"))

##########################
## GRIDSS single sample

millsgr <- bedpe2grmate(paste0(ifelse(as.character(Sys.info())[1] == "Windows", "W:/", "~/"), "na12878/lumpy-Mills2012-call-set.bedpe"))
pacbiogr <- bedpe2grmate(paste0(ifelse(as.character(Sys.info())[1] == "Windows", "W:/", "~/"), "na12878/lumpy-PacBioMoleculo-call-set.bedpe"))
mcols(millsgr) <- NULL
mcols(pacbiogr) <- NULL
write.csv(names(millsgr[overlapsAny(millsgr, pacbiogr),]), "W:/na12878/millsoverlappingpacbio.csv")

millsgr <- bedpe2grmate(paste0(ifelse(as.character(Sys.info())[1] == "Windows", "W:/", "~/"), "na12878/lumpy-Mills2012-call-set.bedpe"))
pacbiogr <- bedpe2grmate(paste0(ifelse(as.character(Sys.info())[1] == "Windows", "W:/", "~/"), "na12878/lumpy-PacBioMoleculo-call-set.bedpe"))
vcf <- readVcf("~/na12878/platinum-na12878.vcf", "hg19")
gvcf <- readVcf("~/garvan-na12878.vcf")

vcf <- vcf[rowRanges(vcf)$QUAL >= 100,]
vcf <- vcf[rowRanges(vcf)$FILTER %in% c(".", "PASS"),]
vcf <- vcf[gridss.vcftodf(vcf)$assembly %in% c("Both"), ] # high qual
vcf <- vcf[(str_detect(as.character(alt(vcf)), "\\[") & str_detect(names(rowRanges(vcf)), "o$")) |
             (str_detect(as.character(alt(vcf)), "]") & str_detect(names(rowRanges(vcf)), "h$")),]
# more bases deleted than added
vcf <- gridss.removeUnpartnerededBreakend(vcf)
vcf <- vcf[abs(start(rowRanges(vcf)) - start(rowRanges(vcf)[as.character(info(vcf)$MATEID),])) > nchar(gridss.vcftodf(vcf)$INSSEQ),]
vcf <- gridss.removeUnpartnerededBreakend(vcf)
vcf <- vcf[seqnames(rowRanges(vcf)) == seqnames(rowRanges(vcf)[as.character(info(vcf)$MATEID),]),] # intra-chromsomal events
vcf <- gridss.removeUnpartnerededBreakend(vcf)

#seqlevels(millsgr) <- sub("chr", "", seqlevels(millsgr))
mills <- gridss.annotateBreakpoints(millsgr, millsgr[millsgr$mate,], vcf, maxgap=100, ignore.strand=TRUE)
mills$bed$truth <- "Mills"
mills$gridss$truth <- "Mills"

#seqlevels(pacbiogr) <- sub("chr", "", seqlevels(pacbiogr))
pb <- gridss.annotateBreakpoints(pacbiogr, pacbiogr[pacbiogr$mate,], vcf, maxgap=100)
pb$bed$truth <- "PacBioOverlap"
pb$gridss$truth <- "PacBioOverlap"

mills$bed$V11 <- NULL
out <- list(bed=c(mills$bed, pb$bed), gridss=rbind(mills$gridss, pb$gridss))
out$bed$size <- abs(start(out$bed) - start(out$bed[out$bed$mate,]))

out$gridss <- out$gridss[order(-out$gridss$QUAL),]
roc <- data.table(out$gridss)[, list(cumtp=cumsum(!is.na(bedid)), cumfp=cumsum(is.na(bedid))), by=list(truth)]
ggplot(roc) + aes(x=cumfp/2, y=cumtp/2, color=truth) + geom_point() + geom_line()
ggsave("na12878_lumpy_fig3_roc.png")

sum(!is.na(out$bed$gridssid)) / length(out$bed) # sens
sum(!is.na(out$gridss$bedid)) / nrow(out$gridss) # prec
ggplot(as.data.frame(out$bed)) + aes(x=size, fill=ifelse(is.na(gridssid), "FN", "TP")) + geom_histogram() + scale_x_log10() + facet_grid(truth ~ . ) +
  labs(fill="", title="Gridss NA12878 Sensitivity")
ggsave("na12878_sens.png")
ggplot(as.data.frame(out$gridss)) + aes(x=size, fill=ifelse(is.na(bedid), "FP", "TP")) + geom_histogram() + scale_x_log10() + facet_grid(truth ~ . ) +
  labs(fill="", title="Gridss NA12878 Precision by event size")
ggsave("na12878_prec_size.png")
ggplot(as.data.frame(out$gridss)) + aes(x=QUAL, fill=!is.na(bedid)) + geom_histogram() + scale_x_log10() + facet_grid(truth ~ . ) +
  labs(fill="", title="Gridss NA12878 Precision by called qual")
ggsave("na12878_prec_qual.png")

ggplot(aes(data=as.data.frame(out$bed), x=size, fill=!is.na(bedid))) + geom_histogram() + scale_x_log10()


# buffer sizes
library(reshape2)
dtbuffer <- read.csv(file="~/na12878/platinum/gridss_vis/positional-chr12-Forward.csv", header=TRUE)
dtbuffers <- dtbuffer[,c("trackerActive","supportProcessedSize","aggregateProcessedSize","aggregateQueueSize","aggregateActiveSize","pathNodeProcessedSize","pathNodeActiveSize","pathNodeEdgeLookupSize","pathNodePathLookupSize","collapseProcessedSize","collapseUnprocessedSize","simplifyProcessedSize","simplifyLookupSize","simplifyUnprocessedSize","trackerLookupSize")]
#dtbuffers <- dtbuffers[sample(1:nrow(dtbuffers), 100000),]
dtbuffermelt <- melt(dtbuffers)
ggplot(dtbuffermelt, aes(x=variable, y=value)) + geom_boxplot() + scale_y_log10() + labs(x="Data Structure", y="Size") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("na12878_internal_buffer_size.png")
dtops <- dtbuffer[,c("contigSize","contigFrontierSize","contigMemoizedSize","contigUnprocessedSize","assemblyActiveSize")]
#dtops <- dtops[sample(1:nrow(dtops), 100000),]
dtopmelt <- melt(dtops)
ggplot(dtopmelt, aes(x=variable, y=value)) + geom_boxplot() + scale_y_log10() + labs(x="", y="") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("na12878_internal_assembly_size.png")
ggplot(rbind(dtopmelt,dtbuffermelt), aes(x=variable, y=value)) + geom_boxplot() + scale_y_log10() + labs(x="", y="") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("na12878_assembly_buffers.png")

mean(dtbuffer$trackerActive)
mean(dtbuffer$trackerLookupSize)

