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

pwd <- getwd()
setwd(paste0(ifelse(as.character(Sys.info())[1] == "Windows", "W:/", "~/"), "i/data.na12878"))
metadata <- LoadMetadata()
vcfs <- LoadVcfs(metadata)
setwd(pwd)

vcfs <- lapply(vcfs, function(vcf) {
  attr(vcf, "metadata")$CX_REFERENCE_VCF <- "00000000000000000000000000000000.reference.vcf"
  vcf <- TrimInterChromosmalEvents(vcf)
  return(vcf)
})
garvan <- readVcf("W:/na12878/garvan/garvan-na12878.vcf", "unknown")
seqlevels(garvan) <- paste0("chr", seqlevels(garvan))
attr(garvan, "metadata") <- attr(vcfs[["a136cfee2d40e62b4fd4be194366f291"]], "metadata")
attr(garvan, "metadata")$CX_READ_LENGTH <- 150
attr(garvan, "metadata")$CX_CALLER <- "gridss_garvan"
attr(garvan, "metadata")$Id <- "10000000000000000000000000000000"
attr(garvan, "metadata")$File <- "10000000000000000000000000000000"
rownames(attr(garvan, "metadata")) <- "10000000000000000000000000000000"
vcfs[["10000000000000000000000000000000"]] <- garvan

vcfs <- lapply(vcfs, function(vcf) {
  caller <- attr(vcf, "metadata")$CX_CALLER
  if (all(is.na(rowRanges(vcf)$QUAL)) && !is.na(caller)) {
    # use total read support as a qual proxy
    if (caller %in% c("delly", "delly/0.6.8")) {
      rowRanges(vcf)$QUAL <- ifelse(is.na(info(vcf)$PE), 0, info(vcf)$PE) + ifelse(is.na(info(vcf)$SR), 0, info(vcf)$SR)
    } else if (caller %in% c("crest")) {
      rowRanges(vcf)$QUAL <- ifelse(is.na(info(vcf)$right_softclipped_read_count), 0, info(vcf)$right_softclipped_read_count) + ifelse(is.na(info(vcf)$left_softclipped_read_count), 0, info(vcf)$left_softclipped_read_count)
    } else if (caller %in% c("pindel")) {
      rowRanges(vcf)$QUAL <- geno(vcf)$AD[,1,2]
    } else {
      warning(paste("No QUAL scores for ", caller))
    }
  }
  return(vcf)
})

truthlist_all <- CalculateTruthSummary(vcfs, maxerrorbp=100, maxerrorpercent=NULL, ignoreFilters=FALSE, ignore.strand=TRUE) # breakdancer does not specify strand

dtroc <- TruthSummaryToROC(truthlist_all, bylist=c("CX_CALLER"))

ggplot(dtroc) + 
  aes(y=tp/2, x=fp/2, color=CX_CALLER) +
  geom_line() + 
  geom_point(size=0.25) +
  scale_x_log10() + 
  scale_x_continuous(limits=c(0, 1000)) + 
  labs(y="tp", x="fp", title="All intrachromosomal calls vs Mills deletion calls")


##########################
## GRIDSS single sample

millsgr <- bedpe2grmate("~/na12878/lumpy-Mills2012-call-set.bedpe")
pacbiogr <- bedpe2grmate("~/na12878/lumpy-PacBioMoleculo-call-set.bedpe")
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

