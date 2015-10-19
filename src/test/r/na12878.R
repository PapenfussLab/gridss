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
rootdir <- ifelse(as.character(Sys.info())[1] == "Windows", "W:/", "~/")

vcfs <- NULL
pwd <- getwd()
setwd(paste0(rootdir, "i/data.na12878"))
metadata <- LoadMetadata()
vcfs <- LoadVcfs(metadata, existingVcfs=vcfs)
setwd(pwd)

ucscgapblacklist <- import(paste0(rootdir, "projects/reference_genomes/human/blacklist_annotations/hg19_ucsc_gap_table.bed"))
mappabilityblacklist <- import(paste0(rootdir, "projects/reference_genomes/human/blacklist_annotations/wgEncodeDukeMapabilityRegionsExcludable.bed"))
encodeblacklist <- import(paste0(rootdir, "projects/reference_genomes/human/blacklist_annotations/wgEncodeDacMapabilityConsensusExcludable.bed"))
altassembly <- GRanges(seqnames=c("chr1_gl000191_random", "chr1_gl000191_random", "chr1_gl000192_random", "chr11_gl000202_random", "chr17_ctg5_hap1", "chr17_gl000203_random", "chr17_gl000204_random", "chr17_gl000205_random", "chr17_gl000206_random", "chr18_gl000207_random", "chr19_gl000208_random", "chr19_gl000209_random", "chr21_gl000210_random", "chr4_ctg9_hap1", "chr4_gl000193_random", "chr4_gl000194_random", "chr6_apd_hap1", "chr6_cox_hap2", "chr6_dbb_hap3", "chr6_mann_hap4", "chr6_mcf_hap5", "chr6_qbl_hap6", "chr6_ssto_hap7", "chr7_gl000195_random", "chr8_gl000196_random", "chr8_gl000197_random", "chr9_gl000198_random", "chr9_gl000199_random", "chr9_gl000200_random", "chr9_gl000201_random", "chrM", "chrUn_gl000211", "chrUn_gl000212", "chrUn_gl000213", "chrUn_gl000214", "chrUn_gl000215", "chrUn_gl000216", "chrUn_gl000217", "chrUn_gl000218", "chrUn_gl000219", "chrUn_gl000220", "chrUn_gl000221", "chrUn_gl000222", "chrUn_gl000223", "chrUn_gl000224", "chrUn_gl000225", "chrUn_gl000226", "chrUn_gl000227", "chrUn_gl000228", "chrUn_gl000229", "chrUn_gl000230", "chrUn_gl000231", "chrUn_gl000232", "chrUn_gl000233", "chrUn_gl000234", "chrUn_gl000235", "chrUn_gl000236", "chrUn_gl000237", "chrUn_gl000238", "chrUn_gl000239", "chrUn_gl000240", "chrUn_gl000241", "chrUn_gl000242", "chrUn_gl000243", "chrUn_gl000244", "chrUn_gl000245", "chrUn_gl000246", "chrUn_gl000247", "chrUn_gl000248", "chrUn_gl000249", "GL000191.1", "GL000192.1", "GL000193.1", "GL000194.1", "GL000195.1", "GL000196.1", "GL000197.1", "GL000198.1", "GL000199.1", "GL000200.1", "GL000201.1", "GL000202.1", "GL000203.1", "GL000204.1", "GL000205.1", "GL000206.1", "GL000207.1", "GL000208.1", "GL000209.1", "GL000210.1", "GL000211.1", "GL000212.1", "GL000213.1", "GL000214.1", "GL000215.1", "GL000216.1", "GL000217.1", "GL000218.1", "GL000219.1", "GL000220.1", "GL000221.1", "GL000222.1", "GL000223.1", "GL000224.1", "GL000225.1", "GL000226.1", "GL000227.1", "GL000228.1", "GL000229.1", "GL000230.1", "GL000231.1", "GL000232.1", "GL000233.1", "GL000234.1", "GL000235.1", "GL000236.1", "GL000237.1", "GL000238.1", "GL000239.1", "GL000240.1", "GL000241.1", "GL000242.1", "GL000243.1", "GL000244.1", "GL000245.1", "GL000246.1", "GL000247.1", "GL000248.1", "GL000249.1", "hs37d5", "MT", "NC_007605", "chrGL000191.1", "chrGL000192.1", "chrGL000193.1", "chrGL000194.1", "chrGL000195.1", "chrGL000196.1", "chrGL000197.1", "chrGL000198.1", "chrGL000199.1", "chrGL000200.1", "chrGL000201.1", "chrGL000202.1", "chrGL000203.1", "chrGL000204.1", "chrGL000205.1", "chrGL000206.1", "chrGL000207.1", "chrGL000208.1", "chrGL000209.1", "chrGL000210.1", "chrGL000211.1", "chrGL000212.1", "chrGL000213.1", "chrGL000214.1", "chrGL000215.1", "chrGL000216.1", "chrGL000217.1", "chrGL000218.1", "chrGL000219.1", "chrGL000220.1", "chrGL000221.1", "chrGL000222.1", "chrGL000223.1", "chrGL000224.1", "chrGL000225.1", "chrGL000226.1", "chrGL000227.1", "chrGL000228.1", "chrGL000229.1", "chrGL000230.1", "chrGL000231.1", "chrGL000232.1", "chrGL000233.1", "chrGL000234.1", "chrGL000235.1", "chrGL000236.1", "chrGL000237.1", "chrGL000238.1", "chrGL000239.1", "chrGL000240.1", "chrGL000241.1", "chrGL000242.1", "chrGL000243.1", "chrGL000244.1", "chrGL000245.1", "chrGL000246.1", "chrGL000247.1", "chrGL000248.1", "chrGL000249.1", "chrhs37d5", "chrMT"),
                       ranges=IRanges(start=0, end=1000000000))
altassembly$name <- as.character(seqnames(altassembly))
altassembly$score <- 1000
ucscgapblacklist$score <- 1000
blacklist <- as(rbind(
  as(ucscgapblacklist, "RangedData"),
  as(mappabilityblacklist, "RangedData"),
  as(encodeblacklist, "RangedData"),
  as(altassembly, "RangedData")), "GRanges")




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
  caller <- str_extract(attr(vcf, "metadata")$CX_CALLER, "^[^/]+")
  if (!is.na(caller) && !is.null(caller) && caller %in% c("breakdancer")) {
    # strip all BND events since they are non-deletion events such as translocations
    vcf <- vcf[info(vcf)$SVTYPE == "DEL",]
  }
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
  scale_color_brewer(palette="Set2") + 
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



#####################################
# extract deletion calls from VCFs
dfdelcalls <- rbindlist(lapply(names(vcfs)[!(names(vcfs) %in% c("00000000000000000000000000000000", "00000000000000000000000000000001", "00000000000000000000000000000002"))], function(id) {
  gr <- vcftobpgr(vcfs[[id]])
  gro <- gr[seq_along(gr) < gr$mateIndex,]
  grh <- gr[gro$mateIndex,]
  df <- data.frame(
    chrom1=seqnames(gro), start1=start(gro), end1=end(gro),
    chrom2=seqnames(grh), start2=start(grh), end2=end(grh),
    name=rep(attr(vcfs[[id]], "metadata")$CX_CALLER, length(gro)), score=gro$QUAL,
    strand1=rep("+", length(gro)), strand2=rep("-", length(gro)),
    length=gro$size,
    pos1=gro$callPosition,
    pos2=grh$callPosition,
    id=gro$vcfid,
    filtered=!(gr$FILTER %in% c(".", "PASS")),
    row.names=NULL)
  return(df)
}))
write.table(dfdelcalls, paste0(rootdir, "i/data.na12878/tovalidate.bedpe"), sep='\t', quote=FALSE, row.names=FALSE)

delbed <- GRanges(seqnames=dfdelcalls$chrom2,ranges=IRanges(start=pmin(dfdelcalls$pos1, dfdelcalls$pos2), end=pmax(dfdelcalls$pos1, dfdelcalls$pos2)), caller=dfdelcalls$name, QUAL=dfdelcalls$score, length=dfdelcalls$length, filtered=dfdelcalls$filtered)
delbed <- delbed[abs(delbed$length) >= minsize,]
# should we blacklist if either breakend is within the interval? Should we add error margin to blacklisted regions?
delbed$blacklisted <- overlapsAny(delbed, blacklist, type="within")

srpacbio <- c(
  withChr(import.bed(con=paste0(rootdir, "na12878/longread/chemistry1.sorted.bam.sr.bam.sr.bed"))),
  withChr(import.bed(con=paste0(rootdir, "na12878/longread/chemistry_2_picard.bam.sr.bam.sr.bed"))),
  withChr(import.bed(con=paste0(rootdir, "na12878/longread/chemistry_3_picard.bam.sr.bam.sr.bed"))),
  withChr(import.bed(con=paste0(rootdir, "na12878/longread/NA12878.pacbio_fr_MountSinai.bwa-sw.20140211.bam.sr.bam.sr.bed"))))
srmoleculo <- withChr(import.bed(con=paste0(rootdir, "na12878/longread/NA12878.moleculo.bwa-mem.20140110.bam.sr.bam.sr.bed")))
sppacbio <- c(
  withChr(import.bed(con=paste0(rootdir, "na12878/longread/chemistry1.sorted.bam.sp.bed"))),
  withChr(import.bed(con=paste0(rootdir, "na12878/longread/chemistry_2_picard.bam.sp.bed"))),
  withChr(import.bed(con=paste0(rootdir, "na12878/longread/chemistry_3_picard.bam.sp.bed"))),
  withChr(import.bed(con=paste0(rootdir, "na12878/longread/NA12878.pacbio_fr_MountSinai.bwa-sw.20140211.bam.sp.bed"))))
spmoleculo <- withChr(import.bed(con=paste0(rootdir, "na12878/longread/NA12878.moleculo.bwa-mem.20140110.bam.sp.bed")))
# Split read validation
delbed$srmoleculo <- countOverlaps(delbed, delbedmoleculo, type="equal", maxgap=2.5*72, algorithm="intervaltree")
delbed$srpacbio <- countOverlaps(delbed, delbedpacbio, type="equal", maxgap=2.5*72, algorithm="intervaltree")
delbed$sr <- delbed$srmoleculo + delbed$srpacbio

# filter blacklisted

delroc <- data.table(as.data.frame(mcols(delbed[order(-delbed$QUAL),])))
delroc <- delroc[!delroc$blacklisted,]
delroc$tp <- ifelse(delroc$sr >= 3, 1, 0)
delroc$fp <- 1 - delroc$tp
delroc[,`:=`(tp=cumsum(tp), fp=cumsum(fp), QUAL=cummin(QUAL)), by=c("caller")]
delroc <- delroc[!duplicated(delroc[, c("caller", "QUAL"), with=FALSE], fromLast=TRUE),] # take only one data point per QUAL
ggplot(delroc) + aes(x=fp, y=tp, color=caller) + geom_line() + scale_x_continuous(limits=c(0, 1000))
ggplot(delroc) + aes(x=tp, y=tp/(tp+fp), color=caller) + geom_line() + geom_point()


# run na12878_validate.sh to extract moleculo & pacbio support
dtann <- read.delim(paste0(rootdir, "i/data.na12878/tovalidate.annotated.bedpe"), stringsAsFactors=FALSE)
minSpanningSizeMultiple <- 0.5
maxSpanningSizeMultiple <- 1.5
maxSpanningLengthDifference
unpack <- function(charvector) {
  charlist <- str_split(charvector, "[\\[\\],]")
  df <- data.frame(ordinal=rep(seq_along(charlist), times=elementLengths(charlist)), value=as.integer(unlist(charlist)))
  df <- df[!is.na(df$value),]
  return(df)
}
countSpanningDeletions <- function(dtann, spans, columnName) {
  spans$calllength <- abs(dtann$length[spans$ordinal])
  spans$name <- dtann$name[spans$ordinal]
  spans$spanlength <- spans$value
  spans$offby <- spans$spanlength / spans$calllength
  # lax matching requirements
  hitspans <- spans[spans$offby >= minSpanningSizeMultiple & spans$offby <= maxSpanningSizeMultiple & abs(spans$calllength - spans$spanlength) <= maxSpanningLengthDifference,]
  hitspans <- data.table(hitspans)[, list(spanCount=.N), by=c("ordinal", "name")] # TODO: split out pacbio & moleculo"src", 
  hitspans[, columnName] <- hitspans$spanCount
  hitspans$spanCount <- NULL
  dtann$ordinal <- seq_len(nrow(dtann))
  dtann <- merge(data.table(dtann), hitspans, by=c("ordinal", "name"), all.x=TRUE)
  dtann$ordinal <- NULL
  return(dtann)
}
# spans <- rbind(
#   unpack(dtann$SpanningDeletionSize, "moleculo"),
#   unpack(dtann$SpanningDeletionSize.1, "pacbio chemistry1"),
#   unpack(dtann$SpanningDeletionSize.2, "pacbio chemistry2"),
#   unpack(dtann$SpanningDeletionSize.3, "pacbio chemistry3"),
#   unpack(dtann$SpanningDeletionSize.4, "pacbio MountSinai"))
# spans$calllength <- abs(dtann$length[spans$ordinal])
# spans$caller <- dtann$name[spans$ordinal]
# spans$spanlength <- spans$value
# spans$offby <- spans$spanlength / spans$calllength
# spans2 <- spans[spans$offby > 0.4,]
# ggplot(spans[sample.int(nrow(spans), size=10000),]) + aes(y=calllength, x=spanlength, color=src, size=0.1) + geom_point() + scale_x_log10(limits=c(10, 10000)) + scale_y_log10(limits=c(10, 10000))
# ggplot(spans[sample.int(nrow(spans), size=10000),]) + aes(y=spanlength, x=calllength, color=src, shape=caller) + geom_point()
# ggplot(spans2[sample.int(nrow(spans2), size=100000),]) + aes(y=spanlength, x=offby, color=caller, alpha=0.1) + geom_point(size=0.1) + scale_x_continuous(limits=c(0, 1.7)) + scale_y_log10() + facet_wrap(~ src)
# ggplot(dtann) + aes(color=name, x=-length) + geom_histogram() + scale_x_log10()
#hitspans <- spans[spans$offby >= 0.5 & spans$offby <= 1.5 & abs(spans$calllength - spans$spanlength) <= 212,]
#hitspans <- data.table(hitspans)[, list(spanCount=.N), by=c("ordinal", "caller")] # TODO: split out pacbio & moleculo"src", 
#ggplot(hitspans) + aes(x=spanCount) + geom_histogram(binwidth=1)

dtann <- countSpanningDeletions(dtann, unpack(dtann$SpanningDeletionSize), "spanning_moleculo")
dtann <- countSpanningDeletions(dtann, c(
  unpack(dtann$SpanningDeletionSize.1), # chem1
  unpack(dtann$SpanningDeletionSize.2), # chem2
  unpack(dtann$SpanningDeletionSize.3), # chem3
  unpack(dtann$SpanningDeletionSize.4)), # Mount Sinai
  "spanning_pacbio")













##########################
## GRIDSS single sample
# 
# millsgr <- bedpe2grmate(paste0(rootdir, "na12878/lumpy-Mills2012-call-set.bedpe"))
# pacbiogr <- bedpe2grmate(paste0(rootdir, "na12878/lumpy-PacBioMoleculo-call-set.bedpe"))
# mcols(millsgr) <- NULL
# mcols(pacbiogr) <- NULL
# write.csv(names(millsgr[overlapsAny(millsgr, pacbiogr),]), "W:/na12878/millsoverlappingpacbio.csv")
# 
# millsgr <- bedpe2grmate(paste0(rootdir, "na12878/lumpy-Mills2012-call-set.bedpe"))
# pacbiogr <- bedpe2grmate(paste0(rootdir, "na12878/lumpy-PacBioMoleculo-call-set.bedpe"))
# vcf <- readVcf("~/na12878/platinum-na12878.vcf", "hg19")
# gvcf <- readVcf("~/garvan-na12878.vcf")
# 
# vcf <- vcf[rowRanges(vcf)$QUAL >= 100,]
# vcf <- vcf[rowRanges(vcf)$FILTER %in% c(".", "PASS"),]
# vcf <- vcf[gridss.vcftodf(vcf)$assembly %in% c("Both"), ] # high qual
# vcf <- vcf[(str_detect(as.character(alt(vcf)), "\\[") & str_detect(names(rowRanges(vcf)), "o$")) |
#              (str_detect(as.character(alt(vcf)), "]") & str_detect(names(rowRanges(vcf)), "h$")),]
# # more bases deleted than added
# vcf <- gridss.removeUnpartnerededBreakend(vcf)
# vcf <- vcf[abs(start(rowRanges(vcf)) - start(rowRanges(vcf)[as.character(info(vcf)$MATEID),])) > nchar(gridss.vcftodf(vcf)$INSSEQ),]
# vcf <- gridss.removeUnpartnerededBreakend(vcf)
# vcf <- vcf[seqnames(rowRanges(vcf)) == seqnames(rowRanges(vcf)[as.character(info(vcf)$MATEID),]),] # intra-chromsomal events
# vcf <- gridss.removeUnpartnerededBreakend(vcf)
# 
# #seqlevels(millsgr) <- sub("chr", "", seqlevels(millsgr))
# mills <- gridss.annotateBreakpoints(millsgr, millsgr[millsgr$mate,], vcf, maxgap=100, ignore.strand=TRUE)
# mills$bed$truth <- "Mills"
# mills$gridss$truth <- "Mills"
# 
# #seqlevels(pacbiogr) <- sub("chr", "", seqlevels(pacbiogr))
# pb <- gridss.annotateBreakpoints(pacbiogr, pacbiogr[pacbiogr$mate,], vcf, maxgap=100)
# pb$bed$truth <- "PacBioOverlap"
# pb$gridss$truth <- "PacBioOverlap"
# 
# mills$bed$V11 <- NULL
# out <- list(bed=c(mills$bed, pb$bed), gridss=rbind(mills$gridss, pb$gridss))
# out$bed$size <- abs(start(out$bed) - start(out$bed[out$bed$mate,]))
# 
# out$gridss <- out$gridss[order(-out$gridss$QUAL),]
# roc <- data.table(out$gridss)[, list(cumtp=cumsum(!is.na(bedid)), cumfp=cumsum(is.na(bedid))), by=list(truth)]
# ggplot(roc) + aes(x=cumfp/2, y=cumtp/2, color=truth) + geom_point() + geom_line()
# ggsave("na12878_lumpy_fig3_roc.png")
# 
# sum(!is.na(out$bed$gridssid)) / length(out$bed) # sens
# sum(!is.na(out$gridss$bedid)) / nrow(out$gridss) # prec
# ggplot(as.data.frame(out$bed)) + aes(x=size, fill=ifelse(is.na(gridssid), "FN", "TP")) + geom_histogram() + scale_x_log10() + facet_grid(truth ~ . ) +
#   labs(fill="", title="Gridss NA12878 Sensitivity")
# ggsave("na12878_sens.png")
# ggplot(as.data.frame(out$gridss)) + aes(x=size, fill=ifelse(is.na(bedid), "FP", "TP")) + geom_histogram() + scale_x_log10() + facet_grid(truth ~ . ) +
#   labs(fill="", title="Gridss NA12878 Precision by event size")
# ggsave("na12878_prec_size.png")
# ggplot(as.data.frame(out$gridss)) + aes(x=QUAL, fill=!is.na(bedid)) + geom_histogram() + scale_x_log10() + facet_grid(truth ~ . ) +
#   labs(fill="", title="Gridss NA12878 Precision by called qual")
# ggsave("na12878_prec_qual.png")
# 
# ggplot(aes(data=as.data.frame(out$bed), x=size, fill=!is.na(bedid))) + geom_histogram() + scale_x_log10()


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

