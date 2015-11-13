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
library(testthat)
source("libgridss.R")
source("libneochromosome.R")
source("libvcf.R")

theme_set(theme_bw())
power4 <- function(x) x^4
power4_trans <- function () {
  trans_new("power4", "power4", function(x) x^(1/4), domain = c(0, Inf))
}
scale_y_power4 <- function(...) {
  scale_y_continuous(..., trans = power4_trans())
}

minsize <- 51 # minimum Mills 2012 event size
# Strict matching
maxLengthRelativeDifference <- 0.25
maxPositionDifference_LengthOrEndpoint <- 25
# Permissive matching
libstddev <- 72
maxLengthRelativeDifference <- 0.5
maxPositionDifference_LengthOrEndpoint <- 2.5*libstddev

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
#  as(ucscgapblacklist, "RangedData"),
#  as(mappabilityblacklist, "RangedData"),
  as(encodeblacklist, "RangedData"),
  as(altassembly, "RangedData")), "GRanges")




truth <- "lumpyPacBioMoleculo"
vcfs <- lapply(vcfs, function(vcf) {
  attr(vcf, "metadata")$CX_REFERENCE_VCF <- "00000000000000000000000000000001.reference.vcf"
  return(vcf)
})
truth <- "merged"
vcfs <- lapply(vcfs, function(vcf) {
  attr(vcf, "metadata")$CX_REFERENCE_VCF <- "00000000000000000000000000000000.reference.vcf"
  return(vcf)
})
truth <- "Mills"
vcfs <- lapply(vcfs, function(vcf) {
  attr(vcf, "metadata")$CX_REFERENCE_VCF <- "00000000000000000000000000000002.reference.vcf"
  return(vcf)
})

# garvan <- readVcf("W:/na12878/garvan/garvan-na12878.vcf", "unknown")
# seqlevels(garvan) <- paste0("chr", seqlevels(garvan))
# attr(garvan, "metadata") <- attr(vcfs[["a136cfee2d40e62b4fd4be194366f291"]], "metadata")
# attr(garvan, "metadata")$CX_READ_LENGTH <- 150
# attr(garvan, "metadata")$CX_CALLER <- "gridss_garvan"
# attr(garvan, "metadata")$Id <- "11111111111111111111111111111111"
# attr(garvan, "metadata")$File <- "11111111111111111111111111111111"
# rownames(attr(garvan, "metadata")) <- "11111111111111111111111111111111"
# vcfs[["11111111111111111111111111111111"]] <- garvan
# Separate out GRIDSS confidence levels
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
# split out gridss subgraph
vcfs <- lapply(vcfs, function(vcf) {
  if (is.na(vcf@metadata$CX_CALLER_ARGS) || is.null(vcf@metadata$CX_CALLER_ARGS) || vcf@metadata$CX_CALLER_ARGS != "assembly.method") return (vcf)
  attr(vcf, "metadata")$CX_CALLER <- "gridss/0.9.0/subgraph"
  return(vcf)
})
# gridss evidence breakdown
gridssvcfs <- c(
  lapply(vcfs, function(vcf) {
    if (is.na(vcf@metadata$CX_CALLER_ARGS) || is.null(vcf@metadata$CX_CALLER_ARGS) || vcf@metadata$CX_CALLER_ARGS != "assembly.method") return (NULL)
    attr(vcf, "metadata")$CX_CALLER <- "Windowed de Bruijn graph assembly"
    df <- gridss.vcftodf(vcf)
    rowRanges(vcf)$FILTER <- "."
    rowRanges(vcf)$QUAL <- df$ASQ + df$RASQ
    #vcf <- vcf[df$RP + df$SR + df$RSR > 0,] # filter assembly-only variants due to high FP rate
    vcf <- vcf[df$ASCSR / df$ASSR < 0.8 & df$ASCRP / df$ASRP < 0.8,] # filter assembly variants with little actual support due to high FP rate
    return(vcf)
  }), unlist(recursive=FALSE, lapply(vcfs, function(vcf) {
    if (is.null(attr(vcf, "metadata")) || is.na(vcf@metadata$CX_CALLER) || str_split(vcf@metadata$CX_CALLER, "/")[[1]][1] !="gridss") return(NULL)
    if (!is.na(vcf@metadata$CX_CALLER_ARGS) && !is.null(vcf@metadata$CX_CALLER_ARGS) && vcf@metadata$CX_CALLER_ARGS == "assembly.method") return (NULL) # exclude subgraph assembly
    df <- gridss.vcftodf(vcf)
    return (lapply(list(
      list(f=function(df) df$QUAL, label="Combined"),
      list(f=function(df) df$ASQ + df$RASQ, label="Positional de Bruijn assembly"),
      list(f=function(df) df$RPQ + df$SRQ + df$RSRQ, label="RP SR"),
      list(f=function(df) df$RPQ, label="RP"),
      list(f=function(df) df$SRQ + df$RSRQ, label="SR"),
      list(f=function(df) df$ASRP + df$ASSR, label="Positional de Bruijn assembly$count"),
      list(f=function(df) df$RP + df$SR + df$RSR, label="RP SR$count"),
      list(f=function(df) df$RP, label="RP$count"),
      list(f=function(df) df$SR + df$RSR, label="SR$count")),
      function(scoring) {
        ovcf <- vcf
        rowRanges(ovcf)$QUAL <- scoring$f(df)
        rowRanges(ovcf)$FILTER <- "."
        attr(ovcf, "metadata")$CX_CALLER <- scoring$label
        return(ovcf)
      }))
  })))
gridssvcfs[sapply(gridssvcfs, is.null)] <- NULL

truthlist_filtered <- CalculateTruthSummary(vcfs, blacklist=blacklist, maxerrorbp=maxPositionDifference_LengthOrEndpoint, maxerrorpercent=maxLengthRelativeDifference, ignoreFilters=FALSE, ignore.strand=TRUE) # breakdancer does not specify strand
truthlist_all <- CalculateTruthSummary(vcfs, blacklist=blacklist, maxerrorbp=maxPositionDifference_LengthOrEndpoint, maxerrorpercent=maxLengthRelativeDifference, ignoreFilters=TRUE, ignore.strand=TRUE) # breakdancer does not specify strand

dtroc_filtered <- TruthSummaryToROC(truthlist_filtered, bylist=c("CX_CALLER"))
dtroc_filtered$Filter <- "Default"
dtroc_all <- TruthSummaryToROC(truthlist_all, bylist=c("CX_CALLER"))
dtroc_all$Filter <- "Including Filtered"
dtroc <- rbind(dtroc_all, dtroc_filtered)

ggplot(dtroc) + 
  aes(y=tp/2, x=fp/2, color=CX_CALLER, linetype=Filter) +
  geom_line() + 
#  scale_color_brewer(palette="Set2") + 
  coord_cartesian(xlim=c(0, 1000), ylim=c(0, 3000)) + 
  labs(y="tp", x="fp", title=paste("ROC curve NA12878 deletions", truth))
ggsave(paste0("na12878_tp_fp_", truth, "_error_", maxLengthRelativeDifference, "_", maxPositionDifference_LengthOrEndpoint, ".pdf"), width=7, height=5)
ggplot(dtroc) + 
  aes(y=prec, x=tp/2, color=CX_CALLER, linetype=Filter) +
#  scale_color_brewer(palette="Set2") + 
  geom_line() + 
  labs(y="Precision", x="True positives", title=paste("Precision-Recall curve NA12878 deletions", truth))
ggsave(paste0("na12878_prec_", truth, "_error_", maxLengthRelativeDifference, "_", maxPositionDifference_LengthOrEndpoint, ".pdf"), width=7, height=5)

dtroc[dtroc$prec > 0.75,][!duplicated(paste(dtroc[dtroc$prec > 0.75,]$CX_CALLER, dtroc[dtroc$prec > 0.75,]$Filter), fromLast=TRUE),]


#####################################
# Moleculo/PacBio truth
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
countBedHits <- function(gr, bed) {
  hits <- data.table(data.frame(findOverlaps(gr, bed, type="equal", maxgap=maxPositionDifference_LengthOrEndpoint, algorithm="intervaltree")))
  hits$reflength <- end(bed[hits$subjectHits]) - start(bed[hits$subjectHits])
  hits$calllength <- end(gr[hits$queryHits]) - start(gr[hits$queryHits])
  hits$percentsize <- abs(hits$calllength / hits$reflength)
  hits <- hits[abs(hits$reflength - hits$calllength) <= maxPositionDifference_LengthOrEndpoint,]
  if (!is.null(maxLengthRelativeDifference)) {
    hits <- hits[hits$percentsize >= 1 - maxLengthRelativeDifference & hits$percentsize <= 1 + maxLengthRelativeDifference,]
  }
  counts <- hits[, list(hitcount=.N), by=c("queryHits")]
  result <- rep(0, length(gr))
  result[counts$queryHits] <- counts$hitcount
  return(result)
}
longReadRoc <- function(vcflist) {
  # extract deletion calls from VCFs
  dfdelcalls <- rbindlist(lapply(vcflist, function(vcf) {
    caller <- attr(vcf, "metadata")$CX_CALLER
    if (is.na(caller) || is.null(caller)) return(NULL)
    gr <- vcftobpgr(vcf)
    gro <- gr[seq_along(gr) < gr$mateIndex,]
    grh <- gr[gro$mateIndex,]
    df <- data.frame(
      chrom1=seqnames(gro), start1=start(gro), end1=end(gro),
      chrom2=seqnames(grh), start2=start(grh), end2=end(grh),
      name=rep(caller, length(gro)), score=gro$QUAL,
      strand1=rep("+", length(gro)), strand2=rep("-", length(gro)),
      length=gro$size,
      pos1=gro$callPosition,
      pos2=grh$callPosition,
      id=gro$vcfid,
      filtered=!(gro$FILTER %in% c(".", "PASS")),
      row.names=NULL)
    return(df)
  }))
  #write.table(dfdelcalls, paste0(rootdir, "i/data.na12878/tovalidate.bedpe"), sep='\t', quote=FALSE, row.names=FALSE)
  delbed <- GRanges(seqnames=dfdelcalls$chrom2,ranges=IRanges(start=pmin(dfdelcalls$pos1, dfdelcalls$pos2), end=pmax(dfdelcalls$pos1, dfdelcalls$pos2)), caller=dfdelcalls$name, QUAL=dfdelcalls$score, length=dfdelcalls$length, filtered=dfdelcalls$filtered)
  delbed <- delbed[abs(delbed$length) >= minsize,]
  # match CalculateTruth blacklisting
  # any overlap
  #delbed$blacklisted <- overlapsAny(GRanges(seqnames=seqnames(delbed), ranges=IRanges(start=start(delbed), end=end(delbed))), blacklist, type="any", maxgap=maxPositionDifference_LengthOrEndpoint)
  # breakend overlap
  delbed$blacklisted <- overlapsAny(GRanges(seqnames=seqnames(delbed), ranges=IRanges(start=start(delbed), width=1)), blacklist, type="any", maxgap=maxPositionDifference_LengthOrEndpoint) |
    overlapsAny(GRanges(seqnames=seqnames(delbed), ranges=IRanges(end=end(delbed), width=1)), blacklist, type="any", maxgap=maxPositionDifference_LengthOrEndpoint)
  
  delbed$srmoleculo <- countBedHits(delbed, srmoleculo)
  delbed$srpacbio <- countBedHits(delbed, srpacbio)
  delbed$sr <- delbed$srmoleculo + delbed$srpacbio
  delbed$spmoleculo <- countBedHits(delbed, spmoleculo)
  delbed$sppacbio <- countBedHits(delbed, sppacbio)
  delbed$sp <- delbed$spmoleculo + delbed$sppacbio
  delbed$tp <- ifelse(delbed$sr >= 3 | delbed$sp >= 7, 1, 0)
  delbed$fp <- 1 - delbed$tp
  delbed <- delbed[order(-delbed$QUAL),]
  delroc <- data.table(as.data.frame(mcols(delbed)))
  delroc <- delroc[!delroc$blacklisted,]
  tmp <- delroc[!delroc$filtered & delroc$caller %in% unique(delroc[delroc$filtered,]$caller),]
  tmp$filtered <- TRUE
  delroc <- rbind(delroc, tmp)
  delroc <- delroc[order(-delroc$QUAL),]
  delroc[,`:=`(tp=cumsum(tp), fp=cumsum(fp), QUAL=cummin(QUAL)), by=c("caller", "filtered")]
  delroc <- delroc[!duplicated(delroc[, c("caller", "filtered", "QUAL"), with=FALSE], fromLast=TRUE),] # take only one data point per QUAL
  delroc$Filter <- ifelse(delroc$filtered, "Including Filtered", "Default")
  delroc$precision <- delroc$tp/(delroc$tp+delroc$fp)
  return (delroc)
}
delroc <- longReadRoc(vcfs)
ggplot(delroc) + aes(x=fp, y=tp, color=caller, linetype=Filter) + geom_line() +
  coord_cartesian(xlim=c(0, 1000), ylim=c(0, 3000))
  labs(x="False Positives", y="True Positives", title="NA12878 deletions PacBio/Moleculo validated")
#  scale_color_brewer(palette="Set2") +
ggsave(paste0("na12878_tp_fp_pacbiomoleculo", "_error_", maxLengthRelativeDifference, "_", maxPositionDifference_LengthOrEndpoint, ".pdf"), width=7, height=5)
ggplot(delroc) + aes(x=tp, y=precision, color=caller, linetype=Filter) + geom_line() +
  labs(x="True Positives", y="Precision", title="NA12878 deletions PacBio/Moleculo validated")
#  scale_color_brewer(palette="Set2") +
ggsave(paste0("na12878_prec_pacbiomoleculo", "_error_", maxLengthRelativeDifference, "_", maxPositionDifference_LengthOrEndpoint, ".pdf"), width=7, height=5)

gdelroc <- longReadRoc(gridssvcfs)
gdelroc$Scoring <- ifelse(str_detect(gdelroc$caller, "[$]"), "Read Count", "Bayesian")
gdelroc$caller <- str_extract(gdelroc$caller, "[^$]*")
gdelroc <- gdelroc[gdelroc$QUAL > 0,] # filter out calls that wouldn't have been called according to that particular scoring scheme
ggplot(gdelroc) + aes(x=fp, y=tp, color=caller, linetype=Scoring) + geom_line(size=3) +
  coord_cartesian(xlim=c(0, 1000), ylim=c(0, 3000)) + 
  scale_color_brewer(palette="Set2") +
  labs(x="False Positives", y="True Positives", title="gridss support breakdown")
ggsave(paste0("na12878_gridss_tp_fp_pacbiomoleculo", "_error_", maxLengthRelativeDifference, "_", maxPositionDifference_LengthOrEndpoint, ".pdf"), width=7, height=5)

ggplot(gdelroc) + aes(x=tp, y=precision, color=caller, linetype=Scoring) + geom_line(size=2) + 
  scale_color_brewer(palette="Set2") +
  labs(x="True Positives", y="Precision", title="gridss support breakdown")
ggsave(paste0("na12878_gridss_prec_pacbiomoleculo", "_error_", maxLengthRelativeDifference, "_", maxPositionDifference_LengthOrEndpoint, ".pdf"), width=7, height=5)


# precision thresholds
delroc[delroc$precision > 0.95,][!duplicated(paste(delroc[delroc$precision > 0.95,]$caller, delroc[delroc$precision > 0.95,]$Filter), fromLast=TRUE),]

ggplot(as.data.frame(delbed[delbed$caller=="gridss" & !delbed$filtered,])) + aes(x=QUAL, fill=ifelse(tp, "_tp", "fp")) + geom_histogram(binwidth=100) + scale_x_continuous(limits=c(0, 3000))
ggsave("na12878_gridss_hq_treshold.png")




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

