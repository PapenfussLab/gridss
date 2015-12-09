##################
# Processing steps (approximate runtime: 186h wall , 243h CPU)
##################
# - Follow na12878.R processing steps
# - Create symbolic links from ~/i/data.na12878/00000000000000000000000000000000.* to ~/i/data.assembly/00000000000000000000000000000000.*
# - ./call_gridss.sh assembly
# - Run this script

library(data.table)
library(ggplot2)
library(VariantAnnotation)
library(GenomicAlignments)
library(Rsamtools)
source("libvcf.R")
rootdir <- ifelse(as.character(Sys.info())[1] == "Windows", "W:/", "~/")
pwd <- getwd()
setwd(paste0(rootdir, "i/data.assembly"))
kmetadata <- LoadMetadata()
setwd(pwd)

theme_set(theme_bw())

kmetadata$breakendassembly <- paste0(rootdir, "i/data.assembly/", kmetadata$Id, "/breakend.vcf.idsv.working/breakend.vcf.idsv.breakend.bam")
kmetadata$assembly <- paste0(rootdir, "i/data.assembly/", kmetadata$Id, "/breakend.vcf.idsv.working/breakend.vcf.idsv.assembly.bam")

# load repeatmasker annotations
rmdt <- rbindlist(lapply(
  list.files(paste0(rootdir, "Papenfuss_lab/projects/reference_genomes/human/hg19/UCSC/repeatmasker/"), pattern="*.fa.out", recursive=TRUE, full.names=TRUE),
  function (file) {
    read.table(file=file, sep="", quote="", skip=3)
}))
grrm <- GRanges(
  seqnames=rmdt$V5,
  ranges=IRanges(start=rmdt$V6 + 1, end=rmdt$V7),
  repeatType=rmdt$V10,
  repeatClass=rmdt$V11)
grrm$repeatClass <- str_replace(str_replace(grrm$repeatClass, "[?]", ""), "/.*", "")

# Finds the subject interval with the largest overlap with the query sequence
findBestOverlap <- function(query, subject, ...) {
  hits <- findOverlaps(query, subject, ...)
  hits <- data.table(queryHits=queryHits(hits), subjectHits=subjectHits(hits))
  # only take the best (longest) hit for each assembly
  hits$width <- pmin(end(query)[hits$queryHits], end(subject)[hits$subjectHits]) - pmax(start(query)[hits$queryHits], start(subject)[hits$subjectHits])
  hits <- hits[order(-hits$width),]
  hits <- hits[!duplicated(hits$queryHits),]
  return(hits)
}
annotateRepeat <- function(gr, grrm) {
  gr$repeatType <- "background"
  gr$repeatClass <- "background"
  hits <- findBestOverlap(gr, grrm, ignore.strand=TRUE, type="any", select="all")
  gr$repeatType[hits$queryHits] <- as.character(grrm$repeatType[hits$subjectHits])
  gr$repeatType <- as.factor(gr$repeatType)
  gr$repeatClass[hits$queryHits] <- grrm$repeatClass[hits$subjectHits]
  return(gr)
}
getassembly <- function(file) {
  print(paste("Processing", file))
  bam <- BamFile(file=file)
  bamrecords <- scanBam(bam, param=ScanBamParam(flag=scanBamFlag(isFirstMateRead=TRUE), what=c("rname", "pos", "qwidth", "flag", "cigar"), tag=c("qp", "qs", "sp", "dp", "sc")))[[1]]
  grbam <- GRanges(seqnames=bamrecords$rname, ranges=IRanges(start=bamrecords$pos, width=bamrecords$qwidth))
  grbam <- annotateRepeat(grbam, grrm)
  dt <- data.table(
    score=unlist(bamrecords$tag$qp) + unlist(bamrecords$tag$qs),
    realigned=!bamFlagTest(bamrecords$flag, "hasUnmappedMate"),
    spanning=!is.na(bamrecords$tag$sp) & bamrecords$tag$sp == "y",
    unanchored=str_detect(bamrecords$cigar, "X"),
    readCount=unlist(bamrecords$tag$sc) + unlist(bamrecords$tag$dp),
    repeatClass=grbam$repeatClass,
    repeatType=grbam$repeatType)
  return(dt)
}
getbreakendassembly <- function(file) {
  print(paste("Processing", file))
  bam <- BamFile(file=file)
  bamrecords <- scanBam(bam, param=ScanBamParam(what=c("rname", "pos", "qwidth", "cigar"), tag=c("oc")))[[1]]
  grbam <- GRanges(seqnames=bamrecords$rname, ranges=IRanges(start=bamrecords$pos, width=bamrecords$qwidth))
  grbam <- annotateRepeat(grbam, grrm)
  dt <- data.table(
    spanning=str_detect(bamrecords$cigar, "^[0-9]+M.*M$"),
    spanningPreRealignment=str_detect(ifelse(is.na(bamrecords$tag$oc), bamrecords$cigar, bamrecords$tag$oc),  "^[0-9]+M.*M$"),
    repeatClass=grbam$repeatClass,
    repeatType=grbam$repeatType)
  return(dt)
}
refbam <- BamFile(file=kmetadata[!is.na(kmetadata$CX_CALLER) & file.exists(kmetadata$assembly),]$assembly[1])
grgenome <- GRanges(seqnames=seqnames(seqinfo(refbam)), ranges=IRanges(start=1, width=seqlengths(seqinfo(refbam))))
grgenome$repeatType <- "genome"
grgenome$repeatClass <- "genome"
gr <- c(grgenome, grrm)
dtrc <- data.table(repeatClass=gr$repeatClass, length=width(gr))[, list(repeatClassLength=sum(length)), by="repeatClass"]
dtrc <- rbind(dtrc, data.table(repeatClass="background", repeatClassLength=2 * dtrc$repeatClassLength[dtrc$repeatClass=="genome"] - sum(dtrc$repeatClassLength)))


dtass <- data.table(kmetadata[!is.na(kmetadata$CX_CALLER) & file.exists(kmetadata$assembly),])[, getassembly(.SD$assembly), by=c("CX_ASSEMBLY_METHOD", "CX_K")]
dtass <- merge(dtass, dtrc, by="repeatClass")
dtass$Method <- ifelse(dtass$CX_ASSEMBLY_METHOD == "Subgraph", "Windowed", "Positional")
dtass$kmer <- dtass$CX_K
dtass$type <- ifelse(dtass$spanning, "Flanking anchors", ifelse(dtass$unanchored, "Unanchored", "Single ancor"))

dtbe <- data.table(kmetadata[!is.na(kmetadata$CX_CALLER) & file.exists(kmetadata$breakendassembly),])[, getbreakendassembly(.SD$breakendassembly), by=c("CX_ASSEMBLY_METHOD", "CX_K")]
dtbe <- merge(dtbe, dtrc, by="repeatClass")
dtbe$Method <- ifelse(dtbe$CX_ASSEMBLY_METHOD == "Subgraph", "Windowed", "Positional")
dtbe$kmer <- dtbe$CX_K
dtbe$type <- ifelse(dtbe$spanning, "Flanking anchors", "Single-sided")
dtbe$typePreRealignment <- ifelse(dtbe$spanningPreRealignment, "Flanking anchors", "Single-sided")


ggplot(dtass[,list(rate=.N/max(repeatClassLength)), by=c("repeatClass", "Method", "kmer")]) +
  aes(x=as.numeric(as.character(kmer)), y=rate*1000000, color=Method) +
  geom_line() +
  geom_point() +
  facet_wrap(~ repeatClass) +
  scale_y_log10() + 
  labs(title="Assembly rate by repeat class", x="kmer", y="Assemblies per Megabase")
ggsave("na12878_gridss_assembly_rate_repeat.png", width=10, height=7.5)

dtrep <- dtass[,list(rate=.N/max(repeatClassLength)), by=c("repeatClass", "Method", "kmer")]
dtrepBase <- dtrep[dtrep$repeatClass=="background",]
dtrepBase$baserate <- dtrepBase$rate
dtrepBase$rate <- NULL
dtrepBase$repeatClass <- NULL
dtrep <- merge(dtrep, dtrepBase, by=c("Method", "kmer"))
dtrep$enrichment <- dtrep$rate / dtrep$baserate
ggplot(dtrep[dtrep$repeatClass!="background",]) +
  aes(x=as.numeric(as.character(kmer)), y=enrichment, color=Method) +
  geom_line() +
  geom_point() +
  facet_wrap(~ repeatClass) +
  labs(title="Increase in assembly rate in repetative regions", x="kmer", y="Incease in assembly rate")
ggsave("na12878_gridss_assembly_rate_repeat_enrichment.png", width=10, height=7.5)

ggplot(dtass[,list(rate=.N/max(repeatClassLength)), by=c("type", "repeatClass", "Method", "kmer")]) +
  aes(x=as.numeric(as.character(kmer)), y=rate*1000000, color=Method, linetype=type) +
  geom_line() +
  geom_point() +
  facet_wrap(~ repeatClass) +
  scale_y_log10() + 
  labs(title="Assembly rate by repeat class", x="kmer", y="Assemblies per Megabase")

ggplot(dtass[,list(rate=.N/dtrc$repeatClassLength[dtrc$repeatClass=="genome"]), by=c("type", "Method", "kmer")]) +
  aes(x=as.numeric(as.character(kmer)), y=rate*1000000, color=Method, linetype=type) +
  geom_line() +
  geom_point() +
  scale_y_log10() + 
  labs(title="Assembly rate", x="kmer", y="Assemblies per Megabase")
ggsave("na12878_gridss_assembly_rate.png", width=10, height=7.5)

ggplot(dtass, aes(x=score, color=kmer, linetype=Method)) +
  geom_density() +
  facet_wrap(~ repeatClass) + 
  scale_color_brewer(type="div", palette="RdYlBu") + 
  scale_x_log10() + 
  labs(title="Assembly quality score distribution by repeat class", x="Gridss Assembly Quality Score")
ggsave("na12878_gridss_assembly_score_repeat.png", width=10, height=7.5)

ggplot(dtass) + aes(x=score, color=kmer) +
  geom_density() +
  facet_grid(type ~ Method) + 
  scale_color_brewer(type="div", palette="RdYlBu") + 
  coord_cartesian(ylim=c(0, 10)) + 
  scale_x_log10(limits=c(50,2000)) + 
  labs(title="Assembly quality score distribution by assembly type", x="Gridss Assembly Quality Score")
ggsave("na12878_gridss_assembly_score_type.png", width=10, height=7.5)

ggplot(dtass[sample(nrow(dtass), 1000000)]) +
  aes(x=score, y=readCount, color=kmer) +
  geom_density2d() +
  facet_grid(Method ~ type) + 
  scale_x_log10(limits=c(50,2000)) + scale_y_log10() + 
  scale_color_brewer(type="div", palette="RdYlBu") + 
  labs(title="Assembly rate", x="Gridss Assembly Quality Score", y="Assembly Read Count")
ggsave("na12878_gridss_assembly_qual_v_read_count_by_type.png", width=10, height=7.5)

for (rc in unique(dtass$repeatClass)) {
  ggplot(dtass[dtass$repeatClass==rc,]) +
    aes(x=score, y=readCount, color=kmer) +
    geom_density2d() +
    facet_grid(repeatClass + Method ~ type) + 
    scale_x_log10(limits=c(50,2000)) + scale_y_log10() + 
    scale_color_brewer(type="div", palette="RdYlBu") + 
    labs(title=paste(rc, "Assembly rate"), x="Gridss Assembly Quality Score", y="Assembly Read Count")
  ggsave(paste0("na12878_gridss_assembly_qual_v_read_count_by_type_rc", rc, ".png"), width=10, height=7.5)
}


prerealignment <- dtbe[,list(rate=.N/max(repeatClassLength), stage="Post realignment"), by=c("type", "Method", "kmer")]
postrealignment <- dtbe[,list(rate=.N/max(repeatClassLength), stage="Pre realignment"), by=c("typePreRealignment", "Method", "kmer")]
postrealignment$type <- postrealignment$typePreRealignment
postrealignment$typePreRealignment <- NULL
ggplot(rbind(prerealignment,postrealignment)) + 
  aes(x=as.numeric(as.character(kmer)), y=rate*1000000, color=Method, linetype=type) +
  geom_line() +
  geom_point() +
  facet_grid(~ stage) +
  scale_y_log10() + 
  labs(title="Assembly rate by anchor type", x="kmer", y="Assemblies per Megabase")
ggsave("na12878_gridss_assembly_rate_anchor_type.png", width=10, height=4)






