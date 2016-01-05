library(ggplot2)
library(rtracklayer)
library(data.table)
library(stringr)
library(scales)
library(parallel)
library(foreach)
library(doParallel)
library(testthat)
source("common.R")


#####################################
# NA12878 matching criteria

minsize <- 51 # minimum Mills 2012 event size

# Strict matching
maxLengthRelativeDifference <- 0.25
maxPositionDifference_LengthOrEndpoint <- 25
# Permissive matching
libstddev <- 72
maxLengthRelativeDifference <- 0.5
maxPositionDifference_LengthOrEndpoint <- 2.5*libstddev

#####################################
# hg19 blacklist
if (!exists("blacklist")) {
  #ucscgapblacklist <- import(paste0(rootdir, "projects/reference_genomes/human/blacklist_annotations/hg19_ucsc_gap_table.bed"))
  #mappabilityblacklist <- import(paste0(rootdir, "projects/reference_genomes/human/blacklist_annotations/wgEncodeDukeMapabilityRegionsExcludable.bed"))
  encodeblacklist <- import(paste0(rootdir, "projects/reference_genomes/human/blacklist_annotations/wgEncodeDacMapabilityConsensusExcludable.bed"))
  altassembly <- GRanges(seqnames=c("chr1_gl000191_random", "chr1_gl000191_random", "chr1_gl000192_random", "chr11_gl000202_random", "chr17_ctg5_hap1", "chr17_gl000203_random", "chr17_gl000204_random", "chr17_gl000205_random", "chr17_gl000206_random", "chr18_gl000207_random", "chr19_gl000208_random", "chr19_gl000209_random", "chr21_gl000210_random", "chr4_ctg9_hap1", "chr4_gl000193_random", "chr4_gl000194_random", "chr6_apd_hap1", "chr6_cox_hap2", "chr6_dbb_hap3", "chr6_mann_hap4", "chr6_mcf_hap5", "chr6_qbl_hap6", "chr6_ssto_hap7", "chr7_gl000195_random", "chr8_gl000196_random", "chr8_gl000197_random", "chr9_gl000198_random", "chr9_gl000199_random", "chr9_gl000200_random", "chr9_gl000201_random", "chrM", "chrUn_gl000211", "chrUn_gl000212", "chrUn_gl000213", "chrUn_gl000214", "chrUn_gl000215", "chrUn_gl000216", "chrUn_gl000217", "chrUn_gl000218", "chrUn_gl000219", "chrUn_gl000220", "chrUn_gl000221", "chrUn_gl000222", "chrUn_gl000223", "chrUn_gl000224", "chrUn_gl000225", "chrUn_gl000226", "chrUn_gl000227", "chrUn_gl000228", "chrUn_gl000229", "chrUn_gl000230", "chrUn_gl000231", "chrUn_gl000232", "chrUn_gl000233", "chrUn_gl000234", "chrUn_gl000235", "chrUn_gl000236", "chrUn_gl000237", "chrUn_gl000238", "chrUn_gl000239", "chrUn_gl000240", "chrUn_gl000241", "chrUn_gl000242", "chrUn_gl000243", "chrUn_gl000244", "chrUn_gl000245", "chrUn_gl000246", "chrUn_gl000247", "chrUn_gl000248", "chrUn_gl000249", "GL000191.1", "GL000192.1", "GL000193.1", "GL000194.1", "GL000195.1", "GL000196.1", "GL000197.1", "GL000198.1", "GL000199.1", "GL000200.1", "GL000201.1", "GL000202.1", "GL000203.1", "GL000204.1", "GL000205.1", "GL000206.1", "GL000207.1", "GL000208.1", "GL000209.1", "GL000210.1", "GL000211.1", "GL000212.1", "GL000213.1", "GL000214.1", "GL000215.1", "GL000216.1", "GL000217.1", "GL000218.1", "GL000219.1", "GL000220.1", "GL000221.1", "GL000222.1", "GL000223.1", "GL000224.1", "GL000225.1", "GL000226.1", "GL000227.1", "GL000228.1", "GL000229.1", "GL000230.1", "GL000231.1", "GL000232.1", "GL000233.1", "GL000234.1", "GL000235.1", "GL000236.1", "GL000237.1", "GL000238.1", "GL000239.1", "GL000240.1", "GL000241.1", "GL000242.1", "GL000243.1", "GL000244.1", "GL000245.1", "GL000246.1", "GL000247.1", "GL000248.1", "GL000249.1", "hs37d5", "MT", "NC_007605", "chrGL000191.1", "chrGL000192.1", "chrGL000193.1", "chrGL000194.1", "chrGL000195.1", "chrGL000196.1", "chrGL000197.1", "chrGL000198.1", "chrGL000199.1", "chrGL000200.1", "chrGL000201.1", "chrGL000202.1", "chrGL000203.1", "chrGL000204.1", "chrGL000205.1", "chrGL000206.1", "chrGL000207.1", "chrGL000208.1", "chrGL000209.1", "chrGL000210.1", "chrGL000211.1", "chrGL000212.1", "chrGL000213.1", "chrGL000214.1", "chrGL000215.1", "chrGL000216.1", "chrGL000217.1", "chrGL000218.1", "chrGL000219.1", "chrGL000220.1", "chrGL000221.1", "chrGL000222.1", "chrGL000223.1", "chrGL000224.1", "chrGL000225.1", "chrGL000226.1", "chrGL000227.1", "chrGL000228.1", "chrGL000229.1", "chrGL000230.1", "chrGL000231.1", "chrGL000232.1", "chrGL000233.1", "chrGL000234.1", "chrGL000235.1", "chrGL000236.1", "chrGL000237.1", "chrGL000238.1", "chrGL000239.1", "chrGL000240.1", "chrGL000241.1", "chrGL000242.1", "chrGL000243.1", "chrGL000244.1", "chrGL000245.1", "chrGL000246.1", "chrGL000247.1", "chrGL000248.1", "chrGL000249.1", "chrhs37d5", "chrMT"),
                         ranges=IRanges(start=0, end=1000000000))
  altassembly$name <- as.character(seqnames(altassembly))
  altassembly$score <- 1000
  #ucscgapblacklist$score <- 1000
  blacklist <- as(rbind(
    #  as(ucscgapblacklist, "RangedData"),
    #  as(mappabilityblacklist, "RangedData"),
    as(encodeblacklist, "RangedData"),
    as(altassembly, "RangedData")), "GRanges")
}

#####################################
# na12878 Moleculo/PacBio truth
if (!exists("srpacbio") || !exists("srmoleculo") || !exists("sppacbio") || !exists("spmoleculo")) {
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
}
# Split read validation
countBedHits <- function(gr, bed) {
  hits <- data.table(findOverlaps_type_equal_df(gr, bed, maxgap=maxPositionDifference_LengthOrEndpoint))
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
longReadBed <- function(vcflist) {
  # extract deletion calls from VCFs
  dfdelcalls <- rbindlist(lapply(vcflist, function(vcf) {
    caller <- attr(vcf, "sourceMetadata")$CX_CALLER
    method <- as.character(attr(vcf, "sourceMetadata")$CX_ASSEMBLY_METHOD)
    kmer <- as.character(attr(vcf, "sourceMetadata")$CX_K)
    model <- as.character(attr(vcf, "sourceMetadata")$CX_MODEL)
    exclusion <- as.character(attr(vcf, "sourceMetadata")$CX_EXCLUSION)
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
      assembly=as.character(rep(ifelse(is.null(method), "", method), length(gro))),
      kmer=as.character(rep(ifelse(is.null(kmer), "", kmer), length(gro))),
      model=as.character(rep(ifelse(is.null(model), "", model), length(gro))),
      exclusion=as.character(rep(ifelse(is.null(exclusion), "", exclusion), length(gro))),
      row.names=NULL)
    return(df)
  }))
  #write.table(dfdelcalls, paste0(rootdir, "i/data.na12878/tovalidate.bedpe"), sep='\t', quote=FALSE, row.names=FALSE)
  delbed <- GRanges(seqnames=dfdelcalls$chrom2,
                    ranges=IRanges(start=pmin(dfdelcalls$pos1, dfdelcalls$pos2), end=pmax(dfdelcalls$pos1, dfdelcalls$pos2)),
                    caller=dfdelcalls$name,
                    QUAL=dfdelcalls$score,
                    length=dfdelcalls$length,
                    filtered=dfdelcalls$filtered,
                    id=dfdelcalls$id,
                    assembly=dfdelcalls$assembly,
                    kmer=dfdelcalls$kmer,
                    model=dfdelcalls$model,
                    exclusion=dfdelcalls$exclusion)
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
  return(delbed)
}
bedToROC <- function(delbed) {
  longReadBed_delroc <- data.table(as.data.frame(mcols(delbed)))
  longReadBed_delroc <- longReadBed_delroc[!longReadBed_delroc$blacklisted,]
  tmp <- longReadBed_delroc[!longReadBed_delroc$filtered & longReadBed_delroc$caller %in% unique(longReadBed_delroc[longReadBed_delroc$filtered,]$caller),]
  tmp$filtered <- TRUE
  longReadBed_delroc <- rbind(longReadBed_delroc, tmp)
  longReadBed_delroc <- longReadBed_delroc[order(-longReadBed_delroc$QUAL),]
  longReadBed_delroc_bylist = c("caller", "filtered", "assembly", "kmer", "model", "exclusion")
  longReadBed_delroc[,`:=`(tp=cumsum(tp), fp=cumsum(fp), QUAL=cummin(QUAL)), longReadBed_delroc_bylist]
  longReadBed_delroc <- longReadBed_delroc[!duplicated(longReadBed_delroc[, c(longReadBed_delroc_bylist, "QUAL"), with=FALSE], fromLast=TRUE),] # take only one data point per QUAL
  longReadBed_delroc$Filter <- ifelse(longReadBed_delroc$filtered, "Including Filtered", "Default")
  longReadBed_delroc$precision <- longReadBed_delroc$tp/(longReadBed_delroc$tp+longReadBed_delroc$fp)
  return (longReadBed_delroc)
}