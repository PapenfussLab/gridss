#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)

library(argparser)
argp = arg_parser("Annotates insertion sequence with the longest alignment overlaping a repeatmasker repeat.")
argp = add_argument(argp, "--input", help="Input GRIDSS VCF")
argp = add_argument(argp, "--output", help="Output GRIDSS VCF")
argp = add_argument(argp, "--repeatmasker", help="RepeatMasker .fa.out file")
argp = add_argument(argp, "--scriptdir", default=ifelse(sys.nframe() == 0, "./", dirname(sys.frame(1)$ofile)), help="Path to libgridss.R script")
# argv = parse_args(argp, c("--input", "D:/hartwig/down/pre.vcf", "--output", "D:/hartwig/down/testrm.vcf","--repeatmasker", "D:/hartwig/hg19.fa.out"))
argv = parse_args(argp)

repeat_notes = data.frame(
  repeatType=c("HSATII", "(GAATG)n", "(CATTC)n", "ALR/Alpha", "(CCCTAA)n", "(TTAGGG)n", "TAR1"),
  annotation=c("centromere", "centromere", "centromere", "centromere", "telomere", "telomere", "telomere"))

if (!file.exists(argv$input)) {
  msg = paste(argv$input, "not found")
  write(msg, stderr())
  print(argp)
  stop(msg)
}
if (!file.exists(argv$repeatmasker)) {
  msg = paste(argv$repeatmasker, "not found")
  write(msg, stderr())
  print(argp)
  stop(msg)
}
libgridssfile = paste0(argv$scriptdir, "/", "libgridss.R")
if (file.exists(libgridssfile)) {
  tmpwd = getwd()
  setwd(argv$scriptdir)
  source("libgridss.R")
  setwd(tmpwd)
} else {
  msg = paste("Could not find libgridss.R in", argv$scriptdir, " - please specify a --scriptdir path to a directory containing the required scripts")
  write(msg, stderr())
  print(argp)
  stop(msg)
}
library(VariantAnnotation)
library(tidyverse)
library(stringr)

# from http://github.com/PapenfussLab/sv_benchmark
import.repeatmasker.fa.out <- function(repeatmasker.fa.out) {
  rmdt <- read_table2(repeatmasker.fa.out, col_names=FALSE, skip=3)
  grrm <- GRanges(
    seqnames=rmdt$X5,
    ranges=IRanges(start=rmdt$X6 + 1, end=rmdt$X7),
    strand=ifelse(rmdt$X9=="C", "-", "+"),
    repeatType=rmdt$X10,
    repeatClass=rmdt$X11)
  return(grrm)
}
# wget http://www.repeatmasker.org/genomes/hg19/RepeatMasker-rm330-db20120124/hg19.fa.out.gz
cache_filename = paste0(argv$repeatmasker, ".grrm.rds")
if (file.exists(cache_filename)) {
  grrm = readRDS(cache_filename)
} else {
  grrm = import.repeatmasker.fa.out(argv$repeatmasker)
  saveRDS(grrm, file=cache_filename)
}
seqlevelsStyle(grrm) = "NCBI"

####
# Start processing

vcf = readVcf(argv$input)
# Add new fields to header
info(header(vcf)) = unique(as(rbind(as.data.frame(info(header(vcf))), data.frame(
  row.names=c("INSRMRT", "INSRMRC", "INSRMRO", "INSRMP"),
  Number=c("1", "1", "1", "1"),
  Type=c("String", "String", "String", "Float"),
  Description=c(
    "Inserted sequence repeatmasker repeat type",
    "Inserted sequence repeatmasker repeat class",
    "Inserted sequence repeatmasker repeat orientation",
    "Portion of inserted sequence whose alignment overlaps the repeatmasker repeat. 1.0 indicates the inserted sequence entirely mapping to the repeat."))), "DataFrame"))

insseqgr = with(data.frame(
    sourceId=rep(names(rowRanges(vcf)), lengths(info(vcf)$BEALN)),
    BEALN=unlist(info(vcf)$BEALN)) %>%
  separate(BEALN, sep="[:|]", into=c("chr", "start", "orientation", "cigar", "maqp")) %>%
  mutate(
    start=as.integer(start),
    end=start+GenomicAlignments::cigarWidthAlongReferenceSpace(cigar)),
  GRanges(seqnames=chr, ranges=IRanges(start=start, end=end), strand=orientation, sourceId=sourceId))
hits = findOverlaps(insseqgr, grrm, select="all", ignore.strand=TRUE)
hits = hits %>% as.data.frame() %>%
  mutate(
    overlap=pmin(end(insseqgr[queryHits]), end(grrm[subjectHits])) - pmax(start(insseqgr[queryHits]), start(grrm[subjectHits])),
    repeatOverlapPercentage=overlap/(end(insseqgr[queryHits])-start(insseqgr[queryHits])),
    repeatType=grrm[subjectHits]$repeatType,
    repeatClass=grrm[subjectHits]$repeatClass,
    sourceId=insseqgr[queryHits]$sourceId,
    ALT=as.character(rowRanges(vcf)[sourceId]$ALT),
    repeat_orientation=as.character(strand(grrm[subjectHits])))
insrmdf = hits %>%
  group_by(sourceId) %>%
  top_n(1, repeatOverlapPercentage) %>%
  distinct(sourceId, .keep_all=TRUE) %>%
  ungroup() %>%
  dplyr::select(sourceId, repeatType, repeatClass, repeat_orientation, repeatOverlapPercentage)

info(vcf)$INSRMRT=NA_character_
info(vcf)$INSRMRC=NA_character_
info(vcf)$INSRMRO=NA_character_
info(vcf)$INSRMP=NA_real_

info(vcf[insrmdf$sourceId])$INSRMRT=insrmdf$repeatType
info(vcf[insrmdf$sourceId])$INSRMRC=insrmdf$repeatClass
info(vcf[insrmdf$sourceId])$INSRMRO=insrmdf$repeat_orientation
info(vcf[insrmdf$sourceId])$INSRMP=insrmdf$repeatOverlapPercentage

writeVcf(vcf, argv$output, index=TRUE)





