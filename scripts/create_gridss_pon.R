#!/usr/bin/env Rscript
#
# Incorporates the given VCF files into the Panel Of Normals
#
library(argparser, quietly=TRUE)
argp = arg_parser("Filters a raw GRIDSS VCF into somatic call subsets.")
argp = add_argument(argp, "--pondir", default=NA, help="Directory to write PON to.")
argp = add_argument(argp, "--scriptdir", default=ifelse(sys.nframe() == 0, "./", dirname(sys.frame(1)$ofile)), help="Path to libgridss.R script")
argp = add_argument(argp, "--cache-only", flag=TRUE, help="Only generate cache objects for input files. Useful for parallel processing of inputs.")
argp = add_argument(argp, "--normalordinal", type="integer", default=1, help="Ordinal of normal sample in the VCF")
argp = add_argument(argp, "--batchsize", default=100, help="Number of VCFs to process")
argp = add_argument(argp, "--input", nargs=Inf, help="Input VCFs normal")
# argv = parse_args(argp, argv=c("--pondir", "D:/hartwig/pon/", "--scriptdir", "D:/hartwig/scripts/gridss", "--batchsize", "1", "--input", "D:/hartwig/down/COLO829R_COLO829T.gridss.vcf", "D:/hartwig/down/DRUP01010047R_DRUP01010047T.gridss.vcf"))
argv = parse_args(argp)

if (!dir.exists(argv$pondir)) {
  write("PON directory not found", stderr())
  q(save="no", status=1)
}
if (!all(file.exists(argv$input))) {
  write("VCF input file not found", stderr())
  q(save="no", status=1)
}
libgridssfile = paste0(argv$scriptdir, "/", "libgridss.R")
if (file.exists(libgridssfile)) {
  tmpwd = getwd()
  setwd(argv$scriptdir)
  library(tidyverse, quietly=TRUE)
  library(readr, quietly=TRUE)
  library(stringr, quietly=TRUE)
  library(rtracklayer, quietly=TRUE)
  library(R.cache, quietly=TRUE)
  source("libgridss.R")
  setwd(tmpwd)
} else {
  msg = paste("Could not find libgridss.R in", argv$scriptdir, " - please specify a --scriptdir path to a directory containing the required scripts")
  write(msg, stderr())
  print(argp)
  stop(msg)
}

samplefile = paste(argv$pondir, "gridss_pon_samples.txt", sep="/")
befile = paste(argv$pondir, "gridss_pon_single_breakend.bed", sep="/")
beimpfile = paste(argv$pondir, "gridss_pon_single_breakend_imprecise.bed", sep="/")
bpfile = paste(argv$pondir, "gridss_pon_breakpoint.bedpe", sep="/")

setCacheRootPath(paste0(argv$pondir, "/Rcache"))
options("R.cache.compress"=TRUE)
load_germline_pon_calls = function(vcf_file, sampleId) {
  vcf_file = normalizePath(vcf_file)
  key=list(vcf_file=vcf_file)
  cached = loadCache(key=key)
  if (!is.null(cached)) {
    write(paste("Using cached data for ", vcf_file), stderr())
    return(cached)
  }
  write(paste("Start load", vcf_file), stderr())
  full_vcf = readVcf(vcf_file, "hg19")
  bpgr = breakpointRanges(full_vcf, unpartneredBreakends=FALSE)
  begr = breakpointRanges(full_vcf, unpartneredBreakends=TRUE)

  bpgr = bpgr[geno(full_vcf[bpgr$sourceId])$QUAL[,argv$normalordinal] > gridss.pon.min_normal_qual | geno(full_vcf[bpgr$partner])$QUAL[,argv$normalordinal] > gridss.pon.min_normal_qual]
  begr = begr[geno(full_vcf[begr$sourceId])$BQ[,argv$normalordinal] > gridss.pon.min_normal_qual * gridss.single_breakend_multiplier]

  minimal_bpgr = bpgr
  mcols(minimal_bpgr) = NULL
  minimal_bpgr$vcf = rep(sampleId, length(bpgr))
  minimal_bpgr$IMPRECISE = info(full_vcf[names(minimal_bpgr)])$IMPRECISE
  names(minimal_bpgr) = paste(minimal_bpgr$vcf, names(minimal_bpgr), sep="_")
  minimal_bpgr$partner = paste(minimal_bpgr$vcf, bpgr$partner, sep="_")

  minimal_begr = begr
  mcols(minimal_begr) = NULL
  minimal_begr$vcf = rep(sampleId, length(begr))
  minimal_begr$IMPRECISE = info(full_vcf[names(minimal_begr)])$IMPRECISE
  names(minimal_begr) = NULL
  result = list(bp=minimal_bpgr, be=minimal_begr)
  saveCache(result, key=key)
  write(paste("End load", vcf_file), stderr())
  return(result)
}
load_pon = function() {
  for (f in c(samplefile, befile, beimpfile, bpfile)) {
    if (!file.exists(f)) {
      write(paste("Missing pon file", f, "creating new empty file."), stderr())
      file.create(f)
    }
  }
  pon = list(
    samples = read_tsv(samplefile, col_names=c("sample"), col_types=c("c")),
    bedf = read_tsv(befile, col_names=c("chrom", "chromStart", "chromEnd", "name", "score", "strand"), col_types=c("ciicic")),
    beimpdf = read_tsv(beimpfile, col_names=c("chrom", "chromStart", "chromEnd", "name", "score", "strand"), col_types=c("ciicic")),
    bpdf = read_tsv(bpfile, col_names=c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2", "IMPRECISE"), col_types=c("ciiciiciccl")))
  return(pon)
}
save_pon = function(pon) {
  withr::with_options(c(scipen = 10), write.table(pon$samples, file=samplefile, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE))
  withr::with_options(c(scipen = 10), write.table(
    pon$bedf %>% arrange(chrom, chromStart, chromEnd, strand),
    file=befile, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE))
  withr::with_options(c(scipen = 10), write.table(
    pon$beimpdf %>% arrange(chrom, chromStart, chromEnd, strand),
    file=beimpfile, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE))
  withr::with_options(c(scipen = 10), write.table(
    pon$bpdf %>% arrange(chrom1, start1, end1, chrom2, start2, end2, strand1, strand2),
    file=bpfile, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE))
}
as_beddf = function(begr) {
  data.frame(
    chrom = as.character(seqnames(begr)),
    chromStart = start(begr) - 1,
    chromEnd = end(begr),
    name = ".",
    score = 1,
    strand = as.character(strand(begr)),
    stringsAsFactors=FALSE)
}
as_bedpedf = function(bpgr) {
  lowergr = bpgr[str_detect(names(bpgr), "o$")]
  uppergr = bpgr[lowergr$partner]
  data.frame(
    chrom1=as.character(GenomeInfoDb::seqnames(lowergr)),
    start1=start(lowergr) - 1,
    end1=end(lowergr),
    chrom2=as.character(GenomeInfoDb::seqnames(uppergr)),
    start2=start(uppergr) - 1,
    end2=end(uppergr),
    name=".",
    score=1,
    strand1=as.character(strand(lowergr)),
    strand2=as.character(strand(uppergr)),
    IMPRECISE=lowergr$IMPRECISE,
    stringsAsFactors=FALSE)
}
merge_pon_be = function(...) {
  bind_rows(...) %>%
    group_by(chrom, chromStart, chromEnd, strand) %>%
    summarise(score=sum(score)) %>%
    ungroup() %>%
    mutate(name=".") %>%
    dplyr::select(chrom, chromStart, chromEnd, name, score, strand)
}
merge_pon_bp = function(...) {
  bind_rows(...) %>%
    group_by(chrom1, start1, end1, chrom2, start2, end2, strand1, strand2, IMPRECISE) %>%
    summarise(score=sum(score)) %>%
    ungroup() %>%
    mutate(name=".") %>%
    dplyr::select(chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2, IMPRECISE)
}
merge_pon_calls = function(pon, samples, bedf, beimpdf, bpdf) {
  assertthat::assert_that(!any(pon$samples$samples %in% samples))
  pon$samples = bind_rows(pon$samples, samples)
  pon$bedf = merge_pon_be(pon$bedf, bedf)
  pon$beimpdf = merge_pon_be(pon$beimpdf, beimpdf)
  pon$bpdf = merge_pon_bp(pon$bpdf, bpdf)
  return(pon)
}

###
# Start script
###
if (!argv$cache_only) {
  pon = load_pon()
}
for (chunk_files in split(argv$input, ceiling(seq_along(argv$input)/argv$batchsize))) {
  full_bp = list()
  full_be = list()
  samples = list()
  for (vcf_file in chunk_files) {
    sampleId = str_replace(str_replace(basename(vcf_file), ".gridss.vcf.gz", ""), ".gridss.vcf", "")
    if (!is.null(pon$samples$samples) & sampleId %in% pon$samples$samples) {
      write(paste("Skipping", sampleId, "already in PON."), stderr())
    } else {
      calls = load_germline_pon_calls(vcf_file, sampleId)
      samples[[sampleId]] = data.frame(samples=sampleId, stringsAsFactors=FALSE)
      full_bp[[sampleId]] = calls$bp
      full_be[[sampleId]] = calls$be
    }
  }
  if (!argv$cache_only && length(full_bp) > 0) {
    write(paste("Updating PON with", length(full_bp), "new samples"), stderr())
    pon = merge_pon_calls(
      pon,
      bind_rows(samples),
      bind_rows(lapply(full_be, function(x) {as_beddf(x[!x$IMPRECISE])})),
      bind_rows(lapply(full_be, function(x) {as_beddf(x[x$IMPRECISE])})),
      bind_rows(lapply(full_bp, as_bedpedf)))
    save_pon(pon)
  }
}
write("Complete", stderr())
