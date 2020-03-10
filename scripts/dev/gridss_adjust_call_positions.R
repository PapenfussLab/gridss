#!/usr/bin/env Rscript
library(argparser)
argp = arg_parser("Filters a raw GRIDSS VCF into somatic call subsets.")
argp = add_argument(argp, "--input", help="GRIDSS VCF")
argp = add_argument(argp, "--output", help="High confidence somatic subset")
argp = add_argument(argp, "--scriptdir", default=ifelse(sys.nframe() == 0, "./", dirname(sys.frame(1)$ofile)), help="Path to libgridss.R script")
argv = parse_args(argp)

if (!file.exists(argv$input)) {
  msg = paste(argv$input, "not found")
  write(msg, stderr())
  print(argp)
  stop(msg)
}
library(tidyverse)
library(readr)
library(stringr)
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
write(paste0("Reading ", argv$input), stderr())
full_vcf = readVcf(argv$input)
full_vcf = align_breakpoints(full_vcf)
write(paste0("Writing ", argv$output), stderr())
writeVcf(full_vcf, argv$output, index=TRUE)



