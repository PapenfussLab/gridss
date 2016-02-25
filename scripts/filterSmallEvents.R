#!/usr/bin/env Rscript
#
# Subsets the calls by removing all calls with an implied event size smaller than the given size
# This corresponds to the GRIDSS varantcalling.minSize parameter
#
usage <- "Usage: filterSmallEvents.R <minimum size> <input VCF> <output VCF>"
# TODO: fix this by creating a proper R package
libgridssPath <- "~/dev/gridss/src/test/r/libgridss.R"
if (!file.exists(libgridssPath)) {
  write(paste("Unable to find", libgridssPath, " Please update excludeCategory.R with the location of your GRIDSS source."), stderr())
  q(save="no", status=1)
}
args <- commandArgs(TRUE)
if (length(args) != 3) {
  write(usage, stderr())
  q(save="no", status=1)
}

source(libgridssPath, chdir=TRUE)

if (!exists("gridss.vcftodf")) {
  write(paste("Unable to find libgridss.R. This should be located at", scriptPath, "/", relativeLocation), stderr())
  q(save="no", status=1)
}
vcf <- readVcf(args[2], "unknown")
bpgr <- vcftobpgr(vcf)
writeVcf(vcf[unique(bpgr$vcfIndex[is.na(bpgr$size) | bpgr$size>=32])], args[3])

#q(save="no", status=0)