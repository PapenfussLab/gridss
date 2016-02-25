source("libvcf.R")


load_grdGV <- function(file) {
  dtdGV <- read.csv(file, sep="\t", stringsAsFactors=FALSE)
  grdGV <- GRanges(seqnames=dtdGV$chr, ranges=IRanges(start=dtdGV$start, end=dtdGV$end, name=dtdGV$variantaccession))
  return(grdGV)
}

findOvleraps_bpgr_dGV <- function(bpgr, grdGV) {
  if (class(grdGV)[1] != "GRanges") {
    stop("dGV must be supplied in GRanges format")
  } else if (class(bpgr)[1] != "GRanges") {
    stop("bpgr not a GRanges object")
  }
}

findOvleraps_bpgr_interval <- function(bpgr, gr) {
  if (class(gr)[1] != "GRanges") {
    stop("gr not a GRanges object")
  } else if (class(bpgr)[1] != "GRanges") {
  }
}