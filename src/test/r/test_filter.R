setwd("W:/dev/gridss/src/test/r")
source("libfilter.R")
library("GenomeInfoDb")

grdGV <- load_grdGV("../../../../../projects/reference_datasets/human/dGV.GRCh37_hg19_variants_2014-10-16.txt")
grdGV <- renameSeqlevels(grdGV, mapSeqlevels(seqlevels(grdGV), "UCSC"))
