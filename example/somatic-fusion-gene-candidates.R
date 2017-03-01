#
# Very basic example R script that demonstrates
# how GRIDSS output can be combined with the 
# StructuralVariantAnnotation package and BioConductor
# to perform useful analyses.
#
# This script performs a very basic check for somatic
# gene fusions that could result in a fusion transcript
# It does not check that the transcript in in-frame, nor
# does it check that the resultant fusion actually involves
# one or more exons from each gene.
#
# CRAN packages
library(devtools)
library(stringr)
library(dplyr)
library(ggplot2)
# bioconductor packages
library(VariantAnnotation)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(StructuralVariantAnnotation) # install_github("d-cameron/StructuralVariantAnnotation")

vcf <- readVcf("normal-tumour.sv.vcf", "hg19")

# Remove low confidence variants
vcf <- vcf[rowRanges(vcf)$QUAL >= 500,]

# very simple somatic filter: normal support less than 5% that of the tumour
# this assumes that the first INPUT was the normal
vcf <- vcf[geno(vcf)$QUAL[,1] < 0.05 * rowRanges(vcf)$QUAL,]

# convert to breakend GRanges
gr <- breakpointRanges(vcf)

# remove events under 10kbp
#gr <- gr[seqnames(gr) != seqnames(partner(gr)) | abs(start(gr) - start(partner(gr))) > 10000]

# retain only primary chromosomes
seqlevelsStyle(gr) <- "UCSC"
gr <- gr[seqnames(gr) %in% paste0("chr", c(1:22, "X", "Y"))]
seqlevels(gr) <- paste0("chr", c(1:22, "X", "Y"))

# remove breakends that now don't have a partner (eg: chr1 -> chrMT)
gr <- gr[gr$partner %in% names(gr)]

# annotate breakends with gene names and gene orientation
gns <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
hits <- as.data.frame(findOverlaps(gr, gns, ignore.strand=TRUE))
hits$SYMBOL <- biomaRt::select(org.Hs.eg.db, gns[hits$subjectHits]$gene_id, "SYMBOL")$SYMBOL
hits$gene_strand <- as.character(strand(gns[hits$subjectHits]))
hits <- hits %>%
  group_by(queryHits) %>%
  summarise(SYMBOL=paste(SYMBOL, collapse=","), gene_strand=paste0(gene_strand, collapse=""))
gr$SYMBOL <- ""
gr$geneStrand <- ""
gr$SYMBOL[hits$queryHits] <- hits$SYMBOL
gr$geneStrand[hits$queryHits] <- hits$gene_strand

# require the breakpoint to be between different genes
gr <- gr[gr$SYMBOL != partner(gr)$SYMBOL,]
gr <- gr[gr$SYMBOL != "" & partner(gr)$SYMBOL != "",]

# require the breakpoint to possibly generate a fusion transcript
gr$couldBeThreePrimeStart <- str_detect(gr$geneStrand, stringr::fixed(as.character(strand(gr))))
gr$couldBeFivePrimeEnd <- str_detect(gr$geneStrand, stringr::fixed(ifelse(as.character(strand(gr))=="+", "-", "+")))
gr <- gr[(gr$couldBeThreePrimeStart & partner(gr)$couldBeFivePrimeEnd) |
		(gr$couldBeFivePrimeEnd & partner(gr)$couldBeThreePrimeStart),]

# return highest scoring variants first
gr <- gr[order(-gr$QUAL)]


