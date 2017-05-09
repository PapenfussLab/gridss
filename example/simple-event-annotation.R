#source("https://bioconductor.org/biocLite.R")
#biocLite("VariantAnnotation")
#install_github("PapenfussLab/StructuralVariantAnnotation")
#install.packages("stringr")
library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(stringr)
#' Simple SV type classifier
simpleEventType <- function(gr) {
  return(ifelse(seqnames(gr) != seqnames(partner(gr)), "ITX", # inter-chromosomosal
          ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS", # TODO: improve classification of complex events
           ifelse(strand(gr) == strand(partner(gr)), "INV",
            ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "DEL",
             "DUP")))))
}
# using the example in the GRIDSS /example directory
vcf <- readVcf("chr12.1527326.DEL1024.sv.vcf", "hg19")
gr <- breakpointRanges(vcf)
svtype <- simpleEventType(gr)
info(vcf)$SIMPLE_TYPE <- NA_character_
info(vcf[gr$vcfId])$SIMPLE_TYPE <- svtype
info(vcf[gr$vcfId])$SVLEN <- gr$svLen
writeVcf(vcf, "chr12.1527326.DEL1024.sv.annotated.vcf")

# TODO: perform event filtering here
# By default, GRIDSS is very sensitive but this comes at the cost of a high false discovery rate
gr <- gr[gr$FILTER == "." & partner(gr)$FILTER == "."] # Remove low confidence calls

simplegr <- gr[simpleEventType(gr) %in% c("INS", "INV", "DEL", "DUP")]
simplebed <- data.frame(
    chrom=seqnames(simplegr),
	# call the centre of the homology/inexact interval
    start=as.integer((start(simplegr) + end(simplegr)) / 2),
    end=as.integer((start(partner(simplegr)) + end(partner(simplegr))) / 2),
    name=simpleEventType(simplegr),
    score=simplegr$QUAL,
    strand="."
    )
# Just the lower of the two breakends so we don't output everything twice
simplebed <- simplebed[simplebed$start < simplebed$end,]
write.table(simplebed, "chr12.1527326.DEL1024.simple.bed", quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
