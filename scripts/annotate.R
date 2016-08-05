#library(StructuralVariantAnnotation)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
#library(BSgenome.Hsapiens.UCSC.hg19)


annotate_genes <- function(gr) {
  var <- locateVariants(gr, TxDb.Hsapiens.UCSC.hg19.knownGene, AllVariants())
  var <- var[!is.na(var$GENEID)] # remove intergenic hits
  var$SYMBOL <- select(org.Hs.eg.db, na.omit(var$GENEID), c("SYMBOL"), keytype="ENTREZID")$SYMBOL
  genes <- rep(NA, length(gr))
  lookup <- aggregate(SYMBOL ~ QUERYID, mcols(var), function(x) paste0(unique(x)))
  genes[var$QUERYID] <- var$SYMBOL
  return(genes)
}


vcf <- readVcf("W:/projects/DREAM/samples/synthetic6-bwa/synthetic6.gridss.vcf", "hg19")
# apply filtering criteria to reduce the noise
vcf <- vcf[rowRanges(vcf)$QUAL >= 500,]

# filter out alternate contigs
gr <- rowRanges(vcf)
seqlevelsStyle(gr) <- "UCSC"
vcf <- vcf[seqnames(gr) %in% paste0("chr", c(1:22, "X", "Y")),]
vcf <- vcf[names(vcf) %in% info(vcf)$PARID,] # remove the dangling primary breakend for primary-alt breakpoint calls

# annotate gene names (if breakend is anywhere in the gene: promoter, both UTRs, split, or coding)
gr <- rowRanges(vcf)
seqlevelsStyle(gr) <- "UCSC"
seqlevels(gr) <- paste0("chr", c(1:22, "X", "Y")) # hack for hs37d5 having different MT length
genenames <- annotate_genes(gr)
genenames[is.na(genenames)] <- ""
info(vcf)$GENE <- genenames
info(vcf)$FUSION <- paste0(info(vcf)$GENE, "/", info(vcf)[info(vcf)$PARID,]$GENE)

