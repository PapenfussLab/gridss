library(openxlsx) # install.packages("openxlsx")
library(VariantAnnotation)

# extracts read pair breakpoint calls from the supplementary table
# creates a breakpoint GRanges object containing
# breakend positions
# .mate column containing name of corresponding breakend in GRanges object
getrpcalls <- function(xlsx = "mmc3.xlsx", tabname) {
  dt <- read.xlsx(xlsx, tabname)
  row.names(dt) <- paste0(dt$chrom1, ":", dt$start1, dt$strand1, dt$chrom2, ":", dt$start2, dt$strand2)
  gr <- GRanges(seqnames=c(dt$chrom1, dt$chrom2),
                ranges=IRanges(start=c(dt$start1, dt$start2), width=1),
                strand=c(as.character(dt$strand1), as.character(dt$strand2)),
                mate=c(paste0(row.names(dt),"/1"), paste0(row.names(dt),"/2")),
                nreads=dt$nreads)
  names(gr) <- c(paste0(row.names(dt),"/2"), paste0(row.names(dt),"/1"))
  return(gr)
}
# extract CGR regions from the supplementary table
getcgr <- function(xlsx = "mmc4.xlsx", tabname) {
  dt <- read.xlsx(xlsx, tabname)
  gr <- GRanges(seqnames=dt$Chromosome, ranges=IRanges(start=dt$Start, end=dt$End))
  names(gr) <- dt$CGR.name
  return(gr)
}
# extract copy number from the supplementary table
getcn <- function(xlsx = "mmc4.xlsx", tabname) {
  dt <- read.xlsx(xlsx, tabname)
  gr <- GRanges(seqnames=dt$Chromosome,
                ranges=IRanges(start=dt$Start, end=dt$End),
                cgrname=dt$CGR.name,
                cn=dt$Copy.number)
  return(gr)
}

# counts overlapping breakpoints
countBreakpointHits <- function(queryGr, mateQueryGr, subjectGr, mateSubjectGr, ...) {
  dfhits <- rbind(as.data.frame(findOverlaps(queryGr, subjectGr, ...)),
                  as.data.frame(findOverlaps(mateQueryGr, mateSubjectGr, ...)))
  dfhits <- dfhits[duplicated(dfhits),] # both breakends match
  queryHits <- rep(0, length(queryGr))
  queryHits[count(dfhits, "queryHits")$queryHits] <- count(dfhits, "queryHits")$freq
  subjectHits <- rep(0, length(subjectGr))
  subjectHits[count(dfhits, "subjectHits")$subjectHits] <- count(dfhits, "subjectHits")$freq
  return(list(queryHits=queryHits, subjectHits=subjectHits))
}
# counts overlapping breakpoints
# VCF records MUST each be a valid symbolic breakend alleles
# VCF MUST have a single valid MATEID for every record
# GRanges MUST have a $mate set to the name of the other breakend for that breakpoint
countVcfGrBreakpointHits <- function(vcf, gr, ...) {
  vcfgr <- rowData(vcf)
  # set strand
  strand(vcfgr) <- ifelse(str_detect(as.character(fixed(vcf)$ALT), "[[:alpha:]]+(\\[|]).*(\\[|])"), "+", "-")
  return(countBreakpointHits(vcfgr, vcfgr[as.character(info(vcf)$MATEID),],
                             gr, gr[gr$mate,], ...))
}

toroc <- function(qual, istp, totaltp=sum(istp)) {
  orderedistp <- istp[order(-qual)]
  orderedqual <- qual[order(-qual)]
  return(data.frame(qual=orderedqual, sens=cumsum(orderedistp) / totaltp))
}






































