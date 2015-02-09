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
                mate=c(paste0(row.names(dt),"/1"), paste0(row.names(dt),"/2")))
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



zz <- getrpcalls("C:/dev/neochromosome/mmc3.xlsx", "778 (DR)")
zz <- getcgr("C:/dev/neochromosome/mmc4.xlsx", "778_CGRs")
zz <- getcn("C:/dev/neochromosome/mmc4.xlsx", "778_CN")
dt <- read.xlsx2(xlsx, tabname, colClasses=c("character", "numeric", "numeric", "character", "numeric"))

xlsx <- "C:/dev/neochromosome/mmc4.xlsx"
tabname <- "778_CGRs"

