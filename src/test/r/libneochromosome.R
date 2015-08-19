library(openxlsx) # install.packages("openxlsx")
library(VennDiagram) # install.packages("VennDiagram")
library(ggplot2)
library(scales)
library(RColorBrewer)
library(stringr)
library(rtracklayer)
source("../../main/r/libgridss.R")
source("libvcf.R")

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
vreplaceNA <- function(vector, value) {
  vector[is.na(vector)] <- value[is.na(vector)]
  return(vector)
}
# matchmaxgap: 778 #chr1:188377591-chr1:188379026- call location is over 120
go <- function(sample, vcf, rp, filterBed=NULL, minimumEventSize, cgrmaxgap=1000, matchmaxgap=200) {
  if (!is.null(filterBed)) {
    # filter so that both ends are within the filterBed
    rp <- subsetByOverlaps(rp, filterBed, maxgap=cgrmaxgap)
    rp <- rp[rp$mate %in% names(rp)]
    vcf <- vcf[overlapsAny(rowRanges(vcf), filterBed, maxgap=cgrmaxgap),]
  }
  vcf <- gridss.removeUnpartnerededBreakend(vcf)
  rpMate <- rp[rp$mate,]
  matches <- gridss.annotateBreakpoints(rp, rpMate, vcf, maxgap=matchmaxgap)
  vcfdf <- matches$gridss
  rp <- matches$bed
  
  rp$hits <- vcfdf[rp$gridssid,]$confidence
  levels(rp$hits) <- c("Low", "Medium", "High", "None")
  rp$hits[is.na(rp$hits)] <- "None"
  rp$spanning <- FALSE
  
  # filter out small events that don't hit an existing call
  vcf <- vcf[is.na(vcfdf$size) | vcfdf$size >= minimumEventSize,]
  vcfdf <- vcfdf[rownames(vcf),]
  
  spanning <- FindFragmentSpanningEvents(rp, vcftobpgr(vcf))
  rp[spanning$query,]$spanning <- TRUE
  rp[rp[spanning$query,]$mate,]$spanning <- TRUE
  
  rp$distanceToMedHigh <- distanceToClosest(rp, rowRanges(vcf[vcfdf$confidence != "Low"]))
  vcfdf$distanceToCall <- distanceToClosest(rowRanges(vcf), rp)
  vcfdf$found <- !is.na(vcfdf$bedid)
  
  rp$mhcallsWithin1kbp <- countOverlaps(rp, rowRanges(vcf[vcfdf$confidence %in% c("High", "Medium")]), maxgap=1000)
  vcfdf$mhcallsWithin1kbp <- countOverlaps(vcf, rowRanges(vcf[vcfdf$confidence %in% c("High", "Medium")]), maxgap=1000)
  
  ###############
  # Variant calling concordance
  ###############
  # VennDiagram requires lists of ids
  # We'll go with the bedid if it exists, then fall back to the gridssid
  venn.plot <- venn.diagram(
    x = list(
      Published=names(rp),
      High=vreplaceNA(vcfdf[vcfdf$confidence=="High",]$bedid, rownames(vcfdf[vcfdf$confidence=="High",])),
      Medium=vreplaceNA(vcfdf[vcfdf$confidence=="Medium",]$bedid, rownames(vcfdf[vcfdf$confidence=="Medium",])),
      Low=vreplaceNA(vcfdf[vcfdf$confidence=="Low",]$bedid, rownames(vcfdf[vcfdf$confidence=="Low",]))),
    filename = paste0("venn_", sample, ".png"))
    
  if (!file.exists(paste0(sample, ".additional.high.vcf"))) {
    # don't rewrite until writeVcf trailing tab bug is fixed
    writeVcf(vcf[is.na(vcfdf$bedid) & vcfdf$confidence=="High",], paste0(sample, ".additional.high.vcf"))
    writeVcf(vcf[is.na(vcfdf$bedid) & vcfdf$confidence=="Medium",], paste0(sample, ".additional.medium.vcf"))
  }
  write.csv(vcfdf[is.na(vcfdf$bedid) & vcfdf$confidence %in% c("High", "Medium"),], paste0(sample, ".additional.csv"))
  export(rp[rp$hits %in% c("None", "Low"),], paste0(sample, ".misses.bed"))
  write.csv(as.data.frame(rp[rp$hits %in% c("None", "Low"),]), paste0(sample, ".missed.csv"))
  
  
  ###############
  # Read Pair variant calling concordance
  ###############
  #library(pROC) #install.packages("pROC")
#   dfroc <- NULL
#   for (confidence in as.numeric(unique(vcfdf$confidence))) {
#     subset <- as.numeric(vcfdf$confidence) >= confidence
#     dftmp <- vcftoroc(vcf[subset], rp, maxgap=rpmaxgap)
#     dftmp$confidence <- levels(vcfdf$confidence)[confidence]
#     dfroc <- rbind(dfroc, dftmp)  
#   }
#   ggplot(dfroc[dfroc$qual!=0,], aes(x=log(qual+1), y=sens, color=confidence)) +
#     geom_line() +
#     scale_x_reverse() +
#     theme_bw() +
#     scale_y_continuous(limits=c(0,1)) +
#     labs(title=paste0("Sensitivity of curated RP call detection - ", sample))
#   ggsave(paste0("rp_roc_", sample, ".png"), width=10, height=7.5)
  
  ggplot(as.data.frame(mcols(rp)), aes(x=nreads, fill=hits)) +
    geom_histogram() +
    scale_x_log10() + 
    labs(title=paste0("Sensitivity of curated RP call detection - ", sample))
  ggsave(paste0("rp_hist_", sample, ".png"), width=10, height=7.5)
  
  ggplot(vcfdf[is.na(vcfdf$bedid) & vcfdf$confidence!="Low",], aes(x=QUAL, fill=confidence)) +
    geom_histogram() +
    scale_x_log10() + 
    labs(title=paste0("Concordence with curated RP call detection - ", sample))
  ggsave(paste0("uncalled_hist_", sample, ".png"), width=10, height=7.5)
  
  ###############
  # Microhomology size distribution
  ###############
  ggplot(vcfdf, aes(x=HOMLEN, color=confidence)) + geom_density(adjust=2) + scale_x_continuous(limits=c(0, 25))
  ggsave(paste0("homlen_", sample, ".png"), width=10, height=7.5)
  
  ###############
  # Assembly rate
  ###############
  ggplot(vcfdf, aes(x=QUAL-ASQ-RASQ, fill=hasAS)) +
    geom_bar(position='fill') +
    scale_x_log10(limits=c(25, 10000), expand=c(0, 0)) +
    scale_y_continuous(labels=percent, expand=c(0, 0)) +
    labs(title="Assembly rate", x="Quality of read pair and split read evidence", y="")
  ggsave(paste0("assembly_rate_qual_", sample, ".png"), width=10, height=7.5)
  
  ###############
  # debugging
  ###############
  return (list(bed=rp, gridss=vcfdf, vcf=vcf))
}




































