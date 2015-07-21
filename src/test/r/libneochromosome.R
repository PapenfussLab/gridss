library(openxlsx) # install.packages("openxlsx")
library(ggplot2)
library(scales)
library(RColorBrewer)
library(stringr)

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

go <- function(sample, vcf, rp) {
  vcf <- gridss.removeUnpartnerededBreakend(vcf)
  rpmaxgap=200 # 778 #chr1:188377591-chr1:188379026- call location is over 120
  rpMate <- rp[rp$mate,]
  matches <- gridss.annotateBreakpoints(rp, rpMate, vcf, maxgap=rpmaxgap)
  vcfdf <- matches$gridss
  rp <- matches$bed
  
  rp$hits <- vcfdf[rp$gridssid,]$confidence
  levels(rp$hits) <- c("Low", "Medium", "High", "None")
  rp$hits[is.na(rp$hits)] <- "None"
  
  
  ###############
  # Read Pair variant calling concordance
  ###############
  #library(pROC) #install.packages("pROC")
  dfroc <- NULL
  for (confidence in as.numeric(unique(vcfdf$confidence))) {
    subset <- as.numeric(vcfdf$confidence) >= confidence
    dftmp <- vcftoroc(vcf[subset], rp, maxgap=rpmaxgap)
    dftmp$confidence <- levels(vcfdf$confidence)[confidence]
    dfroc <- rbind(dfroc, dftmp)  
  }
  ggplot(dfroc[dfroc$qual!=0,], aes(x=log(qual+1), y=sens, color=confidence)) +
    geom_line() +
    scale_x_reverse() +
    theme_bw() +
    scale_y_continuous(limits=c(0,1)) +
    labs(title=paste0("Sensitivity of curated RP call detection - ", sample))
  ggsave(paste0("rp_roc_", sample, ".png"), width=10, height=7.5)
  
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
  return (list(bed=rp, gridss=vcfdf))
}




































