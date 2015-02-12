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
# lists overlapping breakpoints
breakpointHits <- function(queryGr, subjectGr, mateQueryGr=queryGr[queryGr$mate,], mateSubjectGr=subjectGr[subjectGr$mate,], ...) {
  dfhits <- rbind(as.data.frame(findOverlaps(queryGr, subjectGr, ...)),
                  as.data.frame(findOverlaps(mateQueryGr, mateSubjectGr, ...)))
  dfhits <- dfhits[duplicated(dfhits),] # both breakends match
  return(dfhits)
}
# counts overlapping breakpoints
countBreakpointHits <- function(queryGr, subjectGr, mateQueryGr=queryGr[queryGr$mate,], mateSubjectGr=subjectGr[subjectGr$mate,], ...) {
  dfhits <- breakpointHits(queryGr, subjectGr, mateQueryGr, mateSubjectGr, ...)
  queryHits <- rep(0, length(queryGr))
  queryHits[count(dfhits, "queryHits")$queryHits] <- count(dfhits, "queryHits")$freq
  subjectHits <- rep(0, length(subjectGr))
  subjectHits[count(dfhits, "subjectHits")$subjectHits] <- count(dfhits, "subjectHits")$freq
  return(list(queryHitCount=queryHits, subjectHitCount=subjectHits))
}
# counts overlapping breakpoints
# VCF records MUST each be a valid symbolic breakend alleles
# VCF MUST have a single valid MATEID for every record
# GRanges MUST have a $mate set to the name of the other breakend for that breakpoint
countVcfGrBreakpointHits <- function(vcf, gr, ...) {
  vcfgr <- vcftobpgr(vcf)
  return(countBreakpointHits(vcfgr, gr, vcfgr[vcfgr$mate,], gr[gr$mate,], ...))
}
# converts a breakpoint VCF to a breakpoint GR
vcftobpgr <- function(vcf) {
  vcfgr <- rowData(vcf)
  # set strand
  strand(vcfgr) <- ifelse(str_detect(as.character(fixed(vcf)$ALT), "[[:alpha:]]+(\\[|]).*(\\[|])"), "+", "-")
  vcfgr$mate <- as.character(info(vcf)$MATEID)
  return(vcfgr)
}
vcftoroc <- function(vcf, grtp, ...) {
  hits <- breakpointHits(vcftobpgr(vcf), grtp, ...)
  hits$qual <- fixed(vcf)$QUAL[hits$queryHits]
  hits <- hits[order(hits$qual),]
  tp <- rep(0, length(grtp))
  tp[hits$subjectHits] <- hits$qual
  tp <- tp[order(-tp)]
  df <- data.frame(qual=tp, sens=1:length(tp) / length(tp))
  df <- df[order(-df$sens),]
  df <- df[!duplicated(df$qual),]
  return(df)
}

go <- function(sample, vcf, rp) {
  vcf <- gridss.removeUnpartnerededBreakend(vcf)
  vcfdf <- gridss.truthdetails.processvcf.vcftodf(vcf)
  rpmaxgap=120
  hits <- countVcfGrBreakpointHits(vcf, rp, maxgap=rpmaxgap)
  vcfdf$rpHits <- hits$queryHitCount
  rp$hits_all_gap120 <- hits$subjectHitCount
  rp$hits_hc_gap120 <- countVcfGrBreakpointHits(vcf[vcfdf$confidence=="High",], rp, maxgap=rpmaxgap)$subjectHitCount
  rp$hits_mc_gap120 <- countVcfGrBreakpointHits(vcf[vcfdf$confidence=="High" | vcfdf$confidence=="Medium",], rp, maxgap=rpmaxgap)$subjectHitCount
  rp$hits_gap120 <- as.factor((rp$hits_all_gap120 > 0) + (rp$hits_mc_gap120 > 0) + (rp$hits_hc_gap120 > 0))
  levels(rp$hits_gap120) <- c("None", "Low", "Medium", "High")
  
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
  ggsave(paste0(sample, "_rp_roc.png"), width=10, height=7.5)
  
  ggplot(as.data.frame(mcols(rp)), aes(x=nreads, fill=hits_gap120)) +
    geom_histogram() +
    scale_x_log10() + 
    labs(title=paste0("Sensitivity of curated RP call detection - ", sample))
  ggsave(paste0(sample, "_rp_hist.png"), width=10, height=7.5)
  
}




































