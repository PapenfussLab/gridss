library(openxlsx) # install.packages("openxlsx")
library(VennDiagram) # install.packages("VennDiagram")
library(ggplot2)
library(scales)
library(RColorBrewer)
library(stringr)
library(rtracklayer)
source("libgridss.R")
source("libvcf.R")
source("common.R")

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
go <- function(sample, vcf, rp, filterBed=NULL, minimumEventSize, cgrmaxgap=1000, matchmaxgap=200, graphs=TRUE, minqual=0) {
  if (!is.null(filterBed)) {
    # filter so that both ends are within the filterBed
    rp <- subsetByOverlaps(rp, filterBed, maxgap=cgrmaxgap)
    rp <- rp[rp$mate %in% names(rp)]
    vcf <- vcf[overlapsAny(rowRanges(vcf), filterBed, maxgap=cgrmaxgap),]
    vcf <- vcf[rowRanges(vcf)$QUAL >= minqual,]
  }
  vcf <- gridss.removeUnpartnerededBreakend(vcf)
  rpMate <- rp[rp$mate,]
  matches <- gridss.annotateBreakpoints(rp, rpMate, vcf, maxgap=matchmaxgap)
  vcfdf <- matches$gridss
  rp <- matches$bed
  
  rp$hits <- vcfdf[rp$gridssid,]$assembly
  levels(rp$hits) <- c("Low", "Medium", "High", "None")
  rp$hits[is.na(rp$hits)] <- "None"
  
  rp$assembly <- vcfdf[rp$gridssid,]$assembly
  levels(rp$assembly) <- c("No assembly", "Single", "Both", "None")
  rp$assembly[is.na(rp$assembly)] <- "None"
  
  # filter out small events that don't hit an existing call
  vcf <- vcf[is.na(vcfdf$size) | vcfdf$size >= minimumEventSize,]
  vcfdf <- vcfdf[rownames(vcf),]

  rp$spanning <- FALSE
  spanning <- FindFragmentSpanningEvents(rp, vcftobpgr(vcf))
  rp[spanning$query,]$spanning <- TRUE
  rp[rp[spanning$query,]$mate,]$spanning <- TRUE
  
  rp$distanceToas <- distanceToClosest(rp, rowRanges(vcf[vcfdf$assembly != "No assembly"]))
  vcfdf$distanceToCall <- distanceToClosest(rowRanges(vcf), rp)
  vcfdf$found <- !is.na(vcfdf$bedid)
  
  rp$asWithin1kbp <- countOverlaps(rp, rowRanges(vcf[vcfdf$assembly %in% c("Single", "Both")]), maxgap=1000)
  vcfdf$asWithin1kbp <- countOverlaps(vcf, rowRanges(vcf[vcfdf$assembly %in% c("Single", "Both")]), maxgap=1000)
  
  ###############
  # Variant calling concordance
  ###############
  # VennDiagram requires lists of ids
  # since 
  # We'll go with the bedid if it exists, then fall back to the gridssid
  vennlist=list(
    Published=names(rp),
    Spanning=names(rp[rp$spanning,]),
    Both=unique(c(names(rp)[rp$assembly=="Both"], vreplaceNA(vcfdf[vcfdf$assembly=="Both",]$bedid, rownames(vcfdf[vcfdf$assembly=="Both",])))),
    Single=unique(c(names(rp)[rp$assembly=="Single"], vreplaceNA(vcfdf[vcfdf$assembly=="Single",]$bedid, rownames(vcfdf[vcfdf$assembly=="Single",])))),
    NoAssembly=unique(c(names(rp)[rp$assembly=="No assembly"], vreplaceNA(vcfdf[vcfdf$assembly=="No assembly",]$bedid, rownames(vcfdf[vcfdf$assembly=="No assembly",])))))
  # halve all Venn Diagram counts so we're counting breakpoints, not breakends
  vennlist$Published <- vennlist$Published[str_detect(vennlist$Published, "((/1)|(o))$")]
  vennlist$Spanning <- vennlist$Spanning[str_detect(vennlist$Spanning, "((/1)|(o))$")]
  vennlist$Both <- vennlist$Both[str_detect(vennlist$Both, "((/1)|(o))$")]
  vennlist$Single <- vennlist$Single[str_detect(vennlist$Single, "((/1)|(o))$")]
  vennlist$NoAssembly <- vennlist$NoAssembly[str_detect(vennlist$NoAssembly, "((/1)|(o))$")]
  if (graphs) {
    venn.plot <- venn.diagram(x=vennlist, filename = paste0("neo_maxgap", matchmaxgap, "_", "venn_", sample, ""))
    
    if (!file.exists(paste0("neo_maxgap", matchmaxgap, "_",sample, ".additional.both.vcf"))) {
      # don't rewrite until writeVcf trailing tab bug is fixed
      writeVcf(vcf[is.na(vcfdf$bedid) & vcfdf$assembly=="Both",], paste0("neo_maxgap", matchmaxgap, "_",sample, ".additional.both.vcf"))
      writeVcf(vcf[is.na(vcfdf$bedid) & vcfdf$assembly=="Single",], paste0("neo_maxgap", matchmaxgap, "_",sample, ".additional.single.vcf"))
    }
    write.csv(vcfdf[is.na(vcfdf$bedid) & vcfdf$assembly %in% c("Single", "Both"),], paste0("neo_maxgap", matchmaxgap, "_",sample, ".additional.csv"))
    export(rp[rp$assembly %in% c("None", "No assembly"),], paste0("neo_maxgap", matchmaxgap, "_",sample, ".misses.bed"))
    write.csv(as.data.frame(rp[rp$assembly %in% c("None", "No assembly"),]), paste0("neo_maxgap", matchmaxgap, "_",sample, ".missed.csv"))
  
  
  ###############
  # Read Pair variant calling concordance
  ###############
  #library(pROC) #install.packages("pROC")
#   dfroc <- NULL
#   for (assembly in as.numeric(unique(vcfdf$assembly))) {
#     subset <- as.numeric(vcfdf$assembly) >= assembly
#     dftmp <- vcftoroc(vcf[subset], rp, maxgap=rpmaxgap)
#     dftmp$assembly <- levels(vcfdf$assembly)[assembly]
#     dfroc <- rbind(dfroc, dftmp)  
#   }
#   ggplot(dfroc[dfroc$qual!=0,], aes(x=log(qual+1), y=sens, color=assembly)) +
#     geom_line() +
#     scale_x_reverse() +
#     theme_bw() +
#     scale_y_continuous(limits=c(0,1)) +
#     labs(title=paste0("Sensitivity of curated RP call detection - ", sample))
#   saveplot(paste0("neo_maxgap", matchmaxgap, "_","rp_roc_", sample, ""), width=10, height=7.5)
  
    ggplot(as.data.frame(mcols(rp)), aes(x=nreads, fill=hits)) +
      geom_histogram() +
      scale_x_log10() + 
      labs(title=paste0("Sensitivity of curated RP call detection - ", sample))
    saveplot(paste0("neo_rp_hist_", sample, ""), width=10, height=7.5)
    
    ggplot(vcfdf[is.na(vcfdf$bedid) & vcfdf$assembly!="No assembly",], aes(x=QUAL, fill=assembly)) +
      geom_histogram() +
      scale_x_log10() + 
      labs(title=paste0("Concordence with curated RP call detection - ", sample))
    saveplot(paste0("neo_uncalled_hist_", sample, ""), width=10, height=7.5)
    
    ###############
    # Microhomology size distribution
    ###############
    ggplot(vcfdf, aes(x=HOMLEN, color=assembly)) + geom_density(adjust=2) + scale_x_continuous(limits=c(0, 25))
    saveplot(paste0("neo_homlen_", sample, ""), width=10, height=7.5)
    
    ###############
    # Assembly rate
    ###############
    ggplot(vcfdf, aes(x=QUAL-ASQ-RASQ, fill=assembly)) +
      geom_bar(position='fill') +
      scale_x_log10(limits=c(25, 10000), expand=c(0, 0)) +
      scale_y_continuous(labels=percent, expand=c(0, 0)) +
      labs(title="Assembly rate", x="Quality of read pair and split read evidence", y="")
    saveplot(paste0("neo_assembly_rate_qual_", sample, ""), width=10, height=7.5)
  
    ###############
    # Coverage of gridss-only both assembly
    ###############
    assnocall <- vcfdf[is.na(vcfdf$bedid) & vcfdf$assembly=="Both",]
    ggplot(assnocall) +
      aes(x=RP) + 
      geom_histogram() +
      labs("Read Pair coverage, mutual breakend assembly, no published call")
    saveplot(paste0("neo_both_nocall_rp_", sample, ""), width=10, height=7.5)
    ggplot(assnocall) +
      aes(x=asWithin1kbp-1) +
      geom_histogram() +
      labs("Breakends with 1kbp, mutual breakend assembly, no published call")
    saveplot(paste0("neo_both_nocall_calldensity_", sample, ""), width=10, height=7.5)
  }

  ###############
  # prevalence of microhomology and untemplated sequence
  ###############
  dfboth <- vcfdf[vcfdf$assembly=="Both",]
  insmicro <- data.frame(sample=sample, matchmaxgap=matchmaxgap, total=nrow(dfboth), hom=sum(dfboth$HOMLEN>0), ins=sum(nchar(dfboth$INSSEQ)>0))

  ###############
  # 
  ###############
  #dtsummarised <- rbind(
  #  data.table(rp)[, list(count=.N, hit=TRUE), by=list(sample, matchmaxgap, assembly, spanning)],
  #  data.table(vcfdf[is.na(vcfdf$bedid),])[, list(count=.N, hit=FALSE, spanning=NA), by=list(sample, matchmaxgap, assembly)])
  #dtsummarised$count <- dtsummarised$count / 2 # breakend -> breakpoint conversion
  #ggplot()


  ###############
  # debugging
  ###############
  rp$sample <- sample
  rp$matchmaxgap <- matchmaxgap
  vcfdf$sample <- sample
  vcfdf$matchmaxgap <- matchmaxgap
  return (list(bed=rp, gridss=vcfdf, vcf=vcf, insmicro=insmicro))
}























