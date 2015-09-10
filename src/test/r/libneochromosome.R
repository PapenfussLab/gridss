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
  # We'll go with the bedid if it exists, then fall back to the gridssid
  # TODO: halve all Venn Diagram counts so we're counting breakpoints, not breakends
  venn.plot <- venn.diagram(
    x = list(
      Published=names(rp),
      Spanning=names(rp[rp$spanning,]),
      Both=vreplaceNA(vcfdf[vcfdf$assembly=="Both",]$bedid, rownames(vcfdf[vcfdf$assembly=="Both",])),
      Single=vreplaceNA(vcfdf[vcfdf$assembly=="Single",]$bedid, rownames(vcfdf[vcfdf$assembly=="Single",])),
      NoAssembly=vreplaceNA(vcfdf[vcfdf$assembly=="No assembly",]$bedid, rownames(vcfdf[vcfdf$assembly=="No assembly",]))),
    filename = paste0("neo_", "venn_", sample, ".png"))
    
  if (!file.exists(paste0("neo_", sample, ".additional.both.vcf"))) {
    # don't rewrite until writeVcf trailing tab bug is fixed
    writeVcf(vcf[is.na(vcfdf$bedid) & vcfdf$assembly=="Both",], paste0("neo_", sample, ".additional.both.vcf"))
    writeVcf(vcf[is.na(vcfdf$bedid) & vcfdf$assembly=="Single",], paste0("neo_", sample, ".additional.single.vcf"))
  }
  write.csv(vcfdf[is.na(vcfdf$bedid) & vcfdf$assembly %in% c("Single", "Both"),], paste0("neo_", sample, ".additional.csv"))
  export(rp[rp$assembly %in% c("None", "No assembly"),], paste0("neo_", sample, ".misses.bed"))
  write.csv(as.data.frame(rp[rp$assembly %in% c("None", "No assembly"),]), paste0("neo_", sample, ".missed.csv"))
  
  
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
#   ggsave(paste0("neo_", "rp_roc_", sample, ".png"), width=10, height=7.5)
  
  ggplot(as.data.frame(mcols(rp)), aes(x=nreads, fill=hits)) +
    geom_histogram() +
    scale_x_log10() + 
    labs(title=paste0("Sensitivity of curated RP call detection - ", sample))
  ggsave(paste0("neo_rp_hist_", sample, ".png"), width=10, height=7.5)
  
  ggplot(vcfdf[is.na(vcfdf$bedid) & vcfdf$assembly!="No assembly",], aes(x=QUAL, fill=assembly)) +
    geom_histogram() +
    scale_x_log10() + 
    labs(title=paste0("Concordence with curated RP call detection - ", sample))
  ggsave(paste0("neo_uncalled_hist_", sample, ".png"), width=10, height=7.5)
  
  ###############
  # Microhomology size distribution
  ###############
  ggplot(vcfdf, aes(x=HOMLEN, color=assembly)) + geom_density(adjust=2) + scale_x_continuous(limits=c(0, 25))
  ggsave(paste0("neo_homlen_", sample, ".png"), width=10, height=7.5)
  
  ###############
  # Assembly rate
  ###############
  ggplot(vcfdf, aes(x=QUAL-ASQ-RASQ, fill=assembly)) +
    geom_bar(position='fill') +
    scale_x_log10(limits=c(25, 10000), expand=c(0, 0)) +
    scale_y_continuous(labels=percent, expand=c(0, 0)) +
    labs(title="Assembly rate", x="Quality of read pair and split read evidence", y="")
  ggsave(paste0("neo_assembly_rate_qual_", sample, ".png"), width=10, height=7.5)

  ###############
  # Coverage of gridss-only both assembly
  ###############
  assnocall <- vcfdf[is.na(vcfdf$bedid) & vcfdf$assembly=="Both",]
  ggplot(assnocall) +
    aes(x=RP) + 
    geom_histogram() +
    labs("Read Pair coverage, mutual breakend assembly, no published call")
  ggsave(paste0("neo_both_nocall_rp_", sample, ".png"), width=10, height=7.5)
  ggplot(assnocall) +
    aes(x=asWithin1kbp-1) +
    geom_histogram() +
    labs("Breakends with 1kbp, mutual breakend assembly, no published call")
  ggsave(paste0("neo_both_nocall_calldensity_", sample, ".png"), width=10, height=7.5)
  
  ###############
  # debugging
  ###############
  return (list(bed=rp, gridss=vcfdf, vcf=vcf))
}






# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}





























