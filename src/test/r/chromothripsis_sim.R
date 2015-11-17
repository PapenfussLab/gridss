source("libgridss.R")
#source("libneochromosome.R")
source("libvcf.R")
library(data.table)
library(stringr)
library(ggplot2)
library(scales)
library(parallel)
library(foreach)
library(doParallel) #install.packages("doSNOW")
#cl <- makeCluster(detectCores(), outfile="")
#registerDoParallel(cl)

theme_set(theme_bw())
scale_y_power4 <- scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.1)^4, labels=c("0", "", "", "", "", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"))
scale_y_power5 <- scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.1)^5, labels=c("0", "", "", "", "", "0.5", "", "0.7", "0.8", "0.9", "1.0"))

#########################
# Simulation Comparison Data
########################
pwd <- getwd()
setwd(paste0(ifelse(as.character(Sys.info())[1] == "Windows", "W:/", "~/"), "i/data.fastcompare")) # setwd("W:/i/data.fastcompare")
metadata <- LoadMetadata()
vcfs <- LoadVcfs(metadata)
setwd(pwd)

# Separate out GRIDSS confidence levels
vcfs <- c(
  vcfs,
  lapply(vcfs, function(vcf) {
    if (is.null(attr(vcf, "metadata")) || is.na(vcf@metadata$CX_CALLER) || vcf@metadata$CX_CALLER !="gridss") return(NULL)
    vcf <- vcf[gridss.vcftodf(vcf)$assembly %in% c("Both"),]
    attr(vcf, "metadata")$CX_CALLER <- "gridss Both"
    return(vcf)
  }),
  lapply(vcfs, function(vcf) {
    if (is.null(attr(vcf, "metadata")) || is.na(vcf@metadata$CX_CALLER) || vcf@metadata$CX_CALLER !="gridss") return(NULL)
    vcf <- vcf[gridss.vcftodf(vcf)$assembly %in% c("Both", "Single"),]
    attr(vcf, "metadata")$CX_CALLER <- "gridss Either"
    return(vcf)
  }))
#Separate out GRIDSS confidence levels
# vcfs <- c(
#  vcfs,
#  lapply(vcfs, function(vcf) {
#    if (is.null(attr(vcf, "metadata")) || is.na(vcf@metadata$CX_CALLER) || vcf@metadata$CX_CALLER !="gridss") return(NULL)
#    if (!is.na(vcf@metadata$CX_CALLER_FLAGS) & vcf@metadata$CX_CALLER_FLAGS == "ASSEMBLY_ALGORITHM") {
#      #attr(vcf, "metadata")$CX_CALLER <- "gridss Positional" 
#    } else {
#      return(NULL)
#      #attr(vcf, "metadata")$CX_CALLER <- "gridss Subgraph" 
#    }
#    return(vcf)
#  }))
#   lapply(vcfs, function(vcf) {
#     if (is.null(attr(vcf, "metadata")) || is.na(vcf@metadata$CX_CALLER) || vcf@metadata$CX_CALLER !="gridss") return(NULL)
#     vcf <- vcf[gridss.vcftodf(vcf)$confidence %in% c("High", "Medium"),]
#     attr(vcf, "metadata")$CX_CALLER <- "gridss Medium"
#     return(vcf)
#   }))
vcfs[sapply(vcfs, is.null)] <- NULL
vcfs <- lapply(vcfs, function(vcf) {
  if (all(is.na(rowRanges(vcf)$QUAL))) {
    caller <- vcf@metadata$CX_CALLER
    # use total read support as a qual proxy
    if (caller %in% c("delly")) {
      rowRanges(vcf)$QUAL <- ifelse(is.na(info(vcf)$PE), 0, info(vcf)$PE) + ifelse(is.na(info(vcf)$SR), 0, info(vcf)$SR)
    } else if (caller %in% c("crest")) {
      rowRanges(vcf)$QUAL <- ifelse(is.na(info(vcf)$right_softclipped_read_count), 0, info(vcf)$right_softclipped_read_count) + ifelse(is.na(info(vcf)$left_softclipped_read_count), 0, info(vcf)$left_softclipped_read_count)
    } else if (caller %in% c("pindel")) {
      rowRanges(vcf)$QUAL <- geno(vcf)$AD[,1,2]
    } else {
      warning(paste("No QUAL scores for ", caller))
    }
  }
  return(vcf)
})
vcfs_passfilters <- lapply(vcfs, function(vcf) vcf[rowRanges(vcf)$FILTER %in% c(".", "PASS")])


#TODO: apply filters
#Rprof("C:/dev/Rprof.out")
#truthlist_all <- CalculateTruthSummary(vcfs, maxerrorbp=100, ignore.strand=TRUE)
truthlist_all <- CalculateTruthSummary(vcfs, maxerrorbp=100, ignore.strand=TRUE) # breakdancer does not specify strand
truthlist_passfilters <- CalculateTruthSummary(vcfs_passfilters, maxerrorbp=100, ignore.strand=TRUE) # breakdancer does not specify strand

gridssfirst <- function(x) {
  return (relevel(relevel(relevel(x, "gridss"), "gridss Either"), "gridss Both"))
}
truthlist_all$calls$CX_CALLER <- gridssfirst(truthlist_all$calls$CX_CALLER)
truthlist_all$truth$CX_CALLER <- gridssfirst(truthlist_all$truth$CX_CALLER)
truthlist_passfilters$calls$CX_CALLER <- gridssfirst(truthlist_passfilters$calls$CX_CALLER)
truthlist_passfilters$truth$CX_CALLER <- gridssfirst(truthlist_passfilters$truth$CX_CALLER)

# Sensitivity
dtsenssize_all <- truthlist_all$truth[, list(sens=sum(tp)/.N), by=c("SVLEN", "SVTYPE", "CX_ALIGNER", "CX_ALIGNER_SOFTCLIP", "CX_CALLER", "CX_READ_DEPTH", "CX_READ_FRAGMENT_LENGTH", "CX_READ_LENGTH", "CX_REFERENCE_VCF_VARIANTS")]
dtsenssize_passfilters <- truthlist_passfilters$truth[, list(sens=sum(tp)/.N), by=c("SVLEN", "SVTYPE", "CX_ALIGNER", "CX_ALIGNER_SOFTCLIP", "CX_CALLER", "CX_READ_DEPTH", "CX_READ_FRAGMENT_LENGTH", "CX_READ_LENGTH", "CX_REFERENCE_VCF_VARIANTS")]
dtsenssize_all$Filter <- "Unfiltered"
dtsenssize_passfilters$Filter <- "Pass"
dtsenssize <- rbind(dtsenssize_all, dtsenssize_passfilters)
# Precision
dtprec_all <- truthlist_all$truth[, list(prec=sum(tp)/.N), by=c("CX_ALIGNER", "CX_ALIGNER_SOFTCLIP", "CX_CALLER", "CX_READ_DEPTH", "CX_READ_FRAGMENT_LENGTH", "CX_READ_LENGTH", "CX_REFERENCE_VCF_VARIANTS")]
dtprec_passfilters <- truthlist_passfilters$truth[, list(prec=sum(tp)/.N), by=c("CX_ALIGNER", "CX_ALIGNER_SOFTCLIP", "CX_CALLER", "CX_READ_DEPTH", "CX_READ_FRAGMENT_LENGTH", "CX_READ_LENGTH", "CX_REFERENCE_VCF_VARIANTS")]
dtprec_all$Filter <- "Unfiltered"
dtprec_passfilters$Filter <- "Pass"
dtprec <- rbind(dtprec_all, dtprec_passfilters)

plot_sens <- ggplot(dtsenssize[!(dtsenssize$CX_REFERENCE_VCF_VARIANTS %in% c("homBP", "homBP_SINE"))]) + # & dtsenssize$CX_CALLER %in% c("gridss Positional", "gridss Subgraph")
  aes(y=sens, x=SVLEN, shape=factor(CX_READ_FRAGMENT_LENGTH), color=Filter) +
  geom_line() + 
  facet_grid(CX_ALIGNER + CX_CALLER ~ CX_READ_FRAGMENT_LENGTH + CX_READ_LENGTH + CX_READ_DEPTH + CX_REFERENCE_VCF_VARIANTS) +
  scale_x_log10(breaks=2**(0:16)) +  
  theme(panel.grid.major = element_blank()) +
  labs(y="sensitivity", x="Event size", title="Sensitivity")
ggsave("sens_size.png", plot_sens, width=10, height=7.5)
# rescaled to show close to 1
ggsave("sens_size_nonlinearscale.png", plot_sens + aes(y=sens**5) + scale_y_continuous(breaks=seq(0,1,0.1)**5, labels=seq(0,1,0.1)), width=10, height=7.5)

ggplot(truthlist_all$truth) + aes(x=errorsize) + facet_grid(CX_REFERENCE_VCF_VARIANTS~CX_CALLER) + geom_histogram() + scale_y_log10()
ggplot(truthlist_all$calls) + aes(x=partialtp) + facet_grid(CX_REFERENCE_VCF_VARIANTS~CX_CALLER) + geom_histogram()

for (rl in unique(dtsenssize$CX_READ_LENGTH)) {
for (rd in unique(dtsenssize$CX_READ_DEPTH)) {
for (fragsize in unique(dtsenssize$CX_READ_FRAGMENT_LENGTH)) {
  dt <- dtsenssize
  dt <- dt[dt$CX_READ_DEPTH==rd & dt$CX_READ_FRAGMENT_LENGTH==fragsize & dt$CX_READ_LENGTH==rl,]
  dt[dt$CX_REFERENCE_VCF_VARIANTS=="hetDEL",]$CX_REFERENCE_VCF_VARIANTS <- "Deletion"
  dt[dt$CX_REFERENCE_VCF_VARIANTS=="hetINS",]$CX_REFERENCE_VCF_VARIANTS <- "Insertion"
  dt[dt$CX_REFERENCE_VCF_VARIANTS=="hetDUP",]$CX_REFERENCE_VCF_VARIANTS <- "Tandem Duplication"
  dt[dt$CX_REFERENCE_VCF_VARIANTS=="hetINV",]$CX_REFERENCE_VCF_VARIANTS <- "Inversion"
  ggplot(dt[!is.na(dt$SVLEN),]) +
    aes(y=sens, x=SVLEN, color=CX_CALLER, shape=CX_CALLER, linetype=Filter) +
    geom_line() + 
    geom_point() + 
    facet_wrap( ~ CX_REFERENCE_VCF_VARIANTS) +
    #facet_grid(CX_CALLER ~ CX_REFERENCE_VCF_VARIANTS) +
    scale_x_log10(breaks=2**(0:16), labels=c("1", "2", "4", "8", "16", "32", "64", "128", "256", "512", "1k", "2k", "4k", "8k", "16k", "32k", "64k")) +
    aes(y=sens**4) + scale_y_power4 + 
    scale_color_brewer(type="qual", palette=2) +
    theme(panel.grid.major = element_blank()) +
    labs(y="Sensitivity", x="Event size", title=paste0("Sensitivity ", rd, "x coverage"))
  ggsave(file=paste0("sens_simple_rl", rl, "_", fragsize, "_", rd, "x.pdf"), width=300, height=150, units=c("mm"))
}}}
for (rl in unique(dtsenssize$CX_READ_LENGTH)) {
for (rd in unique(dtsenssize$CX_READ_DEPTH)) {
for (fragsize in unique(dtsenssize$CX_READ_FRAGMENT_LENGTH)) {
  # bar graph for precisions
  dt <- dtprec_passfilters
  dt <- dt[dt$CX_READ_DEPTH==rd & dt$CX_READ_FRAGMENT_LENGTH==fragsize & dt$CX_READ_LENGTH==rl,]
  dt <- dt[dt$CX_REFERENCE_VCF_VARIANTS %in% c("hetDEL", "hetINS", "hetDUP", "hetINV")]
  dt[dt$CX_REFERENCE_VCF_VARIANTS=="hetDEL",]$CX_REFERENCE_VCF_VARIANTS <- "Deletion"
  dt[dt$CX_REFERENCE_VCF_VARIANTS=="hetINS",]$CX_REFERENCE_VCF_VARIANTS <- "Insertion"
  dt[dt$CX_REFERENCE_VCF_VARIANTS=="hetDUP",]$CX_REFERENCE_VCF_VARIANTS <- "Tandem Duplication"
  dt[dt$CX_REFERENCE_VCF_VARIANTS=="hetINV",]$CX_REFERENCE_VCF_VARIANTS <- "Inversion"
  ggplot(dt) +
    aes(y=prec, x=CX_CALLER, fill=CX_CALLER) +
    geom_bar(stat="identity") + 
    facet_wrap( ~ CX_REFERENCE_VCF_VARIANTS) +
    scale_color_brewer(type="qual", palette=2) +
    labs(y="Precision", title=paste0("Precision ", rd, "x coverage"))
  ggsave(file=paste0("prec_simple_rl", rl, "_", fragsize, "_", rd, "x.pdf"), width=300, height=150, units=c("mm"))
}}}

############
# Error in called position
#
ggplot(truthlist_all$truth[!is.na(truthlist_all$truth$errorsize),]) +
  aes(x=errorsize, color=factor(SVLEN)) +
  geom_density(aes(y=..density..)) +
  facet_grid(CX_ALIGNER + CX_CALLER ~ CX_READ_FRAGMENT_LENGTH + CX_READ_LENGTH + CX_READ_DEPTH, scale="free")


############
# ROC
#
bylist <- c("CX_ALIGNER", "CX_ALIGNER_SOFTCLIP", "CX_CALLER", "CX_READ_DEPTH", "CX_READ_FRAGMENT_LENGTH", "CX_READ_LENGTH", "CX_REFERENCE_VCF_VARIANTS")
dtrocunfiltered <- TruthSummaryToROC(truthlist_all, bylist=bylist)
dtrocunfiltered$Filter <- "Unfiltered"
dtrocfiltered <- TruthSummaryToROC(truthlist_passfilters, bylist=bylist)
dtrocfiltered$Filter <- "Pass"
dtrocall <- rbind(dtrocunfiltered, dtrocfiltered)

plot_roc <- ggplot(dtrocall[dtrocall$CX_REFERENCE_VCF_VARIANTS %in% c("homBP", "homBP_SINE"),]) + 
  aes(y=sens, x=fp+1, color=log(QUAL), shape=Filter) +
  scale_colour_gradientn(colours = rainbow(11)) +
  geom_line() + 
  geom_point(size=0.25) + 
  facet_grid(CX_ALIGNER + CX_CALLER ~ CX_REFERENCE_VCF_VARIANTS + CX_READ_FRAGMENT_LENGTH + CX_READ_LENGTH + CX_READ_DEPTH) +
  scale_x_log10() + 
  labs(y="sensitivity", title="ROC by call quality threshold")
ggsave("roc_full.png", plot_roc, width=10, height=7.5)
ggsave("roc_full_nonlinear_scale.png", plot_roc + scale_x_continuous(limits=c(0, 100)) + aes(x=fp) + aes(y=sens**5) + scale_y_continuous(breaks=seq(0,1,0.1)**5, labels=seq(0,1,0.1)), width=10, height=7.5)

for (variant in unique(dtrocall$CX_REFERENCE_VCF_VARIANTS)) {
  dtroc <- dtrocall[dtrocall$CX_REFERENCE_VCF_VARIANTS == variant,]
  plot_roc_all <- ggplot(dtroc) +
    aes(y=sens, shape=factor(CX_READ_FRAGMENT_LENGTH), color=factor(CX_READ_DEPTH), linetype=Filter) +
    geom_point(size=1) + 
    facet_grid(CX_ALIGNER + CX_CALLER ~ CX_REFERENCE_VCF_VARIANTS + CX_READ_LENGTH) +
    labs(y="sensitivity", title="ROC by call quality threshold")
  plot_roc_all + aes(x=prec) + labs(x="precision") + scale_x_reverse()
  ggsave(paste0("roc_full_precision_", variant, ".png"), width=10, height=7.5)
  plot_roc_all + aes(x=fdr) + labs(x="false discovery rate")
  ggsave(paste0("roc_full_fdr_", variant, ".png"), width=10, height=7.5)
  plot_roc_all + aes(x=fp) + labs(x="total false positives")
  ggsave(paste0("roc_full_fp_", variant, ".png"), width=10, height=7.5)
  
  ggplot(dtroc[dtroc$CX_READ_FRAGMENT_LENGTH==300 & dtroc$CX_READ_LENGTH==100,]) +
    aes(x=fdr, y=sens, color=log10(QUAL), shape=CX_ALIGNER) +
    geom_point() +
    facet_grid(CX_CALLER + CX_REFERENCE_VCF_VARIANTS + Filter ~ CX_READ_DEPTH) + 
    scale_colour_gradientn(colours = rainbow(11)) +
    labs(x="false discovery rate", y="sensitivity", title="ROC curve by read depth 2x100bp 300bp fragment")
  ggsave(paste0("roc_depth_", variant, ".png"), width=10, height=7.5)
  
  ggplot(dtroc[dtroc$CX_READ_FRAGMENT_LENGTH==300 & dtroc$CX_READ_DEPTH==100,]) +
    aes(x=fdr, y=sens, color=log10(QUAL)) +
    geom_point() +
    facet_grid(CX_ALIGNER + Filter ~ CX_READ_LENGTH) + 
    scale_colour_gradientn(colours = rainbow(11)) +
    labs(x="false discovery rate", y="sensitivity", title="ROC curve by read length 100x 300bp fragment")
  ggsave(paste0("roc_length_", variant, ".png"), width=10, height=7.5)
  
  ggplot(dtroc[dtroc$CX_READ_DEPTH==100 & dtroc$CX_READ_LENGTH==100,]) +
    aes(x=fdr, y=sens, color=log10(QUAL)) +
    geom_point() +
    facet_grid(CX_ALIGNER + Filter ~ CX_READ_FRAGMENT_LENGTH) + 
    scale_colour_gradientn(colours = rainbow(11)) +
    labs(x="false discovery rate", y="sensitivity", title="ROC curve by fragment size 100x 2x100bp")
  ggsave(paste0("roc_fragsize_", variant, ".png"), width=10, height=7.5)
}
for (rl in unique(dtrocall$CX_READ_LENGTH)) {
for (rd in unique(dtrocall$CX_READ_DEPTH)) {
for (fragsize in unique(dtrocall$CX_READ_FRAGMENT_LENGTH)) {
  dt <- dtrocall
  dt <- dt[dt$CX_READ_DEPTH==rd & dt$CX_READ_FRAGMENT_LENGTH==fragsize & dt$CX_READ_LENGTH==rl,]
  dt <- dt[dt$CX_REFERENCE_VCF_VARIANTS %in% c("homBP","homBP_SINE"),]
  dt[dt$CX_REFERENCE_VCF_VARIANTS=="homBP",]$CX_REFERENCE_VCF_VARIANTS <- "Breakpoint"
  dt[dt$CX_REFERENCE_VCF_VARIANTS=="homBP_SINE",]$CX_REFERENCE_VCF_VARIANTS <- "Breakpoint at SINE/ALU"
  ggplot(dt) + 
    aes(y=sens, x=fp+1, color=CX_CALLER, linetype=Filter) +
    scale_color_brewer(type="qual", palette=2) +
    facet_wrap(~ CX_REFERENCE_VCF_VARIANTS, ncol=1) + 
    geom_line() + 
    geom_point() + 
    scale_x_log10() + 
    aes(y=sens**4) + scale_y_power4 + 
    labs(y="Sensitivity", x="False Positives", title="ROC by call quality threshold")
  ggsave(file=paste0("roc_rl", rl, "_", fragsize, "_", rd, "x.pdf"), width=200, height=150, units=c("mm"))
}}}















