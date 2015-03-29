source("../../../../../dev/ws/indelappraisal/libindelappraisal.R")
source("../../main/r/libgridss.R")
source("libneochromosome.R")
source("libvcf.R")
library(data.table)
library(stringr)
library(ggplot2)
library(scales)
library(parallel)
library(foreach)
library(doSNOW) #install.packages("doSNOW")
#cl <- makeCluster(detectCores())
#registerDoSNOW(cl)  

theme_set(theme_bw())

# Load simulation data
pwd <- getwd()
setwd("W:/i/data.fastcompare")
metadata <- LoadMetadata()
vcfs <- LoadVcfs(metadata)
setwd(pwd)

#TODO: apply filters
#Rprof("C:/dev/Rprof.out")
truthlist <- CalculateTruthSummary(vcfs, maxerrorbp=100, ignore.strand=TRUE)
############
# Sensitivity
#
dtsenssize <- truthlist$truth[, list(sens=sum(tp)/.N), by=c("SVLEN", "SVTYPE", "CX_ALIGNER", "CX_ALIGNER_SOFTCLIP", "CX_CALLER", "CX_READ_DEPTH", "CX_READ_FRAGMENT_LENGTH", "CX_READ_LENGTH", "CX_REFERENCE_VCF_VARIANTS")]
ggplot(dtsenssize) +
  aes(y=sens, x=SVLEN, shape=factor(CX_READ_FRAGMENT_LENGTH), color=CX_REFERENCE_VCF_VARIANTS) +
  geom_line() + 
  facet_grid(CX_ALIGNER + CX_CALLER ~ CX_READ_LENGTH + CX_READ_DEPTH) +
  scale_x_log10(breaks=unique(dtsenssize$SVLEN)) +
  labs(y="sensitivity", x="Event size", title="Sensitivity")

ggplot(truthlist$truth) + aes(x=sizeerror) + facet_grid(CX_REFERENCE_VCF_VARIANTS~CX_CALLER) + geom_histogram() + scale_y_log10()
############
# ROC
#
bylist <- c("CX_ALIGNER", "CX_ALIGNER_SOFTCLIP", "CX_CALLER", "CX_READ_DEPTH", "CX_READ_FRAGMENT_LENGTH", "CX_READ_LENGTH", "CX_REFERENCE_VCF_VARIANTS")
dtrocall <- TruthSummaryToROC(truthlist, bylist=bylist)
for (variant in unique(dtrocall$CX_REFERENCE_VCF_VARIANTS)) {
  dtroc <- dtrocall[dtrocall$CX_REFERENCE_VCF_VARIANTS == variant,]
  plot_roc_all <- ggplot(dtroc) +
    aes(y=sens, shape=factor(CX_READ_FRAGMENT_LENGTH), color=factor(CX_READ_DEPTH)) +
    geom_point(size=1) + 
    facet_grid(CX_ALIGNER + CX_CALLER ~ CX_READ_LENGTH + CX_REFERENCE_VCF_VARIANTS) +
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
    facet_grid(CX_CALLER + CX_REFERENCE_VCF_VARIANTS ~ CX_READ_DEPTH) + 
    scale_colour_gradientn(colours = rainbow(11)) +
    labs(x="false discovery rate", y="sensitivity", title="ROC curve by read depth 2x100bp 300bp fragment")
  ggsave(paste0("roc_depth_", variant, ".png"), width=10, height=7.5)
  
  ggplot(dtroc[dtroc$CX_READ_FRAGMENT_LENGTH==300 & dtroc$CX_READ_DEPTH==100,]) +
    aes(x=fdr, y=sens, color=log10(QUAL)) +
    geom_point() +
    facet_grid(CX_ALIGNER ~ CX_READ_LENGTH) + 
    scale_colour_gradientn(colours = rainbow(11)) +
    labs(x="false discovery rate", y="sensitivity", title="ROC curve by read length 100x 300bp fragment")
  ggsave(paste0("roc_length_", variant, ".png"), width=10, height=7.5)
  
  ggplot(dtroc[dtroc$CX_READ_DEPTH==100 & dtroc$CX_READ_LENGTH==100,]) +
    aes(x=fdr, y=sens, color=log10(QUAL)) +
    geom_point() +
    facet_grid(CX_ALIGNER ~ CX_READ_FRAGMENT_LENGTH) + 
    scale_colour_gradientn(colours = rainbow(11)) +
    labs(x="false discovery rate", y="sensitivity", title="ROC curve by fragment size 100x 2x100bp")
  ggsave(paste0("roc_fragsize_", variant, ".png"), width=10, height=7.5)
}

############
# GRIDSS
#

calls <- truthlist$calls[truthlist$calls$CX_CALLER=="gridss", cbind(gridss.vcftodf(vcfs[[Id]])[vcfIndex,], .SD), by="Id"]

# assembly rate
ggplot(calls[calls$CX_READ_FRAGMENT_LENGTH==300 & calls$CX_READ_DEPTH==100,]) +
  aes(x=QUAL-ASQ-RASQ, fill=paste(pmin(AS, RAS), "/", pmax(AS, RAS), "assemblies")) +
  geom_histogram() +
  facet_wrap(~CX_READ_LENGTH) +
  scale_x_log10(limits=c(20, max(calls$QUAL-calls$ASQ-calls$RASQ))) +
  labs(title="Assembly rate 100x 300bp fragment", x="Evidence QUALity")
ggsave(paste0("roc_fragsize_", variant, ".png"), width=10, height=7.5)
ggsave("assembly_rate_QUAL.png", width=10, height=7.5)
ggplot(calls) +
  aes(x=SC+RSC+RP+BRP+BSC, fill=factor(AS+BAS)) +
  scale_y_continuous(labels=percent, expand=c(0, 0)) +  
  geom_bar(position='fill', binwidth=0.15) + 
  scale_x_log10() + 
  facet_grid(CX_READ_FRAGMENT_LENGTH ~ CX_READ_LENGTH) +
  labs(title="Assembly rate k=25", x="Reads supporting variant", fill="Assembled contigs")
ggsave("assembly_rate_reads.png", width=10, height=7.5)

ggplot(calls[calls$CX_READ_LENGTH==100,]) +
  aes(x=SC+RSC+RP+BRP+BSC, fill=paste(AS+BAS, "assemblies"), alpha=tp) +
  scale_y_continuous(labels=percent, expand=c(0, 0)) +  
  geom_bar(position='fill', binwidth=0.15) + 
  scale_x_log10() +  
  labs(title="Assembly rate 2x100bp", x="Reads supporting variant", y="Assembly rate", fill="Assembly contigs")
ggsave("assembly_rate_100bp.png", width=10, height=7.5)

# call QUAL
ggplot(calls) +
  aes(x=QUAL, y=BQ, color=hasAS, shape=tp, alpha=0.1) +
  geom_point() +
  facet_grid(CX_READ_DEPTH ~ CX_READ_FRAGMENT_LENGTH) +
  scale_x_log10() +
  scale_y_log10() + 
  labs(x="Called QUAL (including assembly)", y="Breakend QUAL", title="Assembly rate")
ggsave("assembly_rate_called_QUAL.png", width=10, height=7.5)






