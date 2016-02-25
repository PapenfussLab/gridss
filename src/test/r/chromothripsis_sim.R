##################
# Processing steps (expected runtime: 191h wall, 205h CPU)
##################
# - Follow common.R processing steps
# - ./gehumanchr12variants.sh fastcompare
# - ./genreads.sh fastcompare
# - ./alignbam.sh fastcompare
# - ./sortbam.sh fastcompare
# - ./gendownsample.sh fastcompare
# - ./call_gridss.sh fastcompare
# - ./call_breakdancer.sh fastcompare
# - ./call_delly.sh fastcompare
# - ./call_pindel.sh fastcompare
# - ./call_socrates.sh fastcompare
# - ./call_lumpy.sh fastcompare
# - Run this script (R requires at least 6GB of memory for this script)


# http://www.nature.com/nmeth/authors/submit/index.html#formats
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
source("common.R")

#########################
# Simulation Comparison Data
########################
pwd <- getwd()
setwd(paste0(ifelse(as.character(Sys.info())[1] == "Windows", "W:/", "~/"), "i/data.fastcompare")) # setwd("W:/i/data.fastcompare")
metadata <- LoadMetadata()
metadata <- WithVcf(metadata)
metadata <- main_callers_subset(metadata)
vcfs <- LoadVcfs(metadata)
setwd(pwd)

vcfs <- lapply(vcfs, function(vcf) {
  caller <- str_extract(attr(vcf, "sourceMetadata")$CX_CALLER, "^[^/]+")
  vcf <- withqual(vcf, caller)
  if (!is.na(caller) && caller == "cortex") {
    vcf <- cleanCortex(vcf)
  }
  if (str_detect(attr(vcf, "sourceMetadata")$CX_REFERENCE_VCF_VARIANTS, "BP")) {
    # filtering small events calls to remove spurious indels caused by sequence homology around breakpoints
    callSize <- eventSize(vcf)
    callSize[is.na(callSize)] <- 1000000000
    vcf <- vcf[callSize > 50,]
  }
  return(vcf)
})
vcfs_passfilters <- lapply(vcfs, function(vcf) vcf[rowRanges(vcf)$FILTER %in% c(".", "PASS")])

#Rprof("C:/dev/Rprof.out")
truthlist_all <- CalculateTruthSummary(vcfs, maxerrorbp=100, ignore.strand=TRUE) # breakdancer does not specify strand
truthlist_passfilters <- CalculateTruthSummary(vcfs_passfilters, maxerrorbp=100, ignore.strand=TRUE) # breakdancer does not specify strand

truthlist_all$calls <- truthlist_all$calls
truthlist_all$truth <- truthlist_all$truth
truthlist_passfilters$calls <- truthlist_passfilters$calls
truthlist_passfilters$truth <- truthlist_passfilters$truth

# Sensitivity
dtsenssize_all <- truthlist_all$truth[, list(sens=sum(tp)/.N), by=c("SVLEN", "SVTYPE", "CX_ALIGNER", "CX_ALIGNER_SOFTCLIP", "CX_CALLER", "CX_READ_DEPTH", "CX_READ_FRAGMENT_LENGTH", "CX_READ_LENGTH", "CX_REFERENCE_VCF_VARIANTS")]
dtsenssize_passfilters <- truthlist_passfilters$truth[, list(sens=sum(tp)/.N), by=c("SVLEN", "SVTYPE", "CX_ALIGNER", "CX_ALIGNER_SOFTCLIP", "CX_CALLER", "CX_READ_DEPTH", "CX_READ_FRAGMENT_LENGTH", "CX_READ_LENGTH", "CX_REFERENCE_VCF_VARIANTS")]
dtsenssize_all$Filter <- prettyFilter(text_all_calls)
dtsenssize_passfilters$Filter <- prettyFilter(text_default_calls)
dtsenssize_all$Caller <- dtsenssize_all$CX_CALLER
dtsenssize_passfilters$Caller <- dtsenssize_passfilters$CX_CALLER
dtsenssize <- rbind(dtsenssize_all, dtsenssize_passfilters)
# Precision
dtprec_all <- truthlist_all$truth[, list(prec=sum(tp)/.N), by=c("CX_ALIGNER", "CX_ALIGNER_SOFTCLIP", "CX_CALLER", "CX_READ_DEPTH", "CX_READ_FRAGMENT_LENGTH", "CX_READ_LENGTH", "CX_REFERENCE_VCF_VARIANTS")]
dtprec_passfilters <- truthlist_passfilters$truth[, list(prec=sum(tp)/.N), by=c("CX_ALIGNER", "CX_ALIGNER_SOFTCLIP", "CX_CALLER", "CX_READ_DEPTH", "CX_READ_FRAGMENT_LENGTH", "CX_READ_LENGTH", "CX_REFERENCE_VCF_VARIANTS")]
dtprec_all$Filter <- prettyFilter(text_all_calls)
dtprec_passfilters$Filter <- prettyFilter(text_default_calls)
dtprec_all$Caller <- dtprec_all$CX_CALLER
dtprec_passfilters$Caller <- dtprec_passfilters$CX_CALLER
dtprec <- rbind(dtprec_all, dtprec_passfilters)

# plot_sens <- ggplot(dtsenssize[!(dtsenssize$CX_REFERENCE_VCF_VARIANTS %in% c("hetBP", "hetBP_SINE"))]) + # & dtsenssize$Caller %in% c("gridss Positional", "gridss Subgraph")
#   aes(y=sens, x=SVLEN, shape=factor(CX_READ_FRAGMENT_LENGTH), color=factor(Filter)) +
#   geom_line() + 
#   facet_grid(CX_ALIGNER + Caller ~ CX_READ_FRAGMENT_LENGTH + CX_READ_LENGTH + CX_READ_DEPTH + CX_REFERENCE_VCF_VARIANTS) +
#   scale_x_log10(breaks=2**(0:16)) +  
#   theme(panel.grid.major = element_blank()) +
#   labs(y="sensitivity", x="Event size", title="Sensitivity")
# saveplot("sim_sens_size", plot_sens, width=10, height=7.5)
# # rescaled to show close to 1
# saveplot("sim_sens_size_nonlinearscale", plot_sens + aes(y=sens**5) + scale_y_continuous(breaks=seq(0,1,0.1)**5, labels=seq(0,1,0.1)), width=10, height=7.5)
# 
# ggplot(truthlist_all$truth) + aes(x=errorsize) + facet_grid(CX_REFERENCE_VCF_VARIANTS~Caller) + geom_histogram() + scale_y_log10()
# ggplot(truthlist_all$calls) + aes(x=partialtp) + facet_grid(CX_REFERENCE_VCF_VARIANTS~Caller) + geom_histogram()

scale_x_svlen <- scale_x_continuous(breaks=2**(0:16),
                                    labels=c("1", "2", "4", "8", "16", "32", "64", "128", "256", "512", "1k", "2k", "4k", "8k", "16k", "32k", "64k"),
                                    minor_breaks=log10(c(1,2,3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 20, 24, 28, 32, 48, 64, 80, 96, 112, 128, 160, 192, 224, 256, 288, 320, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536)),
                                    trans="log10")
dtsss <- dtsenssize
dtsss <- dtsss[!is.na(dtsss$SVLEN),]
dtsss[dtsss$CX_REFERENCE_VCF_VARIANTS=="hetDEL",]$CX_REFERENCE_VCF_VARIANTS <- "Deletion"
dtsss[dtsss$CX_REFERENCE_VCF_VARIANTS=="hetINS",]$CX_REFERENCE_VCF_VARIANTS <- "Insertion"
dtsss[dtsss$CX_REFERENCE_VCF_VARIANTS=="hetDUP",]$CX_REFERENCE_VCF_VARIANTS <- "Tandem Duplication"
dtsss[dtsss$CX_REFERENCE_VCF_VARIANTS=="hetINV",]$CX_REFERENCE_VCF_VARIANTS <- "Inversion"
for (rl in unique(dtsss$CX_READ_LENGTH)) {
for (rd in unique(dtsss$CX_READ_DEPTH)) {
for (fragsize in unique(dtsss$CX_READ_FRAGMENT_LENGTH)) {
  dt <- dtsss
  dt <- dt[dt$CX_READ_DEPTH==rd & dt$CX_READ_FRAGMENT_LENGTH==fragsize & dt$CX_READ_LENGTH==rl,]
  plot_sim_sens_simple <- ggplot(dt[!(dt$Caller=="gridss"& dt$Filter=="All calls")]) +
    aes(y=sens, x=SVLEN, color=Caller, shape=Caller, linetype=Filter) +
    geom_line() + 
    geom_point(data=dt[,.SD[chull(.SD$sens, log(.SD$SVLEN)),], by=c("Caller", "Filter", "CX_REFERENCE_VCF_VARIANTS")]) + 
    geom_point() +
    geom_line(data=dt[dt$Caller=="gridss",], size=1.2) + 
    facet_wrap( ~ CX_REFERENCE_VCF_VARIANTS) +
    scale_x_svlen + 
    aes(y=sens**4) + scale_y_power4 + 
    #scale_color_brewer(type="qual", palette=2) +
    #theme(panel.grid.major = element_blank()) +
    labs(y="Sensitivity", x="Event size", title=paste0("Sensitivity ", rd, "x coverage"))
  plot_sim_sens_simple
  saveplot(paste0("sim_sens_simple_rl", rl, "_", fragsize, "_", rd, "x"), width=300, height=150, units=c("mm"))
  plot_sim_sens_simple + aes(y=sens) + scale_y_continuous(limits=c(0,1))
  saveplot(paste0("sim_sens_simple_linear_rl", rl, "_", fragsize, "_", rd, "x"), width=300, height=150, units=c("mm"))
}}}
for (rl in unique(dtsenssize$CX_READ_LENGTH)) {
for (fragsize in unique(dtsenssize$CX_READ_FRAGMENT_LENGTH)) {
    dt <- dtsss
    dt <- dt[dt$CX_READ_FRAGMENT_LENGTH==fragsize & dt$CX_READ_LENGTH==rl,]
    plot_sim_sens_simple <- ggplot(dt[!(dt$Caller=="gridss"& dt$Filter=="All calls")]) +
      aes(y=sens, x=SVLEN, color=Caller, shape=Caller, linetype=Filter) +
      geom_line() + 
      geom_point(data=dt[,.SD[chull(.SD$sens, log(.SD$SVLEN)),], by=c("Caller", "Filter", "CX_REFERENCE_VCF_VARIANTS")]) + 
      geom_point() +
      geom_line(data=dt[dt$Caller=="gridss",], size=1.2) + 
      
      facet_grid(CX_REFERENCE_VCF_VARIANTS ~ CX_READ_DEPTH) +
      scale_x_svlen +
      aes(y=sens**4) + scale_y_power4 + 
      #scale_color_brewer(type="qual", palette=2) +
      labs(y="Sensitivity", x="Event size", title=paste0("Sensitivity by read depth"))
    plot_sim_sens_simple
    saveplot(paste0("sim_sens_simple_rl", rl, "_", fragsize), width=800, height=300, units=c("mm"))
    plot_sim_sens_simple + aes(y=sens) + scale_y_continuous(limits=c(0,1))
    saveplot(paste0("sim_sens_simple_linear_rl", rl, "_", fragsize), width=800, height=300, units=c("mm"))
}}
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
    aes(y=prec, x=Caller, fill=Caller) +
    geom_bar(stat="identity") + 
    facet_wrap( ~ CX_REFERENCE_VCF_VARIANTS) +
    #scale_color_brewer(type="qual", palette=2) +
    labs(y="Precision", title=paste0("Precision ", rd, "x coverage"))
  saveplot(paste0("sim_prec_simple_rl", rl, "_", fragsize, "_", rd, "x"), width=300, height=150, units=c("mm"))
}}}

############
# Error in called position
#
ggplot(truthlist_all$truth[!is.na(truthlist_all$truth$errorsize),]) +
  aes(x=errorsize, color=factor(SVLEN)) +
  geom_density(aes(y=..density..)) +
  facet_grid(CX_ALIGNER + Caller ~ CX_READ_FRAGMENT_LENGTH + CX_READ_LENGTH + CX_READ_DEPTH, scale="free")


############
# ROC
#
bylist <- c("CX_ALIGNER", "CX_ALIGNER_SOFTCLIP", "CX_CALLER", "CX_READ_DEPTH", "CX_READ_FRAGMENT_LENGTH", "CX_READ_LENGTH", "CX_REFERENCE_VCF_VARIANTS")
dtrocunfiltered <- TruthSummaryToROC(truthlist_all, bylist=bylist)
dtrocunfiltered$Filter <- prettyFilter(text_all_calls)
dtrocunfiltered$Caller <- dtrocunfiltered$CX_CALLER
dtrocfiltered <- TruthSummaryToROC(truthlist_passfilters, bylist=bylist)
dtrocfiltered$Filter <- prettyFilter(text_default_calls)
dtrocfiltered$Caller <- dtrocfiltered$CX_CALLER
dtrocall <- rbind(dtrocunfiltered, dtrocfiltered)

# plot_roc <- ggplot(dtrocall[dtrocall$CX_REFERENCE_VCF_VARIANTS %in% c("hetBP", "hetBP_SINE"),]) + 
#   aes(y=sens, x=fp+1, color=log(QUAL), shape=Filter) +
#   scale_colour_gradientn(colours = rainbow(11)) +
#   geom_line() + 
#   geom_point(size=0.25) + 
#   facet_grid(CX_ALIGNER + Caller ~ CX_REFERENCE_VCF_VARIANTS + CX_READ_FRAGMENT_LENGTH + CX_READ_LENGTH + CX_READ_DEPTH) +
#   scale_x_log10() + 
#   labs(y="sensitivity", title="ROC by call quality threshold")
# saveplot("sim_roc_full", plot_roc, width=10, height=7.5)
# saveplot("sim_roc_full_nonlinear_scale", plot_roc + scale_x_continuous(limits=c(0, 100)) + aes(x=fp) + aes(y=sens**5) + scale_y_continuous(breaks=seq(0,1,0.1)**5, labels=seq(0,1,0.1)), width=10, height=7.5)

# for (variant in unique(dtrocall$CX_REFERENCE_VCF_VARIANTS)) {
#   dtroc <- dtrocall[dtrocall$CX_REFERENCE_VCF_VARIANTS == variant,]
#   plot_roc_all <- ggplot(dtroc) +
#     aes(y=sens, shape=factor(CX_READ_FRAGMENT_LENGTH), color=factor(CX_READ_DEPTH), linetype=Filter) +
#     geom_point(size=1) + 
#     facet_grid(CX_ALIGNER + Caller ~ CX_REFERENCE_VCF_VARIANTS + CX_READ_LENGTH) +
#     labs(y="sensitivity", title="ROC by call quality threshold")
#   plot_roc_all + aes(x=prec) + labs(x="precision") + scale_x_reverse()
#   saveplot(paste0("sim_roc_full_precision_", variant, ""), width=10, height=7.5)
#   plot_roc_all + aes(x=fdr) + labs(x="false discovery rate")
#   saveplot(paste0("sim_roc_full_fdr_", variant, ""), width=10, height=7.5)
#   plot_roc_all + aes(x=fp) + labs(x="total false positives")
#   saveplot(paste0("sim_roc_full_fp_", variant, ""), width=10, height=7.5)
#   
#   ggplot(dtroc[dtroc$CX_READ_FRAGMENT_LENGTH==300 & dtroc$CX_READ_LENGTH==100,]) +
#     aes(x=fdr, y=sens, color=log10(QUAL), shape=CX_ALIGNER) +
#     geom_point() +
#     facet_grid(Caller + CX_REFERENCE_VCF_VARIANTS + Filter ~ CX_READ_DEPTH) + 
#     scale_colour_gradientn(colours = rainbow(11)) +
#     labs(x="false discovery rate", y="sensitivity", title="ROC curve by read depth 2x100bp 300bp fragment")
#   saveplot(paste0("sim_roc_depth_", variant, ""), width=10, height=7.5)
#   
#   ggplot(dtroc[dtroc$CX_READ_FRAGMENT_LENGTH==300 & dtroc$CX_READ_DEPTH==100,]) +
#     aes(x=fdr, y=sens, color=log10(QUAL)) +
#     geom_point() +
#     facet_grid(CX_ALIGNER + Filter ~ CX_READ_LENGTH) + 
#     scale_colour_gradientn(colours = rainbow(11)) +
#     labs(x="false discovery rate", y="sensitivity", title="ROC curve by read length 100x 300bp fragment")
#   saveplot(paste0("sim_roc_length_", variant, ""), width=10, height=7.5)
#   
#   ggplot(dtroc[dtroc$CX_READ_DEPTH==100 & dtroc$CX_READ_LENGTH==100,]) +
#     aes(x=fdr, y=sens, color=log10(QUAL)) +
#     geom_point() +
#     facet_grid(CX_ALIGNER + Filter ~ CX_READ_FRAGMENT_LENGTH) + 
#     scale_colour_gradientn(colours = rainbow(11)) +
#     labs(x="false discovery rate", y="sensitivity", title="ROC curve by fragment size 100x 2x100bp")
#   saveplot(paste0("sim_roc_fragsize_", variant, ""), width=10, height=7.5)
# }
for (rl in unique(dtrocall$CX_READ_LENGTH)) {
for (rd in unique(dtrocall$CX_READ_DEPTH)) {
for (fragsize in unique(dtrocall$CX_READ_FRAGMENT_LENGTH)) {
  dt <- dtrocall
  dt <- dt[dt$CX_READ_DEPTH==rd & dt$CX_READ_FRAGMENT_LENGTH==fragsize & dt$CX_READ_LENGTH==rl,]
  dt <- dt[dt$CX_REFERENCE_VCF_VARIANTS %in% c("hetBP","hetBP_SINE"),]
  dt[dt$CX_REFERENCE_VCF_VARIANTS=="hetBP",]$CX_REFERENCE_VCF_VARIANTS <- "Fusion (15,000 variants)"
  dt[dt$CX_REFERENCE_VCF_VARIANTS=="hetBP_SINE",]$CX_REFERENCE_VCF_VARIANTS <- "Fusion at SINE/ALU (14,554 variants)"
  plot_sim_roc <- ggplot(dt[!(dt$Caller=="gridss"& dt$Filter=="All calls")]) +
    aes(y=sens, x=fp/2+1, color=Caller, shape=Caller, linetype=Filter) +
    #scale_color_brewer(type="qual", palette=2) +
    facet_wrap(~ CX_REFERENCE_VCF_VARIANTS, ncol=1) + 
    geom_line() + 
    geom_point(data=dt[,.SD[chull(log(.SD$fp/2+1), .SD$sens),], by=c("Caller", "Filter", "CX_REFERENCE_VCF_VARIANTS")]) + 
    geom_point() +
    geom_line(data=dt[dt$Caller=="gridss",], size=1.2) + 
    scale_x_log10(breaks=c(1, 11, 101, 1001, 10001), labels=c("0", "10", "100", "1000", "10000"), minor_breaks=c(6, 16, 151, 1501, 15001)) +
    aes(y=sens**4) + scale_y_power4 + 
    labs(y="Sensitivity", x="False Positives", title="ROC by call quality threshold")
  plot_sim_roc
  saveplot(paste0("sim_roc_rl", rl, "_", fragsize, "_", rd, "x"), width=120, height=150, units=c("mm"))
  plot_sim_roc + aes(y=sens) + scale_y_continuous(limits=c(0,1))
  saveplot(paste0("sim_roc_linear_rl", rl, "_", fragsize, "_", rd, "x"), width=120, height=150, units=c("mm"))
}}}
for (rl in unique(dtrocall$CX_READ_LENGTH)) {
for (fragsize in unique(dtrocall$CX_READ_FRAGMENT_LENGTH)) {
  dt <- dtrocall
  dt <- dt[dt$CX_READ_FRAGMENT_LENGTH==fragsize & dt$CX_READ_LENGTH==rl,]
  dt <- dt[dt$CX_REFERENCE_VCF_VARIANTS %in% c("hetBP","hetBP_SINE"),]
  dt[dt$CX_REFERENCE_VCF_VARIANTS=="hetBP",]$CX_REFERENCE_VCF_VARIANTS <- "Fusion (15,000 variants)"
  dt[dt$CX_REFERENCE_VCF_VARIANTS=="hetBP_SINE",]$CX_REFERENCE_VCF_VARIANTS <- "Fusion at SINE/ALU (14,554 variants)"
  dt$Caller <- dt$Caller
  plot_sim_roc <- ggplot(dt[!(dt$Caller=="gridss"& dt$Filter=="All calls")]) +
    aes(y=sens, x=fp/2+1, color=Caller, shape=Caller, linetype=Filter) +
    #scale_color_brewer(type="qual", palette=2) +
    facet_grid(CX_REFERENCE_VCF_VARIANTS ~ CX_READ_DEPTH) + 
    geom_line() + 
    geom_point(data=dt[,.SD[chull(log(.SD$fp/2+1), .SD$sens),], by=c("Caller", "Filter", "CX_REFERENCE_VCF_VARIANTS", "CX_READ_DEPTH")]) + 
    geom_point() +
    geom_line(data=dt[dt$Caller=="gridss",], size=1.2) + 
    scale_x_log10(breaks=c(1, 11, 101, 1001, 10001), labels=c("0", "10", "100", "1000", "10000"), minor_breaks=c(6, 16, 151, 1501, 15001)) +
    aes(y=sens**4) + scale_y_power4 + 
    labs(y="Sensitivity", x="False Positives", title="ROC by read depth")
  plot_sim_roc
  saveplot(paste0("sim_roc_rl", rl, "_", fragsize), width=300, height=150, units=c("mm"))
  plot_sim_roc + aes(y=sens) + scale_y_continuous(limits=c(0,1))
  saveplot(paste0("sim_roc_linear_rl", rl, "_", fragsize), width=300, height=150, units=c("mm"))
}}







#save.image(file="~/chromothripsis_sim.RData")
#load(file="~/chromothripsis_sim.RData")









