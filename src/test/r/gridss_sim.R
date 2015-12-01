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

theme_set(theme_bw())

pwd <- getwd()
vcfs_all <- NULL
setwd(paste0(ifelse(as.character(Sys.info())[1] == "Windows", "W:/", "~/"), "i/data.fullmatrix"))
#setwd(paste0(ifelse(as.character(Sys.info())[1] == "Windows", "W:/", "~/"), "i/data.fullmatrix/test"))
metadata <- LoadMetadata()
vcfs_all <- LoadVcfs(metadata, existingVcfs=vcfs_all)
setwd(pwd)

vcfs_filtered <- lapply(vcfs_all, function(vcf) vcf[rowRanges(vcf)$FILTER %in% c(".", "PASS")])
vcfs_assembly_both <- lapply(vcfs_all, function(vcf) {
	if (is.na(vcf@metadata$CX_CALLER)) return(vcf)
	return(vcf[!is.na(info(vcf)$AS) & !is.na(info(vcf)$RAS),])
})

truth_all <- CalculateTruthSummary(vcfs_all, maxerrorbp=100, ignore.strand=FALSE)
truth_filtered <- CalculateTruthSummary(vcfs_filtered, maxerrorbp=100, ignore.strand=FALSE)
truth_assembly_both <- CalculateTruthSummary(vcfs_assembly_both, maxerrorbp=100, ignore.strand=FALSE)

calls_all <- truth_all$calls[, cbind(gridss.vcftodf(vcfs_all[[Id]])[vcfIndex,], .SD), by="Id"]
calls_filtered <- truth_filtered$calls[, cbind(gridss.vcftodf(vcfs_filtered[[Id]])[vcfIndex,], .SD), by="Id"]
calls_assembly_both <- truth_assembly_both$calls[, cbind(gridss.vcftodf(vcfs_assembly_both[[Id]])[vcfIndex,], .SD), by="Id"]

save.image(file="~/gridss_sim.RData")
############
# GRIDSS
#

calls <- calls_all

# assembly rate
bp_assembly_rate <- ggplot(calls[calls$tp,]) +
  aes(fill=paste(pmin(AS, RAS), "/", pmax(AS, RAS), "assemblies")) +
  facet_grid(CX_READ_FRAGMENT_LENGTH ~ CX_READ_LENGTH) +
  scale_y_continuous(labels=percent, expand=c(0, 0)) +  
  labs(title="Breakpoint assembly rate, k=25", y="", fill="Breakend assembly count")
be_assembly_rate <- ggplot(calls[calls$tp,]) +
  aes(fill=as.factor(AS)) +
  facet_grid(CX_READ_FRAGMENT_LENGTH ~ CX_READ_LENGTH) +
  scale_y_continuous(labels=percent, expand=c(0, 0)) +  
  labs(title="Breakend assembly rate, k=25", y="", fill="Breakend assembly count")
ggsave("gridss_breakpoint_assembly_rate_reads.png", width=10, height=7.5, plot=bp_assembly_rate +
         geom_bar(position='fill', binwidth=1) +
         aes(x=SR+RSR+RP+BUM+BSC) + 
         scale_x_continuous(limits=c(1, 100)) + 
         labs(x="Supporting reads"))
ggsave("gridss_breakpoint_assembly_rate_qual.png", width=10, height=7.5, plot=bp_assembly_rate +
         geom_bar(position='fill') + 
         aes(x=QUAL-ASQ-RASQ) + 
         scale_x_log10(limits=c(25, max(calls$QUAL-calls$ASQ-calls$RASQ))) +
         labs(x="Total evidence quality"))
ggsave("gridss_assembly_rate_reads.png", width=10, height=7.5, plot=be_assembly_rate +
         geom_bar(position='fill', binwidth=1) +
         aes(x=SR+RSR+RP+BUM+BSC) + 
         scale_x_continuous(limits=c(1, 100)) + 
         labs(x="Supporting reads"))
ggsave("gridss_assembly_rate_qual.png", width=10, height=7.5, plot=be_assembly_rate +
         geom_bar(position='fill') + 
         aes(x=QUAL-ASQ-RASQ) + 
         scale_x_log10(limits=c(25, max(calls$QUAL-calls$ASQ-calls$RASQ))) +
         labs(x="Total evidence quality"))
# call QUAL
ggplot(calls) +
  aes(x=QUAL, y=BQ, color=hasAS, shape=tp, alpha=0.1) +
  geom_point() +
  facet_grid(CX_READ_DEPTH ~ CX_READ_FRAGMENT_LENGTH) +
  scale_x_log10() +
  scale_y_log10() + 
  labs(x="Called QUAL (including assembly)", y="Breakend QUAL", title="Assembly rate")
ggsave("assembly_rate_called_QUAL.png", width=10, height=7.5)

### Debugging: what calls are missed by positional assembly, but picked up by subgraph assembly?
# posTruth <- truthlist_passfilters$truth[truthlist_passfilters$truth$CX_CALLER=="gridss Positional",]
# subTruth <- truthlist_passfilters$truth[truthlist_passfilters$truth$CX_CALLER=="gridss Subgraph",]
# posTruth <- posTruth[order(posTruth$vcfIndex, posTruth$CX_REFERENCE_VCF),]
# subTruth <- subTruth[order(subTruth$vcfIndex, subTruth$CX_REFERENCE_VCF),]
# missingTruth <- posTruth[!posTruth$tp & subTruth$tp,]
# dtmissing <- missingTruth[, list(count=.N), by=c("SVLEN", "SVTYPE")]
# ggplot(dtmissing) +
#   aes(y=count, x=SVLEN, color=SVTYPE) +
#   geom_line() + 
#   scale_x_log10(breaks=2**(0:16)) +  
#   theme(panel.grid.major = element_blank()) +
#   labs(y="missed events", x="Event size")
# missingTruth[missingTruth$SVTYPE=="BND",]

plot_f1_roc <- function(truthcalls, title, filename) {
  #sens <- truthcalls$truth[, list(sens=sum(tp)/.N), by=c("SVLEN", "SVTYPE", "CX_ALIGNER", "CX_ALIGNER_SOFTCLIP", "CX_CALLER", "CX_READ_DEPTH", "CX_READ_FRAGMENT_LENGTH", "CX_READ_LENGTH", "CX_REFERENCE_VCF_VARIANTS")]
  bylist <- c("CX_ALIGNER", "CX_ALIGNER_SOFTCLIP", "CX_READ_DEPTH", "CX_READ_FRAGMENT_LENGTH", "CX_READ_LENGTH", "CX_REFERENCE_VCF_VARIANTS")
  roc <- TruthSummaryToROC(truthcalls, bylist=bylist)
  roc$f1 <- 2 * roc$prec * roc$sens / (roc$prec + roc$sens)
  for (rl in unique(roc$CX_READ_LENGTH)) {
    ggplot(roc[roc$CX_READ_LENGTH==rl,]) + 
      aes(y=sens, x=fdr, color=factor(CX_READ_DEPTH)) + #aes(color=log(QUAL)) + scale_colour_gradientn(colours = rainbow(11)) +
      geom_point(size=1) + 
      facet_grid(CX_REFERENCE_VCF_VARIANTS + CX_ALIGNER ~ CX_READ_FRAGMENT_LENGTH) +
      scale_x_continuous(limits=c(0, 0.05)) + 
      labs(y="sensitivity", title=paste("ROC by call quality threshold 2x", rl, "bp", title), color="read depth")
    ggsave(paste0("gridss_roc_rl", rl, "_", filename, ".png"), width=10, height=7.5)
  }
  f1 <- data.table(roc)[, j=list(f1=max(f1)), by=bylist]
  f1$QUAL <- merge(roc, f1, by=c(bylist, "f1"))[, j=list(QUAL=max(QUAL)), by=bylist]$QUAL
  ggplot(f1) +
    aes(y=f1^4, x=CX_READ_DEPTH, color=factor(CX_READ_FRAGMENT_LENGTH), shape=CX_ALIGNER) +
    geom_line() + geom_point() +
    scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.1)^4, labels=c("0", "", "", "", "", "0.5", "", "0.7", "0.8", "0.9", "1.0")) +
    scale_x_continuous(limits=c(0, 100), breaks=c(4,8,15,30,60)) + 
    facet_grid(CX_READ_LENGTH ~ CX_REFERENCE_VCF_VARIANTS) +
    labs(y="Maximum F1 score", x="Read Depth", title=paste("gridss F-score", title), shape="aligner", color="Fragment Size")
  ggsave(paste0("gridss_maxf1_", filename, ".png"), width=10, height=7.5)
}
plot_f1_roc(truth_all, " (all calls)", "all")
plot_f1_roc(truth_filtered, " (assembly support)", "assembly")
plot_f1_roc(truth_assembly_both, " (reciprocal assembly support)", "assembly_both")



# sens/prec by assembly
bylistsensprec <- c("CX_ALIGNER", "CX_ALIGNER_SOFTCLIP", "CX_CALLER", "CX_READ_DEPTH", "CX_READ_FRAGMENT_LENGTH", "CX_READ_LENGTH", "CX_REFERENCE_VCF_VARIANTS") #"SVLEN", "SVTYPE", 
calls_all$assembly <- ifelse(calls_all$AS>0 & calls_all$RAS>0, "Both Breakends", ifelse(calls_all$AS>0, "Single Breakend", "No Assembly"))
dtprecision <- calls_all[, list(precision=sum(tp)/.N), by=c("assembly", bylistsensprec)]
save(dtprecision, file = "dtprecision.RData")

for (rl in unique(dtprecision$CX_READ_LENGTH)) {
  ggplot(dtprecision[dtprecision$CX_READ_LENGTH==rl,]) +
    aes(y=precision, x=CX_READ_DEPTH, color=assembly, shape=CX_ALIGNER, linetype=CX_ALIGNER) +
    geom_line() + geom_point() +
    facet_grid(CX_READ_LENGTH + CX_READ_FRAGMENT_LENGTH ~ CX_REFERENCE_VCF_VARIANTS) +
    labs(title=paste("Precision by assembly status 2x", rl, "bp"), y="Read Depth") + scale_x_log10()
  ggsave(paste0("gridss_assembly_precision_rl", rl, ".png"), width=10, height=7.5)
}

dtsens <- calls_all[, list(count=sum(tp)), by=c("assembly", bylistsensprec)]
dtsens <- merge(dtsens, truth_all$truth[, list(total=.N), bylistsensprec], by=bylistsensprec)
dtsens$sensitivity <- dtsens$count / dtsens$total

for (rl in unique(dtsens$CX_READ_LENGTH)) {
  ggplot(dtsens[dtsens$CX_READ_LENGTH==rl,]) +
    aes(y=sensitivity, x=CX_READ_DEPTH, color=assembly, shape=CX_ALIGNER, linetype=CX_ALIGNER) +
    geom_line() + geom_point() +
    facet_grid(CX_READ_LENGTH + CX_READ_FRAGMENT_LENGTH ~ CX_REFERENCE_VCF_VARIANTS) +
    labs(title=paste("Sensitivity by assembly status 2x", rl, "bp"), y="Read Depth") + scale_x_log10()
  ggsave(paste0("gridss_assembly_sensitivity_rl", rl, ".png"), width=10, height=7.5)
}













