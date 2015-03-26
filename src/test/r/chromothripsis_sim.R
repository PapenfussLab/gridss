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
cl <- makeCluster(detectCores())
registerDoSNOW(cl)  


callingMargin <- 32

theme_set(theme_bw())

# Load simulation data
pwd <- getwd()
setwd("W:/i/data.fastcompare")
metadata <- LoadMetadata()
metadata <- metadata[metadata$CX_REFERENCE_VCF_VARIANTS=="homSINE",]
vcfs <- LoadVcfs(metadata)
#vcfs <- vcfs[sapply(vcfs, function(vcf) { return(!is.na(attr(vcf, "metadata")$Id))})]
setwd(pwd)

Rprof("C:/dev/Rprof.out")
truthlist <- CalculateTruthSummary(vcfs, maxgap=100)
Rprof(NULL)
calls <- data.table(truthlist$calls)
truth <- data.table(truthlist$truth)
dtroc <- NULL
bylist <- c("CX_ALIGNER", "CX_ALIGNER_SOFTCLIP", "CX_CALLER", "CX_READ_DEPTH", "CX_READ_FRAGMENT_LENGTH", "CX_READ_LENGTH", "CX_REFERENCE_VCF_VARIANTS")
rocPoints <- 128
for (minQUAL in c(0.001, exp(0:(rocPoints+1)*log(max(calls$QUAL))/rocPoints))) { # 128 data points
  dttmp <- merge(
    calls[, list(tpcall=sum(tp & QUAL>=minQUAL), fp=sum(!tp & QUAL>=minQUAL)), by=bylist],
    truth[, list(tp=sum(QUAL>=minQUAL), fn=sum(QUAL<minQUAL)), by=bylist],
    by=bylist
  )
  dttmp$QUAL <- minQUAL
  dtroc <- rbind(dtroc, dttmp)
  dttmp <- NULL
}
dtroc$sens <- dtroc$tp/(dtroc$tp+dtroc$fn)
dtroc$prec <- dtroc$tp/(dtroc$tp+dtroc$fp)
dtroc$prec[is.na(dtroc$prec)] <- 1
dtroc$fdr <- dtroc$fp/(dtroc$tp+dtroc$fp)
dtroc$fdr[is.na(dtroc$fdr)] <- 0

# alternate approach would be to sort by sens ; ddply group; accumulate on max(prec or fdr)
dtRocSensPrecHull <- ddply(dtroc, as.quoted(~ CX_CALLER + CX_READ_FRAGMENT_LENGTH + CX_READ_DEPTH + CX_ALIGNER + CX_READ_LENGTH + CX_REFERENCE_VCF_VARIANTS), function(df) df[chull(df$sens, df$prec),])
dtRocSensFdrHull <- ddply(dtroc, as.quoted(~ CX_CALLER + CX_READ_FRAGMENT_LENGTH + CX_READ_DEPTH + CX_ALIGNER + CX_READ_LENGTH + CX_REFERENCE_VCF_VARIANTS), function(df) df[chull(df$sens, df$fdr),])

plot_roc_all <- ggplot(dtroc) +
  aes(y=sens, color=factor(CX_READ_DEPTH), shape=factor(CX_READ_FRAGMENT_LENGTH)) +
  geom_point(size=1) + 
  facet_grid(CX_CALLER + CX_ALIGNER ~ CX_REFERENCE_VCF_VARIANTS + CX_READ_LENGTH) +
  labs(y="sensitivity", title="ROC by call QUALity threshold")
plot_roc_all + aes(x=prec) + labs(x="precision") + scale_x_reverse()
ggsave("roc_all_precision.png", width=10, height=7.5)
plot_roc_all + aes(x=fdr) + labs(x="false discovery rate")
ggsave("roc_all_fdr.png", width=10, height=7.5)
plot_roc_all + aes(x=fp) + labs(x="total false positives")
ggsave("roc_all_fp.png", width=10, height=7.5)


ggplot(dtroc[dtroc$CX_READ_FRAGMENT_LENGTH==300 & dtroc$CX_READ_LENGTH==100,]) +
  aes(x=fdr, y=sens, color=log10(QUAL), shape=CX_ALIGNER) +
  geom_point() +
  facet_grid(CX_CALLER + CX_REFERENCE_VCF_VARIANTS ~ CX_READ_DEPTH) + 
  scale_colour_gradientn(colours = rainbow(11)) +
  labs(x="false discovery rate", y="sensitivity", title="ROC curve by read depth 2x100bp 300bp fragment")
ggsave("roc_read_depth.png", width=10, height=7.5)

ggplot(dtroc[dtroc$CX_READ_FRAGMENT_LENGTH==300 & dtroc$CX_READ_DEPTH==100,]) +
  aes(x=fdr, y=sens, color=log10(QUAL)) +
  geom_point() +
  facet_grid(CX_ALIGNER ~ CX_READ_LENGTH) + 
  scale_colour_gradientn(colours = rainbow(11)) +
  labs(x="false discovery rate", y="sensitivity", title="ROC curve by read length 100x 300bp fragment")
ggsave("roc_read_length.png", width=10, height=7.5)

ggplot(dtroc[dtroc$CX_READ_DEPTH==100 & dtroc$CX_READ_LENGTH==100,]) +
  aes(x=fdr, y=sens, color=log10(QUAL)) +
  geom_point() +
  facet_grid(CX_ALIGNER ~ CX_READ_FRAGMENT_LENGTH) + 
  scale_colour_gradientn(colours = rainbow(11)) +
  labs(x="false discovery rate", y="sensitivity", title="ROC curve by fragment size 100x 2x100bp")
ggsave("roc_fragment_size.png", width=10, height=7.5)


# assembly rate
ggplot(calls[calls$CX_READ_FRAGMENT_LENGTH==300 & calls$CX_READ_DEPTH==100,]) +
  aes(x=QUAL-ASQ-RASQ, fill=paste(pmin(AS, RAS), "/", pmax(AS, RAS), "assemblies")) +
  geom_histogram() +
  facet_wrap(~CX_READ_LENGTH) +
  scale_x_log10(limits=c(20, max(calls$QUAL-calls$ASQ-calls$RASQ))) +
  labs(title="Assembly rate 100x 300bp fragment", x="Evidence QUALity")
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
  aes(x=SC+RSC+RP+BRP+BSC, fill=paste(AS+BAS, "assemblies")) +
  scale_y_continuous(labels=percent, expand=c(0, 0)) +  
  geom_bar(position='fill', binwidth=0.15) + 
  scale_x_log10() +  
  labs(title="Assembly rate 2x100bp", x="Reads supporting variant", y="Assembly rate", fill="Assembly contigs")
ggsave("assembly_rate_100bp.png", width=10, height=7.5)

# call QUAL
ggplot(calls) +
  aes(x=QUAL, y=BQ, color=hasAS, shape=tp) +
  geom_point() +
  facet_grid(CX_READ_DEPTH ~ CX_READ_FRAGMENT_LENGTH) +
  scale_x_log10() +
  scale_y_log10() + 
  labs(x="Called QUAL (including assembly)", y="Breakend QUAL", title="Assembly rate")
ggsave("assembly_rate_called_QUAL.png", width=10, height=7.5)







