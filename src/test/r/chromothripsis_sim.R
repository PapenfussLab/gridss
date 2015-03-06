source("../../main/r/libgridss.R")
source("libneochromosome.R")
source("../../../../../dev/ws/indelappraisal/libindelappraisal.R")
library(data.table)
library(stringr)

callingMargin <- 32

theme_set(theme_bw())

# Load simulation data
pwd <- getwd()
setwd("E:/gridss/data.chromothripsis")
metadata <- LoadMetadata()
vcfs <- LoadBreakpointVcfs(metadata)
setwd(pwd)

truthlist <- LoadBreakpointTruthSummary(vcfs)
calls <- data.table(truthlist$calls)
truth <- data.table(truthlist$truth)
dtroc <- NULL
bylist <- c("CX_ALIGNER", "CX_ALIGNER_SOFTCLIP", "CX_CALLER", "CX_READ_DEPTH", "CX_READ_FRAGMENT_LENGTH", "CX_READ_LENGTH")
for (minQual in 1.11**(0:100)) { # 2**(0:16)) {
  dttmp <- merge(
    calls[, list(tpcall=sum(tp & QUAL>=minQual), fp=sum(!tp & QUAL>=minQual)), by=bylist],
    truth[, list(tp=sum(qual>=minQual), fn=sum(qual<minQual)), by=bylist],
    by=bylist
  )
  dttmp$qual <- minQual
  dtroc <- rbind(dtroc, dttmp)
}

plot_roc_all <- ggplot(dtroc) +
  aes(y=tp/(tp+fn), color=factor(CX_READ_FRAGMENT_LENGTH), shape=factor(CX_READ_DEPTH)) +
  geom_point(size=1) + 
  facet_grid(CX_ALIGNER ~ CX_READ_LENGTH) +
  labs(y="sensitivity", title="ROC by call quality threshold")
plot_roc_all + aes(x=tp/(tp+fp)) + labs(x="precision") + scale_x_reverse()
ggsave("roc_all_precision.png", width=10, height=7.5)
plot_roc_all + aes(x=fp/(tp+fp)) + labs(x="false positive rate")
ggsave("roc_all_fdr.png", width=10, height=7.5)

ggplot(dtroc[dtroc$CX_READ_FRAGMENT_LENGTH==400 & dtroc$==100,]) +
  aes(x=fp/(tp+fp), y=tp/(tp+fn), linetype=factor(CX_ALIGNER), color=log10(qual), shape=CX_ALIGNER) +
  geom_point() +
  facet_wrap( ~ CX_READ_DEPTH) + 
  scale_colour_gradientn(colours = rainbow(11)) +
  labs(x="false positive rate", y="sensitivity", title="ROC curve by read depth 2x100bp 400bp fragment")
ggsave("roc_read_depth.png", width=10, height=7.5)

ggplot(dtroc[dtroc$CX_READ_FRAGMENT_LENGTH==400 & dtroc$CX_READ_DEPTH==100,]) +
  aes(x=fp/(tp+fp), y=tp/(tp+fn), color=log10(qual)) +
  geom_point() +
  facet_grid(CX_ALIGNER ~ CX_READ_LENGTH) + 
  scale_colour_gradientn(colours = rainbow(11)) +
  labs(x="false positive rate", y="sensitivity", title="ROC curve by read length 100x 400bp fragment")
ggsave("roc_read_length.png", width=10, height=7.5)

ggplot(dtroc[dtroc$CX_READ_DEPTH==100 & dtroc$CX_READ_LENGTH==100,]) +
  aes(x=fp/(tp+fp), y=tp/(tp+fn), color=log10(qual)) +
  geom_point() +
  facet_grid(CX_ALIGNER ~ CX_READ_FRAGMENT_LENGTH) + 
  scale_colour_gradientn(colours = rainbow(11)) +
  labs(x="false positive rate", y="sensitivity", title="ROC curve by fragment size 100x 2x100bp")
ggsave("roc_fragment_size.png", width=10, height=7.5)


# assembly rate
ggplot(calls[calls$CX_READ_FRAGMENT_LENGTH==400 & calls$CX_READ_DEPTH==100,]) +
  aes(x=QUAL-ASQ-RASQ, fill=paste(pmin(AS, RAS), "/", pmax(AS, RAS), "assemblies")) +
  geom_histogram() +
  facet_wrap(~CX_READ_LENGTH) +
  scale_x_log10(limits=c(20, max(calls$QUAL-calls$ASQ-calls$RASQ))) +
  labs(title="Assembly rate 100x 400bp fragment", x="Evidence quality")
ggsave("assembly_rate_qual.png", width=10, height=7.5)
ggplot(calls) +
  aes(x=SC+RSC+RP+BRP+BSC, fill=paste(AS+BAS, "assemblies")) +
  scale_y_continuous(labels=percent, expand=c(0, 0)) +  
  geom_bar(position='fill', binwidth=0.15) + 
  scale_x_log10() + 
  facet_wrap(~CX_READ_LENGTH) +
  labs(title="Assembly rate 100x 400bp fragment", x="Reads supporting variant", fill="Assembly contigs")
ggsave("assembly_rate_reads.png", width=10, height=7.5)

ggplot(calls[calls$CX_READ_LENGTH==100,]) +
  aes(x=SC+RSC+RP+BRP+BSC, fill=paste(AS+BAS, "assemblies")) +
  scale_y_continuous(labels=percent, expand=c(0, 0)) +  
  geom_bar(position='fill', binwidth=0.15) + 
  scale_x_log10() +  
  labs(title="Assembly rate 2x100bp", x="Reads supporting variant", y="Assembly rate" fill="Assembly contigs")
ggsave("assembly_rate_100bp.png", width=10, height=7.5)

# call qual
ggplot(calls) +
  aes(x=QUAL, y=BQ, color=hasAS, shape=tp) +
  geom_point() +
  facet_grid(CX_READ_DEPTH ~ CX_READ_FRAGMENT_LENGTH) +
  scale_x_log10() +
  scale_y_log10() + 
  labs(x="Called qual (including assembly)", y="Breakend qual", title="Assembly rate")
ggsave("assembly_rate_called_qual.png", width=10, height=7.5)







