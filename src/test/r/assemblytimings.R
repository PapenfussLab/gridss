##################
# Processing steps
##################
# - Follow assembly-method.R processing steps
# - Run this script (updating to the correct timings.csv location)


library(ggplot2)
require(reshape)
# generate input with
# cut -d";" -f 8 gridss_vis/debruijn.assembly.metrics.* | cut -b 8- | grep -Ev '[^0-9,]' > timings.csv


timings <- c(
  "KmerConstruction",
  "GeneratePathGraph",
  "Shrink",
  "Collapse",
  "SplitOutReferencePaths",
  "CalcNonReferenceSubgraphs",
  "AssemblyNonReferenceContigs",
  "ToAssemblyEvidence",
  "Assembly")
df <- read.csv("W:/na12878/gridss_vis/timings.csv", header=FALSE, col.names=timings, colClasses=rep("integer", length(timings)))
df$total <- rowSums(as.matrix(df))
# filter out the rows that are due to java starting to runn out of memory
df <- df[df$total < 600000, ]

totalByStage <- colSums(as.matrix(df))

mdf <- melt(df, variable_name='stage')
colnames(mdf)[2] <- "time"

ggplot(mdf, aes(x=time)) + geom_histogram() + facet_wrap(~ stage, scales="free") + scale_x_log10()


## buffer sizes
library(reshape2)
dtbuffer <- read.csv(file="~/na12878/platinum/gridss_vis/positional-chr12-Forward.csv", header=TRUE)
dtbuffers <- dtbuffer[,c("trackerActive","supportProcessedSize","aggregateProcessedSize","aggregateQueueSize","aggregateActiveSize","pathNodeProcessedSize","pathNodeActiveSize","pathNodeEdgeLookupSize","pathNodePathLookupSize","collapseProcessedSize","collapseUnprocessedSize","simplifyProcessedSize","simplifyLookupSize","simplifyUnprocessedSize","trackerLookupSize")]
#dtbuffers <- dtbuffers[sample(1:nrow(dtbuffers), 100000),]
dtbuffermelt <- melt(dtbuffers)
ggplot(dtbuffermelt, aes(x=variable, y=value)) + geom_boxplot() + scale_y_log10() + labs(x="Data Structure", y="Size") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
saveplot("na12878_internal_buffer_size")
dtops <- dtbuffer[,c("contigSize","contigFrontierSize","contigMemoizedSize","contigUnprocessedSize","assemblyActiveSize")]
#dtops <- dtops[sample(1:nrow(dtops), 100000),]
dtopmelt <- melt(dtops)
ggplot(dtopmelt, aes(x=variable, y=value)) + geom_boxplot() + scale_y_log10() + labs(x="", y="") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
saveplot("na12878_internal_assembly_size")
ggplot(rbind(dtopmelt,dtbuffermelt), aes(x=variable, y=value)) + geom_boxplot() + scale_y_log10() + labs(x="", y="") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
saveplot("na12878_assembly_buffers")

mean(dtbuffer$trackerActive)
mean(dtbuffer$trackerLookupSize)
