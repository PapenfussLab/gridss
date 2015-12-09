##################
# Processing steps
##################
# - Follow na12878.R processing steps
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

