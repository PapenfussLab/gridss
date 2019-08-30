library(tidyverse)
library(readr)

posdf = bind_rows(lapply(list.files(path="D:/colo829/chr1/bench/visualisation/", pattern="positional-1_.*.csv", full.names=TRUE), function (x) {
#posdf = bind_rows(lapply(list.files(path="D:/colo829/telemetry/", pattern="positional.*.csv", full.names=TRUE), function (x) {
  df = read_csv(x)
  df$file=x
  return(df)
}))
posdf$is_slow = posdf$nsElapsedTime > 1000000
ggplot(posdf) + aes(x=supportProcessedSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=aggregateProcessedSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=aggregateQueueSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=aggregateActiveSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=pathNodeProcessedSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=pathNodeActiveSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=pathNodeEdgeLookupSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=pathNodePathLookupSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=collapseProcessedSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=collapseUnprocessedSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=simplifyProcessedSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=simplifyLookupSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=simplifyUnprocessedSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=trackerLookupSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=contigFrontierSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=contigMemoizedSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=assemblyActiveSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=contigNodeSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=contigStartAnchorNodeSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=contigEndAnchorNodeSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=contigTruncatedNodeSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=memoizedSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=memoizedRemovalSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10() + scale_x_log10()
ggplot(posdf) + aes(x=memoizedPathsRemovalSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10() + scale_x_log10()
ggplot(posdf) + aes(x=descendentPathsRemovalSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10() + scale_x_log10()
ggplot(posdf) + aes(x=memoizedPathsRestartSize, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10() + scale_x_log10()
ggplot(posdf) + aes(x=nsElapsedTime, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=supportConsumed, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=aggregateConsumed, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=pathNodeConsumed, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=collapseConsumed, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=simplifyConsumed, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=assemblerConsumed, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=trackerConsumed, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=trackerActive, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=collapseTraversalCount, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=collapsedBranchCount, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=collapsedLeafCount, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()
ggplot(posdf) + aes(x=simplifiedCount, fill=is_slow) + geom_histogram(bins=200) + scale_y_log10()


ggplot(posdf %>% filter(descendentPathsRemovalSize > 1000)) +
  aes(x=contigStartPosition, y=descendentPathsRemovalSize) +
  geom_point(size=0.1) +
  scale_y_log10()

library(GenomicRanges)
gr = with(posdf %>% filter(descendentPathsRemovalSize > 150, assemblerFirstPosition-assemblerPosition < 1000000000), {
  reduce(GRanges(seqnames="1", IRanges(start=assemblerPosition, end=assemblerFirstPosition), strand="*"))
})
export(with(posdf %>% filter(descendentPathsRemovalSize > 10000, assemblerFirstPosition-assemblerPosition < 1000000000), {
  reduce(GRanges(seqnames="1", IRanges(start=assemblerPosition, end=assemblerFirstPosition), strand="*"))
}), "D:/colo829/chr1/desc10000.bed")



sizedf = posdf[,c("nsElapsedTime", names(posdf)[str_detect(names(posdf), "Size")])]
sizedf = log(sizedf + 1)
sizedf = sizedf %>% mutate(is_slow = nsElapsedTime > 14)
library(GGally)
ggpairs(sizedf[,1:5] %>% sample_n(1000))
               
eventdf = read_csv("D:/colo829/telemetry/out.telemetry.bam.events.csv", col_types="icccdddld", col_names=c("chunk", "direction", "op", "chr", "start", "end", "count", "filtered", "usec"))

ggplot(eventdf) + aes(x=count, y=usec) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10() + facet_wrap(~ op)











  df = read_csv(x)
  df$file=x
  return(df)
}))












