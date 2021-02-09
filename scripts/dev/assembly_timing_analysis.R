library(tidyverse)
library(readr)

posdf = bind_rows(lapply(list.files(path="C:/dev/debug_gridss/oom", pattern="positional-.*.csv", full.names=TRUE), function (x) {
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

bindf = posdf %>%
  mutate(
    msElapsedTime = sum(posdf$nsElapsedTime / 1000000),
    bin=round(pmin(supportPosition, aggregatePosition, pathNodePosition, collapsePosition, simplifyPosition, assemblerPosition), -4)) %>%
  group_by(bin) %>%
  summarise(
    msElapsedTime=sum(msElapsedTime),
    supportConsumed=sum(supportConsumed),
    aggregateConsumed=sum(aggregateConsumed),
    pathNodeConsumed=sum(pathNodeConsumed),
    collapseConsumed=sum(collapseConsumed),
    simplifyConsumed=sum(simplifyConsumed),
    trackerConsumed=sum(trackerConsumed),
    max_trackerActive=sum(trackerActive),
    max_supportProcessedSize = max(supportProcessedSize),
    max_aggregateProcessedSize = max(aggregateProcessedSize),
    max_aggregateQueueSize = max(aggregateQueueSize),
    max_aggregateActiveSize = max(aggregateActiveSize),
    max_pathNodeProcessedSize = max(pathNodeProcessedSize),
    max_pathNodeActiveSize = max(pathNodeActiveSize),
    max_pathNodeEdgeLookupSize = max(pathNodeEdgeLookupSize),
    max_pathNodePathLookupSize = max(pathNodePathLookupSize),
    max_collapseProcessedSize = max(collapseProcessedSize),
    max_collapseUnprocessedSize = max(collapseUnprocessedSize),
    max_simplifyProcessedSize = max(simplifyProcessedSize),
    max_simplifyLookupSize = max(simplifyLookupSize),
    max_simplifyUnprocessedSize = max(simplifyUnprocessedSize),
    max_trackerLookupSize = max(trackerLookupSize),
    max_contigFrontierSize = max(contigFrontierSize),
    max_contigMemoizedSize = max(contigMemoizedSize),
    max_assemblyActiveSize = max(assemblyActiveSize),
    max_contigNodeSize = max(contigNodeSize),
    max_contigStartAnchorNodeSize = max(contigStartAnchorNodeSize),
    max_contigEndAnchorNodeSize = max(contigEndAnchorNodeSize),
    max_contigTruncatedNodeSize = max(contigTruncatedNodeSize),
    max_memoizedSize = max(memoizedSize),
    max_memoizedRemovalSize = max(memoizedRemovalSize),
    max_memoizedPathsRemovalSize = max(memoizedPathsRemovalSize),
    max_descendentPathsRemovalSize = max(descendentPathsRemovalSize),
    max_memoizedPathsRestartSize = max(memoizedPathsRestartSize))
ggplot(bindf) + aes(x=bin, y=msElapsedTime) + geom_point(size=0.1)
ggplot(bindf) + aes(x=bin, y=supportConsumed) + geom_point(size=0.1)
ggplot(bindf) + aes(x=bin, y=max_assemblyActiveSize) + geom_point(size=0.1)

ggplot(bindf) + aes(x=max_assemblyActiveSize, y=max_trackerLookupSize, colour=msElapsedTime > 10000000000) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()  

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












