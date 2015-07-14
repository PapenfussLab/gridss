library(ggplot2)
library(stringr)
library(reshape2)

dtraw <- read.csv("W:/na12878/gridss_vis/positional-4.csv")
#dtraw <- read.csv("W:/i/data.fastcompare/f5ec72c351cf2d12e54125c8a546ea47/gridss_vis/positional-chr12.csv")

dtraw$index <- seq_len(nrow(dtraw))
for (pos in names(dtraw)[str_detect(names(dtraw), "Position")]) {
  dtraw[,paste0(pos,"Delta")] <- dtraw[,pos] - c(0, head(dtraw[,pos], nrow(dtraw) - 1))
}
dtraw$totalSize = rowSums(dtraw[,str_detect(names(dtraw), "Size")])
dtraw$aggregateLag = dtraw$supportPosition - dtraw$aggregatePosition
dtraw$pathNodeLag = dtraw$aggregatePosition - dtraw$pathNodePosition
dtraw$collapseLag = dtraw$pathNodePosition - dtraw$collapsePosition
dtraw$simplifyLag = dtraw$collapsePosition - dtraw$simplifyPosition
dtraw$assemblerLag = dtraw$simplifyPosition - dtraw$assemblerPosition
dtraw$firstContigLag = dtraw$firstContigPosition - dtraw$assemblerPosition
dtraw$calledContigLag = dtraw$calledContigPosition - dtraw$assemblerPosition
dtraw$totalLag = dtraw$supportPosition - dtraw$assemblerPosition

dtlong <- melt(dtraw, id.vars="index")
dtSize <- dtlong[str_detect(dtlong$variable, "Size"),]
dtSize$variable <- str_replace(dtSize$variable, "Size", "")
dtLag <- dtlong[str_detect(dtlong$variable, "Lag"),]
dtLag$variable <- str_replace(dtLag$variable, "Lag", "")
dtLag <- dtLag[abs(dtLag$value) < 1000000,] #remove MAX_VALUE outliers
dtLag$variable <- factor(dtLag$variable, levels=unique(dtLag$variable)) # sort by column order

ggplot(dtSize, aes(x=variable, y=value + 1)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(x="Data Structure Size", y="Record count", title="Size of streaming positional de Bruijn graph data structures") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("bufferSize.png", width=10, height=7.5)

ggplot(dtLag, aes(x=variable, y=value)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(x="Process", y="Base Pairs", title="Genomic width of streaming positional de Bruijn graph data transformation buffers") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("bufferWidth.png", width=10, height=7.5)


tail(dtraw$collapsedBranchCount)
tail(dtraw$collapsedLeafCount)
