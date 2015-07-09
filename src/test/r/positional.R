library(ggplot2)
library(stringr)
library(reshape2)

#dtraw <- read.csv("W:/na12878/gridss_vis/positional-14.csv")
dtraw <- read.csv("W:/i/data.fastcompare/f5ec72c351cf2d12e54125c8a546ea47/gridss_vis/positional-chr12.csv")

dtraw$index <- seq_len(nrow(dtraw))
for (pos in names(dtraw)[str_detect(names(dtraw), "Position")]) {
  dtraw[,paste0(pos,"Delta")] <- dtraw[,pos] - c(0, head(dtraw[,pos], nrow(dtraw) - 1))
}
dtlong <- melt(dtraw, id.vars="index")
dtSize <- dtlong[str_detect(dtlong$variable, "Size"),]

ggplot(dtSize, aes(x=variable, y=value + 1)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(x="Data Structure Size", y="Record count", title="Size of streaming positional de Bruijn graph data structures")



tail(dtraw$collapsedBranchCount)
tail(dtraw$collapsedLeafCount)
