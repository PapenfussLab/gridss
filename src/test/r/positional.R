library(ggplot2)
library(stringr)
library(reshape2)

dtraw <- read.csv("W:/na12878/gridss_vis/positional-13.csv")
#names(dt)
#names(dt)[str_detect(names(dt), "Position")]
#names(dt)[str_detect(names(dt), "Consumed")]
dtlong <- melt(dtraw)
dtSize <- dtlong[str_detect(dtlong$variable, "Size"),]

ggplot(dtSize, aes(x=variable, y=value + 1)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(x="Data Structure Size", y="Record count", title="Size of streaming positional de Bruijn graph data structures")
