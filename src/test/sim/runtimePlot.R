library(ggplot2)
library(plyr)
setwd("W:/i")
#./runtime.sh | tee time.tsv
sw <- read.csv("runtime.tsv", sep="\t")
sw <- sw[sw$caller!="idsv",]

ggplot(ddply(sw, ~ caller + readDepth, summarise, cpu=mean(user) + mean(sys), wall=mean(real))) +
  aes(x=readDepth, group=caller, color=caller) + 
  geom_line(aes(y=cpu)) + 
  geom_point(aes(y=cpu)) + 
  scale_y_log10() + scale_x_log10() +
  #scale_x_continuous() + scale_y_continuous(limits=c(0, 1000)) +
  geom_text(aes(label=caller, color=caller, y=cpu, hjust=1.1)) +
  ylab("CPU time (s)") + 
  theme_bw() +
  ggtitle("Caller Run-time hg19 chr12 simulated variants")
ggsave("runtime-depth.png")

ggplot(ddply(sw[sw$readDepth==100,], ~ caller + readLength, summarise, cpu=mean(user) + mean(sys), real=mean(real))) +
  aes(x=readLength, group=caller, color=caller) + 
  geom_line(aes(y=cpu)) + 
  geom_point(aes(y=cpu)) + 
  scale_y_log10() + 
  geom_text(aes(label=caller, color=caller, y=cpu, hjust=1.1)) +
  ylab("CPU time (s)") + 
  theme_bw() +
  ggtitle("Caller Run-time at 100x")
ggsave("runtime-length.png")


ggplot(ddply(sw[sw$maxResMem > 0,], ~ caller + readDepth, summarise, cpu=max(maxResMem) / 1024 / 4)) +
  aes(x=readDepth, group=caller, color=caller) + 
  geom_line(aes(y=cpu)) + 
  geom_point(aes(y=cpu)) + 
  scale_y_log10() + 
  geom_text(aes(label=caller, color=caller, y=cpu, hjust=1.1)) +
  ylab("Memory Usage (Mb)") + 
  theme_bw() +
  ggtitle("Caller Memory Usage")
ggsave("memory.emf")

