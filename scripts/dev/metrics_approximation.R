library(tidyverse)
library(scales)
theme_set(theme_bw())

isdf = read_tsv("insertsize.tsv", col_names=c("n", "insertsize", "count")) %>%
	group_by(n) %>%
	mutate(rate=count/sum(count)) %>%
	ungroup()

flevels=isdf %>% dplyr::select(n) %>% arrange(n) %>% distinct() %>% pull(n)
flabels=comma(flevels)

ggplot(isdf %>% filter(n > 500000)) +
	aes(x=insertsize, y=rate, colour=factor(n, levels=flevels, labels=flabels)) +
	geom_point(size=0.1) +
	labs(title="Insert size distribution vs #reads")

scdf = read_tsv("softclip.tsv", col_names=c("n", "sclen", "count"))%>%
	group_by(n) %>%
	mutate(rate=count/sum(count)) %>%
	ungroup()

ggplot(scdf) +
	aes(x=sclen, y=rate, colour=factor(n, levels=flevels, labels=flabels)) +
	geom_point() + scale_y_log10() +
	labs(title="Soft clip length distribution vs #reads")

idsvdf = read_tsv("idsv.tsv") %>%
	mutate(
		oea_rate=READ_PAIRS_ONE_MAPPED / READ_PAIRS)

ggplot(idsvdf %>% filter(n > 500000)) +
	aes(x=n, y=oea_rate) +
	geom_point()
	labs(title="mapping rate vs #reads")
