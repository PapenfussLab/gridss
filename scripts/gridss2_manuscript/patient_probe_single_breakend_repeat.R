source("libbenchmark.R")
rmdf = import.repeatmasker.fa.out("probe_repeat_masker/breakend_probe_sequences.fasta.out")
rmdf$tp = str_detect(seqnames(rmdf), "TRUE")

ggplot(as.data.frame(rmdf)) +
	aes(x=repeatClass, fill=tp) +
	geom_bar()

rmdf %>% as.data.frame() %>%
	group_by(repeatClass) %>%
	summarise(tp=sum(tp), n=n(), pct=sum(tp)/n()) %>%
	arrange(desc(pct)) %>%
	View()

ggplot(as.data.frame(rmdf)) +
	aes(x=width, fill=tp) +
	facet_wrap(~ repeatClass) +
	geom_histogram()

as.data.frame(rmdf) %>%
	dplyr::select(seqnames, tp) %>%
	distinct() %>%
	summarise(tp=sum(tp), n=n(), pct=sum(tp)/n())


as.data.frame(rmdf) %>%
	group_by(seqnames, tp) %>%
	summarise(has_simple=any(repeatClass=="Simple_repeat")) %>%
	ungroup() %>%
	group_by(has_simple) %>%
	summarise(tp=sum(tp), n=n(), pct=sum(tp)/n())
