library(StructuralVariantAnnotation)
library(tidyverse)
theme_set(theme_bw())
setwd("../virusbreakend_manuscript/scripts")

#batvi_files = list.files("gen", pattern="final_hits.txt", recursive=TRUE, full.names=TRUE)
getndepth = function(str) {
	matches=str_match(str, "[_/]([0-9]+)[_/]([0-9]+)x")
	return(data.frame(n=as.integer(matches[,2]), depth=as.integer(matches[,3])))
}
batvimsa = read_tsv("../publicdata/sim/batvi_all_final_hits_MSA.txt", col_names=c("lib", "host_chr", "host_pos", "host_ori", "virus_ori", "virus_pos", "MSA", "reads", "sr", "uniquemap", "multimap", "rankonehits", "undocumented1"), col_types="cciccicciiiid")
batvinomsa = read_tsv("../publicdata/sim/batvi_all_final_hits_noMSA.txt", col_names=c("lib", "host_chr", "host_pos", "host_ori", "virus_ori", "virus_pos", "reads", "sr", "uniquemap", "multimap", "rankonehits", "undocumented1", "undocumented2", "undocumented3"), col_types="ccicciciiiiddd")
batvi = bind_rows(batvimsa, batvinomsa)
batvi = bind_cols(batvi, getndepth(batvi$lib))

vbevcf = readVcf("../publicdata/sim/virusbreakend.vcf")
vbevcf = vbevcf[!str_detect(as.character(seqnames(vbevcf)), "adjusted_kraken_taxid")]
vbe = data.frame(
	host_chr=str_match(as.character(seqnames(vbevcf)), ".vcf:(.*)$")[,2],
	host_pos=start(vbevcf),
	virus_chr=str_match(as.character(alt(vbevcf)), "adjusted_kraken_taxid_([0-9]+)_(.*):([0-9]+)")[,3],
	virus_pos=as.integer(str_match(as.character(alt(vbevcf)), "adjusted_kraken_taxid_([0-9]+)_(.*):([0-9]+)")[,4]),
	virus_taxid=str_match(as.character(alt(vbevcf)), "adjusted_kraken_taxid_([0-9]+)_(.*):([0-9]+)")[,2],
	qual=rowRanges(vbevcf)$QUAL
	)
vbe = bind_cols(vbe, getndepth(as.character(seqnames(vbevcf))))

gridssvcf = readVcf("../publicdata/sim/gridss2.vcf")
gridss = data.frame(
	host_chr=str_match(as.character(seqnames(gridssvcf)), ".vcf:(.*)$")[,2],
	host_pos=start(gridssvcf),
	virus_taxid=info(gridssvcf)$INSTAXID,
	qual=rowRanges(gridssvcf)$QUAL
)
gridss = bind_cols(gridss, getndepth(as.character(seqnames(gridssvcf))))


versedf = read_tsv(
		"../publicdata/sim/verse.tsv", #"../publicdata/sim/verse.tsv",
		col_names=c("chr1", "pos1", "ori1", "chr2", "pos2", "ori2", "rpsr", "confidence"),
		col_types=c("cicciccc")) %>%
	mutate(
		file=str_extract(chr1, "^[^:]+"),
		chr1=str_extract(chr1, "[^:]+$"),
		rp=as.integer(str_extract(rpsr, "^[^+]+")),
		sr=as.integer(str_extract(rpsr, "[^+]+$"))) %>%
	mutate(
		host_chr=ifelse(chr1 == "chrVirus", chr2, chr1),
		host_pos=ifelse(chr1 == "chrVirus", pos2, pos1),
		virus_chr=ifelse(chr1 == "chrVirus", chr1, chr2),
		virus_pos=ifelse(chr1 == "chrVirus", pos1, pos2),
		score=ifelse(confidence=="high", 1000000, 0) + rp + sr)
versedf = bind_cols(versedf, getndepth(versedf$file))

vifidf = read_csv("../publicdata/sim/vifi.tsv", col_names=c("file", "host_chr","Min","Max","Split1","Split2"), col_types="ccnnnn")
vifidf = bind_cols(vifidf, getndepth(vifidf$file))


calls_virus_position_not_reported = bind_rows(
		vifidf %>% mutate(host_pos=Split1, ori="+"),
		vifidf %>% mutate(host_pos=Split2, ori="-")) %>%
	mutate(caller="ViFi") %>%
	bind_rows(gridss %>% mutate(caller="GRIDSS2")) %>%
	mutate(
		host_chr_match=host_chr=="chr1",
		host_pos_match=abs((n * 1000000) - host_pos) < 1000000,
		tp=host_chr_match & host_pos_match) %>%
	group_by(caller, n, depth) %>%
	mutate(tp = tp & cumsum(tp) <= 2) %>%
	ungroup() %>%
	mutate(
		fp = !tp,
		hom_tp = FALSE, # Doesn't report viral position so we can't know
		score=0) %>%
	dplyr::select(n, depth, caller, host_chr, host_pos, score, tp, hom_tp, fp)
calls=bind_rows(
		batvi %>% mutate(caller="BATVI", score=as.numeric(reads)),
		vbe %>% mutate(caller="VIRUSBreakend", score=qual),
		versedf %>% mutate(caller="VERSE")) %>%
	replace_na(list(score=0)) %>%
	mutate(
		n=as.numeric(n),
		depth=as.numeric(depth)) %>%
	dplyr::select(n, depth, caller, host_chr, host_pos, virus_pos, score, n) %>%
	mutate(
		host_chr_match=host_chr=="chr1",
		host_pos_match=abs((n * 1000000) - host_pos) < 1000000,
		virus_start_match=abs((n * 8) - virus_pos) < 100,
		virus_end_match=abs((n * 8 + 1000) - virus_pos) < 100) %>%
	group_by(caller, n, depth) %>%
	mutate(
		start_tp=host_chr_match & host_pos_match & virus_start_match,
		end_tp=host_chr_match & host_pos_match & virus_end_match) %>%
	mutate(
		# Only the first call is called a true positive
		start_tp = start_tp & cumsum(start_tp) == 1,
		end_tp = end_tp & cumsum(end_tp) == 1) %>%
	# Assume anything that gets the virus position correct is a homology
	mutate(
		start_hom_tp=sum(start_tp) == 0 & virus_start_match,
		end_hom_tp=sum(end_tp) == 0 & virus_end_match) %>%
	mutate(
		# Only the first call is called a true positive
		start_hom_tp = start_hom_tp & cumsum(start_hom_tp) == 1,
		end_hom_tp = end_hom_tp & cumsum(end_hom_tp) == 1) %>%
	mutate(
		tp=start_tp | end_tp,
		hom_tp = start_hom_tp | end_hom_tp,
		fp=!tp & !hom_tp) %>%
	bind_rows(calls_virus_position_not_reported) %>%
	mutate(caller=factor(caller, levels=c("VIRUSBreakend", "GRIDSS2", "BATVI", "VERSE", "ViFi")))

# sanity check we never have more than 2 true positives
calls %>% group_by(caller, n, depth) %>%
	summarise(tp=sum(tp), hom_tp=sum(hom_tp)) %>%
	group_by(caller, tp, hom_tp) %>%
	filter(tp + hom_tp > 2)

# 2 false negative records for everything then filter based on #TPs
full_calls = expand_grid(n=1:248, depth=c(5, 10, 15, 30, 60), caller=unique(calls$caller), ordinal=c(1,2)) %>%
	mutate(fn=TRUE, tp=FALSE, hom_tp=FALSE, fp=FALSE) %>%
	bind_rows(calls %>% mutate(fn=FALSE)) %>%
	group_by(caller, n, depth) %>%
	filter(!fn | ordinal > sum(tp | hom_tp)) %>%
	mutate(
		classification=factor(
			ifelse(tp, "True positive", ifelse(hom_tp, "Homologous call", ifelse(fn, "False negative", "False positive"))),
			levels=c("False positive",  "False negative", "Homologous call", "True positive")))

ggplot(full_calls) +
	aes(x=as.factor(depth), fill=classification) +
	geom_bar() +
	scale_fill_manual(values=c("#d95f02", "#BBBBBB", "#7570b3", "#1b9e77")) + 
	scale_y_continuous(expand = c(0, 0, 0.01, 0)) + 
	facet_grid( ~ caller) +
	labs(x="Sequencing Depth", fill="", y="Insertions")

callsroc = calls %>%
	group_by(caller, depth, score) %>%
	summarise(tp=sum(tp | hom_tp), fp=sum(fp)) %>%
	group_by(caller, depth) %>%
	arrange(desc(score)) %>%
	mutate(tp=cumsum(tp), fp=cumsum(fp))
ggplot(callsroc) +
	aes(x=fp, y=tp, colour=as.factor(depth), group=as.factor(depth)) +
	geom_line() +
	facet_wrap(~ caller)


# Insert sites unique to VIRUSBreakend
calls %>% group_by(n) %>%
	filter(!any(caller != "VIRUSBreakend")) %>%
	summarise()















