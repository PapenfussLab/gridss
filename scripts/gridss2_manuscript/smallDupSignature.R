source("libbenchmark.R")
lnx_svs = read_csv(paste0(privatedatadir, "/hartwig/LNX_SVS.csv"))
duptypedf = read_csv(paste0(privatedatadir, "/hartwig/sample_dup_type.csv"))
rmdf = read_tsv("../../../ref/refdata/hg19/dbs/repeatmasker/hg19.rm.bedops.bed",
								col_type="ciicdcdddccciic",
								col_names=c("chr", "start", "end", "repeatname", "score", "strand", "pctsub", "pctdel", "pctins", "querybases", "repeatclass", "matchstart", "matchend", "uid"))
rmgr = with(rmdf, GRanges(seqnames=chr, ranges=IRanges(start=start, end=end), repeatname=repeatname, repeatclass=repeatclass))
seqlevelsStyle(rmgr) = "NCBI"

truncated_small_dup_label = function(primaryTumorLocation) {
	ifelse(primaryTumorLocation %in% c("Colon/Rectum", "Colon/Rectum", "Skin", "Breast", "Lung", "Nervous system"), primaryTumorLocation, "Other")
}
truncated_small_dup_label_df = duptypedf %>%
	dplyr::select(sampleId, msIndelsPerMb, primaryTumorLocation) %>%
	filter(sampleId %in% lnx_svs$SampleId) %>%
	mutate(label=truncated_small_dup_label(primaryTumorLocation)) %>%
	group_by(label) %>%
	summarise(count=length(unique(sampleId))) %>%
	mutate(label_text=paste0(label, " (", count, ")"))

smldup_svs = lnx_svs %>%
	filter(Type == "DUP" & PosEnd - PosStart < 5000) %>%
	left_join(duptypedf %>% dplyr::select(sampleId, msIndelsPerMb, primaryTumorLocation), by=c("SampleId"="sampleId")) %>%
	mutate(primaryTumorLocation %na% "Unknown")
simplegr = with(smldup_svs, GRanges(seqnames=ChrStart, ranges=IRanges(start=PosStart, end=PosEnd)))
smldup_svs$rmclass = rmgr$repeatclass[findOverlaps(simplegr, rmgr, select="first")]
smldup_svs$inmicrosatellite = overlapsAny(simplegr, rmgr[rmgr$repeatclass == "Simple_repeat"])

ggplot(smldup_svs) +
	aes(x=PosEnd - PosStart, fill=ifelse(inmicrosatellite, "Overlapping microsatellite", "Not in microsatellite")) +
	geom_histogram() +
	geom_vline(xintercept=50) +
	scale_x_log10(limits=c(32, 1000), breaks=c(32, 50, 100, 200, 500, 1000), expand=c(0,0,0,0)) +
	scale_y_continuous(expand=c(0,0,0.05,0)) +
	labs(x="Duplication Size", fill="")
	

smldup_svs %>%
	mutate(label=truncated_small_dup_label(primaryTumorLocation)) %>%
	group_by(SampleId, label) %>%
	summarise(lt100dup=sum(!inmicrosatellite & PosEnd - PosStart < 100)) %>%
	ungroup() %>%
	arrange(desc(lt100dup)) %>%
	mutate(ordinal=row_number()) %>%
	inner_join(truncated_small_dup_label_df) %>%
	ggplot()  +
	aes(x=ordinal, y=lt100dup, fill=label_text) +
	geom_bar(stat="identity") +
	labs(x="", y="Non-microsatellite duplications less than 100bp", fill="") +
	theme(plot.margin = margin(0,0,0,0)) +
	scale_y_continuous(expand=c(0,0,0.05,0)) +
	scale_x_continuous(expand=c(0,0,0,0), breaks=seq(0, 3800, 50))
figsave("small_dup_samples_all", width=100, height=4, limitsize=FALSE)



# PCAWG comparison
#cd publicdata/pcawgcnsv; grep DUP *.bedpe > dups.bedpe
pcawgdups = read_tsv(paste0(datadir, "pcawgcnsv/dups.bedpe"), col_names=c("sample_chr", "start1", "end1", "chr2", "start2", "end2", "id", "score", "strand1", "strand2", "svtype", "method"))
pcawgdups %>%
	mutate(duplen = abs(start2 - start1)) %>%
	filter(duplen < 1000000) %>%
	ggplot() +
	aes(x=duplen) +
	geom_histogram() +
	scale_x_log10() +
	labs(title="PCAWG small dups")
pcawgdups %>%
	mutate(duplen = abs(start2 - start1)) %>%
	mutate(dupbin=cut(duplen, breaks=c(0, 32, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 5000, 10000, 100000, 1000000, 10000000, 100000000))) %>%
	group_by(dupbin) %>%
	summarise(n=n())













