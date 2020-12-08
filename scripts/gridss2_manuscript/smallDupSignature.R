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
	geom_histogram(bins=150) +
	geom_vline(xintercept=50) +
	scale_x_continuous(limits=c(32, 1000), breaks=c(32, seq(50, 1000, 50)), expand=c(0,0,0,0)) +
	#scale_x_log10(limits=c(32, 1000), breaks=c(32, 50, 100, 200, 500, 1000), expand=c(0,0,0,0)) +
	scale_y_continuous(expand=c(0,0,0.05,0)) +
	scale_fill_manual(values=c("#1b9e77", "#7570b3")) +
	labs(x="Duplication Size", fill="")
figsave("small_dup_signature", width=30, height=3.5)

persamplesmldupdf = smldup_svs %>%
	mutate(label=truncated_small_dup_label(primaryTumorLocation)) %>%
	group_by(SampleId, label) %>%
	summarise(lt100dup=sum(!inmicrosatellite & PosEnd - PosStart < 100)) %>%
	ungroup() %>%
	arrange(desc(lt100dup)) %>%
	mutate(ordinal=row_number()) %>%
	inner_join(truncated_small_dup_label_df)
ggplot(persamplesmldupdf) +
	aes(x=ordinal, y=lt100dup, fill=label_text) +
	geom_bar(stat="identity") +
	labs(x="", y="Non-microsatellite duplications less than 100bp", fill="") +
	theme(plot.margin = margin(0,0,0,0)) +
	scale_y_continuous(expand=c(0,0,0.05,0)) +
	scale_x_continuous(expand=c(0,0,0,0), breaks=seq(0, 3800, 50))
figsave("small_dup_samples_all", width=100, height=4, limitsize=FALSE)

ggplot(persamplesmldupdf) +
	aes(x=label_text, y=lt100dup) +
	coord_cartesian(ylim=c(0,25)) +
	geom_boxplot() +
	geom_violin()


indeldf = read_tsv(paste0(privatedatadir, "/hartwig/indel.tsv"), col_names = c("SampleId", "len", "count"), col_types="cii")
indelbindf = indeldf %>%
	mutate(bin=cut(abs(len), breaks=c(0, 8, 16, 1000))) %>%
	group_by(SampleId, bin) %>%
	summarise(indels=sum(count))

event_count_df =
	persamplesmldupdf %>%
	inner_join(lnx_svs %>% group_by(SampleId) %>% summarise(total_svs=n())) %>%
	inner_join(lnx_svs %>% filter(Type == "DUP") %>% group_by(SampleId) %>% summarise(dups=n())) %>%
	left_join(indeldf %>% filter(len > 0 & abs(len) < 8) %>% group_by(SampleId) %>% summarise(del_0_8=n())) %>%
	left_join(indeldf %>% filter(len > 0 & abs(len) >= 8 & abs(len) < 16) %>% group_by(SampleId) %>% summarise(del_8_16=n())) %>%
	left_join(indeldf %>% filter(len > 0 & abs(len) >= 16) %>% group_by(SampleId) %>% summarise(del_16_=n())) %>%
	left_join(indeldf %>% filter(len < 0 & abs(len) < 8) %>% group_by(SampleId) %>% summarise(ins_0_8=n())) %>%
	left_join(indeldf %>% filter(len < 0 & abs(len) >= 8 & abs(len) < 16) %>% group_by(SampleId) %>% summarise(ins_8_16=n())) %>%
	left_join(indeldf %>% filter(len < 0 & abs(len) >= 16) %>% group_by(SampleId) %>% summarise(ins_16_=n())) %>%
	replace_na(list(ins_0_8=0, ins_8_16=0, ins_16_=0, del_0_8=0, del_8_16=0, del_16_=0))

with(event_count_df, cor(lt100dup, total_svs - lt100dup))
with(event_count_df, cor(lt100dup, dups - lt100dup))
with(event_count_df, cor(lt100dup, del_0_8))
with(event_count_df, cor(lt100dup, del_8_16))
with(event_count_df, cor(lt100dup, del_16_))
with(event_count_df, cor(lt100dup, ins_0_8))
with(event_count_df, cor(lt100dup, ins_8_16))
with(event_count_df, cor(lt100dup, ins_16_))

ggplot(event_count_df) +
	aes(y=lt100dup, x=total_svs - lt100dup) +
	geom_point(size=0.1)

ggplot(event_count_df) +
	aes(y=lt100dup, x=dups - lt100dup) +
	geom_point(size=0.1)

ggplot(event_count_df) +
	aes(y=lt100dup, x=del_16_) +
	geom_jitter(size=0.1)

ggplot(event_count_df) +
	aes(y=lt100dup, x=ins_16_) +
	geom_jitter(size=0.1)


indeldf %>% inner_join(persamplesmldupdf, by="SampleId", suffix=c(".indel", ".lt100dup")) %>%
	mutate(
		type=ifelse(len < 0, "DEL", "INS"),
		lt100dup_gt10=ifelse(lt100dup > 10, "10+ DUPs", "<10 DUPs"),
		clipped_length=pmin(abs(len), 32)) %>%
	ggplot() +
	#aes(x=count.indel, y=lt100dup) +
	#geom_jitter() +
	aes(x=count.indel) +
	geom_histogram(bins=100) +
	facet_grid(type + lt100dup_gt10 ~ clipped_length, scales="free")

ggplot(indeldf %>% mutate(len=sign(len) * pmin(abs(len), 32))) +
	aes(x=count) +
	geom_histogram() +
	facet_wrap(~ len, scales="free") +
	labs(title="Indel length distribution")

indeldf %>%
	mutate(lt100dup_ge20 = SampleId %in% (persamplesmldupdf %>% filter(lt100dup >= 20) %>% pull(SampleId))) %>%
	group_by(len, lt100dup_ge20) %>%
	summarise(count=sum(count)) %>%
	filter(abs(len) < 32, abs(len) > 0) %>%
ggplot() +
	aes(x=len, y=count, fill=lt100dup_ge20) +
	geom_bar(stat="identity") +
	facet_wrap(~ lt100dup_ge20, scales="free") +
	labs(title="Indel length distribution")



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



library(RMySQL)
library(qvalue)
library(R.cache)
library(qvalue)
memoized.fisher.test <- addMemoization(fisher.test)
if (!exists("db")) db = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
driversdf = DBI::dbGetQuery(db, "Select * from driverCatalog where driverlikelihood > 0.5")
cancertypedf = DBI::dbGetQuery(db, "select sampleId, primaryTumorLocation FROM clinical")
cohortdriversdf = driversdf %>% filter(sampleId %in% lnx_svs$SampleId) %>% dplyr::select(sampleId, gene) 
cohortcancertypedf = cancertypedf %>% filter(sampleId %in% lnx_svs$SampleId) 
cohortcancertype_enough_samples = cohortcancertypedf %>%
	group_by(primaryTumorLocation) %>%
	summarise(n=n()) %>%
	filter(n>30) %>%
	pull(primaryTumorLocation)

cohortsmalldupstatus = persamplesmldupdf %>%
	inner_join(cohortcancertypedf, by=c("SampleId"="sampleId")) %>%
	mutate(smalldup10 = lt100dup >= 15)

driver_association = function(cohortstatus, condition) {
	cohortstatus$hasCondition = condition
	tbl = table(cohortstatus$smalldup10, cohortstatus$hasCondition)
	if (nrow(tbl) < 2 || ncol(tbl) < 2) {
		return(data.frame())
	}
	test_obj = memoized.fisher.test(tbl)
	return(data.frame(
		estimate_odds_ratio=test_obj$estimate,
		estimate_odds_ratio_95_lower=test_obj$conf.int[1],
		estimate_odds_ratio_95_upper=test_obj$conf.int[2],
		p_value=test_obj$p.value))
}
gene_driver_association = function(cohortstatus, gene_name) {
	return(driver_association(cohortstatus, cohortstatus$SampleId %in% cohortdriversdf$sampleId[cohortdriversdf$gene == gene_name]) %>%
 		mutate(gene=gene_name))
}
primary_driver_association = function(cohortstatus, primary) {
	return(driver_association(cohortstatus, cohortstatus$primaryTumorLocation == primary) %>%
				 	mutate(primary=primary))
}
primary_gene_driver_association = function(cohortstatus, primary, gene_name) {
	cohortstatus = cohortstatus %>% filter(primaryTumorLocation == primary)
	return(gene_driver_association(cohortstatus, gene_name) %>%
				 	mutate(
				 		primary=primary,
				 		gene=gene_name))
}

pancancer_driver_assoc = bind_rows(lapply(unique(cohortdriversdf$gene), function(driver_gene) gene_driver_association(cohortsmalldupstatus, driver_gene)))
pancancer_driver_assoc$qvalue = qvalue(pmin(pancancer_driver_assoc$p_value, 1))$qvalues

primaryloc_assoc = bind_rows(lapply(unique(cohortsmalldupstatus$primaryTumorLocation), function(primary) primary_driver_association(cohortsmalldupstatus, primary)))
primaryloc_assoc$qvalue = qvalue(pmin(primaryloc_assoc$p_value, 1))$qvalues

primary_gene_driver_association = bind_rows(lapply(cohortcancertype_enough_samples, function(primary) {
	bind_rows(lapply(unique(cohortdriversdf$gene), function(driver_gene) { primary_gene_driver_association(cohortsmalldupstatus, primary, driver_gene)})) } ))
primary_gene_driver_association$qvalue = qvalue(pmin(primary_gene_driver_association$p_value, 1))$qvalues



cohortsmalldupstatus %>%
	inner_join(lnx_svs %>% group_by(SampleId) %>% summarise(totalBreaks=n())) %>%
ggplot() +
	aes(x=totalBreaks, y=lt100dup, colour=smalldup10) +
	geom_jitter()

cohortsmalldupstatus %>%
	ggplot() +
	aes(x=pmin(lt100dup, 20)) +
	geom_histogram(bins=21) +
	facet_grid(primaryTumorLocation ~ .)

View(pancancer_driver_assoc)
View(primaryloc_assoc)
View(primary_gene_driver_association)