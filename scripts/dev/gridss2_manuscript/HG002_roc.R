library(StructuralVariantAnnotation)
library(cowplot)
library(tidyverse)
library(rtracklayer)
theme_set(theme_bw())


### helper functions
#' Distance to Nearest, vectorised return
dtn = function(qgr, sgr) {
	qgr$dne = NA
	hitdf = distanceToNearest(qgr, sgr, select=c("all"), ignore.strand=TRUE) %>%
		as.data.frame() %>%
		group_by(queryHits) %>%
		arrange(desc(distance)) %>%
		top_n(1) %>%
		ungroup()
	qgr$dne[hitdf$queryHits] = hitdf$distance
	return(qgr$dne)
}
ispassing = function(df) {return(df %>% filter((!isSingleBreakend & FILTER=="PASS") | (isSingleBreakend & BUM > 0 & BSC > 0 & BEbias < 0.5 & QUAL >= 1500)))}
hackyqual = function(df, bemultiplier) {return(df %>% mutate(QUAL=QUAL * ifelse(isSingleBreakend, bemultiplier, 1)), qualmultipler=bemultipler) }
withcallset = function(df) { return(bind_rows(df %>% mutate(callset="All calls"), df %>% ispassing() %>% mutate(callset="PASS"))) }


rmbed = read_tsv("C:/dev/ref/refdata/hg19/dbs/repeatmasker/hg19.rm.bedops.ncbi.bed",
	col_types=c("ciicdc----c----"),
	col_names=c("seqnames", "start", "end", "repeatname", "score", "strand", "pctsub", "pctdel", "pctins", "pastmatch", "repeatclass", "compconsensus", "matchstart", "matchend", "uid", "highermatch"))
rmgr = with(rmbed, GRanges(seqnames=seqnames, ranges=IRanges(start=start, end=end), strand=strand, repeatname=repeatname, repeatclass=repeatclass))
gridss_vcf = readVcf(paste0(datadir, "germline/gridss_2.9.1_mapq10.vcf"))
gr = calc_roc(giab_tier1_gr, cgrs$gridss2, filter_to_region_gr=giab_tier1_regions)$caller_gr
info(gridss_vcf)$tp = NA
info(gridss_vcf[gr$sourceId])$tp = gr$tp

df = bind_cols(
	as.data.frame(rowRanges(gridss_vcf)),
	as.data.frame(info(gridss_vcf))) %>%
	filter(!is.na(tp)) %>%
	mutate(
		isPASS=FILTER=="PASS",
		isASSEMBLY_ONLY=str_detect(FILTER, "ASSEMBLY_ONLY"),
		isSINGLE_ASSEMBLY=str_detect(FILTER, "SINGLE_ASSEMBLY"),
		isLOW_QUAL=str_detect(FILTER, "LOW_QUAL"),
		isSingleBreakend=VF == 0)
dfgr = with(df, GRanges(seqnames=seqnames, ranges=IRanges(start, end), strand=strand))
df$repeatclass = rmgr$repeatclass[findOverlaps(dfgr, rmgr, select="arbitrary", ignore.strand=TRUE)]
df$repeatclass = str_replace(df$repeatClass, "/.*$", "")
df$distanceToRegionEdge = dtn(dfgr, gaps(giab_tier1_regions))
df$distanceToNearestVariant = dtn(dfgr, giab_tier1_gr)
df$distanceToNearestGridssTp = dtn(dfgr, gr[gr$tp])
df = df %>% mutate(
	BEbias=(BASSR + BASRP - BSC - BUM) / (BASSR + BASRP + BSC + BUM),
	BPbias=(ASSR + ASRP - SR - RP) / (ASSR + ASRP + SR + RP))
bpdf = df %>% filter(!isSingleBreakend)
bedf = df %>% filter(isSingleBreakend)

rocqual <- function(df, ...) {
	rocby(df, QUAL, ...)
	# groupingCols = quos(...)
	# df %>% dplyr::select(!!!groupingCols, QUAL, tp) %>%
	# 	group_by(!!!groupingCols, QUAL) %>%
	# 	summarise(tp=sum(tp), ncalls=n()) %>%
	# 	group_by(!!!groupingCols) %>%
	# 	arrange(desc(QUAL)) %>%
	# 	mutate(cumtp=cumsum(tp),cumncalls=cumsum(ncalls)) %>%
	# 	mutate(precision=cumtp/cumncalls) %>%
	# 	ungroup()
}
rocby <- function(df, rocField, ...) {
	groupingCols = quos(...)
	rocCol = enquo(rocField)
	df %>% dplyr::select(!!!groupingCols, !!rocCol, tp) %>%
		group_by(!!!groupingCols, !!rocCol) %>%
		summarise(tp=sum(tp), ncalls=n()) %>%
		group_by(!!!groupingCols) %>%
		arrange(desc(!!rocCol)) %>%
		mutate(cumtp=cumsum(tp),cumncalls=cumsum(ncalls)) %>%
		mutate(precision=cumtp/cumncalls) %>%
		ungroup()
}

# Overall summary
ggplot(bind_rows(
		df %>% withcallset() %>% rocqual(callset),
		rocqual(df %>% filter(FILTER=="PASS")) %>% mutate(callset="rawPASS")
	)) +
	aes(y=precision, x=cumtp, color=callset) +
	geom_line()


ggplot(bind_rows(
		rocqual(bpdf %>% mutate(QUAL=round(QUAL, -1))) %>% mutate(name="QUAL"),
		rocqual(bpdf %>% mutate(QUAL=round(ASQ, -1))) %>% mutate(name="ASQ"),
		rocqual(bpdf %>% mutate(QUAL=round(SRQ, -1))) %>% mutate(name="SRQ"),
		rocqual(bpdf %>% mutate(QUAL=round(RPQ, -1))) %>% mutate(name="RPQ"),
		rocqual(bpdf %>% mutate(QUAL=round(MQ, 0))) %>% mutate(name="MQ"),
		rocqual(bpdf %>% mutate(QUAL=round(MQN, 0))) %>% mutate(name="MQN"),
		rocqual(bpdf %>% mutate(QUAL=round(MQX, 0))) %>% mutate(name="MQX"),
		rocqual(bpdf %>% mutate(QUAL=AS)) %>% mutate(name="AS"),
		rocqual(bpdf %>% mutate(QUAL=StructuralVariantAnnotation:::elementExtract(HOMLEN) %na% 0)) %>% mutate(name="HOMLEN")
		)) +
	aes(y=precision, x=cumtp, label=QUAL) +
	geom_line() +
	geom_text() +
	facet_wrap(~ name) +
	labs(title="Non-QUAL bp ROCs")


testdf = withcallset(bedf)
ggplot(bind_rows(
	rocqual(testdf %>% mutate(QUAL=round(QUAL, -1)), callset) %>% mutate(name="QUAL"),
		rocqual(testdf %>% mutate(QUAL=round(BAQ, -1)), callset) %>% mutate(name="BAQ"),
		rocqual(testdf %>% mutate(QUAL=round(BUMQ, -1)), callset) %>% mutate(name="BUMQ"),
		rocqual(testdf %>% mutate(QUAL=round(BSCQ, -1)), callset) %>% mutate(name="BSCQ"),
		rocqual(testdf %>% mutate(QUAL=round(BUMQ, -1)), callset) %>% mutate(name="BUMQ"),
		rocqual(testdf %>% mutate(QUAL=round(BMQ, 0)), callset) %>% mutate(name="BMQ"),
		rocqual(testdf %>% mutate(QUAL=round(BMQN, 0)), callset) %>% mutate(name="BMQN"),
		rocqual(testdf %>% mutate(QUAL=round(BMQX, 0)), callset) %>% mutate(name="BMQX"),
		rocqual(testdf %>% mutate(QUAL=BA), callset) %>% mutate(name="BA"),
		rocqual(testdf %>% mutate(QUAL=BUM), callset) %>% mutate(name="BUM"),
		rocqual(testdf %>% mutate(QUAL=BSC), callset) %>% mutate(name="BSC"))) +
	aes(y=precision, x=cumtp, label=QUAL) +
	geom_line() +
	facet_grid(callset ~ name, scales="free_x") +
	labs(title="Non-QUAL be ROCs")

# TODO: POTENTIAL SUPP FIGURE
ggplot(rocqual(bedf, repeatclass)) +
	aes(y=precision, x=cumtp) +
	geom_line() +
	geom_line(data=rocqual(bedf %>% ispassing(), repeatclass), color="blue") +
	facet_wrap( ~ repeatclass, scales="free_x") +
	labs(title="Single breakend performance by repeat class")

ggplot(withcallset(bedf)) + 
	aes(x=distanceToRegionEdge + 1, fill=tp) +
	geom_histogram() +
	scale_x_log10() +
	facet_wrap(~ callset, scales="free")

ggplot(withcallset(bedf)) + 
	aes(x=distanceToNearestVariant + 1, fill=tp) +
	geom_histogram() +
	scale_x_log10() +
	facet_wrap(~ callset, scales="free")

ggplot(withcallset(bedf)) + 
	aes(x=distanceToNearestGridssTp + 1, fill=tp) +
	geom_histogram() +
	scale_x_log10() +
	facet_wrap(~ callset, scales="free")

ggplot(rocqual(df, isASSEMBLY_ONLY, isPASS)) +
	aes(y=precision, x=cumtp) +
	geom_line() +
	facet_wrap(isPASS ~ isASSEMBLY_ONLY)

ggplot(df %>% filter(FILTER=="PASS")) +
	aes(x=QUAL, fill=tp, color=isSingleBreakend) +
	geom_histogram() +
	facet_wrap(~ isSingleBreakend) +
	scale_x_log10()

ggplot(rocqual(df %>% filter(FILTER=="PASS")))  +
	aes(y=precision, x=cumtp, label=QUAL) +
	geom_line()

ggplot(rocqual(bedf %>% filter(BUM > 1 & BSC > 1)
	))  +
	aes(y=precision, x=cumtp, label=QUAL) +
	geom_line(data=rocqual(bedf %>% filter(FILTER=="PASS")), color="red") +
	geom_line(data=rocqual(bedf %>% filter((BUM > 0 & BSC > 0))), color="blue") +
	geom_line()

ggplot(bedf %>% filter(QUAL > 1500)) +
	aes(x=BVF / (REF1 + REFPAIR + BVF), fill=tp) +
	geom_histogram()

ggplot(withcallset(bedf)) +
	aes(x=str_length(ALT), fill=tp) +
	facet_wrap(~ callset, scales="free") +
	geom_histogram()

ggplot(withcallset(bedf)) +
	aes(x=acss::entropy(str_replace(ALT, "[.]", "")), fill=tp) +
	facet_wrap(~ callset, scales="free") +
	geom_histogram()

# TODO: check how close the FPs are to the regions of uncertainty
# -ve 
# TODO: check I'm not doing a TPDUP for single breakends
# TODO: check BEALN
# - not written yet
# TODO: check breakend sequence entropy
# -ve


ggplot(rocqual(bedf %>% mutate(BAQpctbin=cut(BAQ/QUAL, seq(0, 1, 0.2))), BAQpctbin)) +
	aes(y=precision, x=cumtp, colour=BAQpctbin) +
	geom_line() +
	scale_colour_brewer(palette="Dark2")

ggplot(bedf) +
	aes(x=BEbias, fill= tp) +
	facet_wrap( ~ FILTER, scales="free") + 
	geom_histogram()

ggplot(bpdf) +
	aes(x=BPbias, fill= tp) +
	facet_wrap( ~ FILTER, scales="free") + 
	geom_histogram()


















