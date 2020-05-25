library(StructuralVariantAnnotation)
library(rtracklayer)
library(tidyverse)
library(ggrepel)
library(cowplot)
datadir = "C:/dev/gridss/scripts/dev/gridss2_manuscript/data/"
figdir = "C:/dev/gridss/scripts/dev/gridss2_manuscript/figures/"


figsave = function(figureName, ...) {
	if(!dir.exists(figdir)) { dir.create(figdir) }
	ggsave(paste0(figdir, figureName, ".pdf"), ...)
	ggsave(paste0(figdir, figureName, ".png"), ...)
}
get_legend = function(p) {
	grobs <- ggplotGrob(p)$grobs
	legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
	return(legend)
}
theme_set(theme_bw() + theme(
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank()))



#' Replaces the NA values in a with corresponding values in b
#' @param a,b objects to be tested or coerced.
#' @return The altered object.
'%na%' <- function(a, b) {
	if (is.null(a) || length(a) == 0) return(b)
	if (is.null(b) || length(b) == 0) return(a)
	return(ifelse(is.na(a), b, a))
}

score_for_caller = function(caller_name, vcf) {
	qual = switch(caller_name,
								"crest"=(info(vcf)$right_softclipped_read_count %na% 0) + (info(vcf)$left_softclipped_read_count %na% 0),
								"pindel"=geno(vcf)$AD[,1,2],
								"delly"=(info(vcf)$PE %na% 0) + (info(vcf)$SR %na% 0),
								#"breakdancer_1.4.5"
								#"gridss_1.6.1",
								#"gridss_2.8.1",
								#"hydra_master20160129"=,
								#"manta_1.1.1",
								# svaba
								rowRanges(vcf)$QUAL)
	qual[is.na(qual)] = 0
	return(qual)
}

calc_roc_pass_all = function(truth_gr, caller_gr, ...) {
	all_calls=calc_roc(truth_gr, caller_gr, ...)
	all_calls$roc$subset = "All calls"
	all_calls$roc_by_type$subset = "All calls"
	all_calls$gr$subset = "All calls"
	pass_calls=calc_roc(truth_gr, caller_gr[is.na(caller_gr$FILTER) | caller_gr$FILTER %in% c("PASS", "", ".")], ...)
	pass_calls$roc$subset = "PASS only"
	pass_calls$roc_by_type$subset = "PASS only"
	pass_calls$gr$subset = "PASS only"
	return(list(roc=bind_rows(all_calls$roc, pass_calls$roc), gr=c(all_calls$gr, pass_calls$gr), roc_by_type=bind_rows(all_calls$roc_by_type, pass_calls$roc_by_type)))
}
calc_roc = function(truth_gr, caller_gr, filter_to_region_gr=NULL, bpmaxgap=100, bemaxgap=100, minsize=50, maxsizedifference=25, keepInterchromosomal=FALSE, additional_filter=function(gr) { gr } ) {
	caller_gr$QUAL = caller_gr$QUAL %na% 0
	caller_name=unique(caller_gr$caller)
	write(paste("Processing", caller_name), stderr())
	caller_gr = caller_gr[is.na(caller_gr$partner) | caller_gr$partner %in% names(caller_gr)]
	if (!is.null(filter_to_region_gr)) {
		truth_gr = truth_gr[overlapsAny(truth_gr, filter_to_region_gr, ignore.strand=TRUE)]
		caller_gr$breakendInRegion = overlapsAny(caller_gr, filter_to_region_gr, ignore.strand=TRUE)
		caller_gr$isIntrachromosomal = seqnames(caller_gr) == seqnames(caller_gr[caller_gr$partner %na% names(caller_gr)])
		caller_interval_gr = GRanges(
			seqnames = seqnames(caller_gr),
			ranges = IRanges(
				start=pmin(start(caller_gr), start(caller_gr[caller_gr$partner %na% names(caller_gr)])),
				end=pmax(end(caller_gr), end(caller_gr[caller_gr$partner %na% names(caller_gr)]))))
		caller_gr$fullyContainedInRegion = overlapsAny(caller_interval_gr, filter_to_region_gr, type="within")
		caller_gr = caller_gr[caller_gr$breakendInRegion & caller_gr$isIntrachromosomal & caller_gr$fullyContainedInRegion]
	}
	truth_gr = additional_filter(truth_gr)
	caller_gr = additional_filter(caller_gr)
	truth_gr = truth_gr[is.na(truth_gr$partner) | truth_gr$partner %in% names(truth_gr)]
	caller_bpgr = caller_gr[!is.na(caller_gr$partner)]
	caller_begr = caller_gr[is.na(caller_gr$partner)]
	
	rawbphits = findBreakpointOverlaps(caller_bpgr, truth_gr, maxgap=bpmaxgap)
	rawdihits = findInsDupOverlaps(caller_bpgr, truth_gr, maxgap=bpmaxgap, maxsizedifference=maxsizedifference)
	bphits = data.frame(
		queryHits=c(queryHits(rawbphits), queryHits(rawdihits)),
		subjectHits=c(subjectHits(rawbphits), subjectHits(rawdihits))) %>%
		mutate(QUAL=caller_bpgr$QUAL[queryHits]) %>%
		group_by(subjectHits) %>%
		arrange(desc(QUAL)) %>%
		mutate(isBestHit=row_number() == 1) %>%
		ungroup()
	bestbphits = bphits %>% filter(isBestHit)
	
	behits = findOverlaps(caller_begr, truth_gr, maxgap=bemaxgap) %>%
		as.data.frame() %>%
		mutate(QUAL=caller_begr$QUAL[queryHits]) %>%
		group_by(subjectHits) %>%
		arrange(desc(QUAL)) %>%
		mutate(isBestHit=row_number() == 1) %>%
		ungroup()
	bestbehits = behits %>% filter(isBestHit)
	
	caller_bpgr$tp = FALSE
	caller_bpgr$tpdup = FALSE
	caller_begr$tp = rep(FALSE, length(caller_begr))
	caller_begr$tpdup = rep(FALSE, length(caller_begr))
	truth_gr$tp = FALSE
	truth_gr$QUAL = -1
	truth_gr$betp = FALSE
	truth_gr$beQUAL = -1
	
	# Match to breakpoint hits
	truth_gr$tp[bphits$subjectHits] = TRUE
	truth_gr$QUAL[bestbphits$subjectHits] = bestbphits$QUAL
	caller_bpgr$tpdup[bphits$queryHits] = TRUE
	caller_bpgr$tp[bestbphits$queryHits] = TRUE
	caller_bpgr$tpdup = caller_bpgr$tpdup & !caller_bpgr$tp
	
	caller_begr$tpdup[behits$queryHits] = TRUE
	caller_begr$tp[bestbehits$queryHits] = TRUE
	caller_begr$tpdup = caller_begr$tpdup & !caller_begr$tp
	truth_gr$tp[behits$subjectHits] = TRUE
	truth_gr$QUAL[bestbehits$subjectHits] = pmax(truth_gr$QUAL[bestbehits$subjectHits], bestbehits$QUAL)
	
	# match breakpoint if either breakend is called
	truth_gr$tp = truth_gr$tp | partner(truth_gr)$tp
	truth_gr$QUAL = pmax(truth_gr$QUAL, partner(truth_gr)$QUAL)
	# use max qual
	
	# wait to size filter here as this allows us to match calls close to the min event size
	if (!keepInterchromosomal) {
		truth_gr = truth_gr[abs(truth_gr$svLen) >= minsize]
		caller_bpgr = caller_bpgr[!is.na(caller_bpgr$svLen) & abs(caller_bpgr$svLen) >= minsize]
	} else {
		truth_gr = truth_gr[is.na(truth_gr$svLen) | abs(truth_gr$svLen) >= minsize]
		caller_bpgr = caller_bpgr[is.na(caller_bpgr$svLen) | abs(caller_bpgr$svLen) >= minsize]
	}
	
	truth_gr$fn = !truth_gr$tp
	roc_gr = c(truth_gr, caller_bpgr[!(caller_bpgr$tp | caller_bpgr$tpdup)], caller_begr[!(caller_begr$tp | caller_begr$tpdup)])
	roc_gr$caller = caller_name
	
	roc = roc_gr %>% as.data.frame() %>%
		replace_na(list(fn=FALSE, tpdup=FALSE)) %>%
		mutate(QUAL=round(QUAL)) %>%
		dplyr::select(QUAL, tp, tpdup, fn, caller) %>%
		group_by(caller, QUAL) %>%
		summarise(tp=sum(tp), tpdup=sum(tpdup), fn=sum(fn), ncalls=n()) %>%
		group_by(caller) %>%
		arrange(desc(QUAL)) %>%
		mutate(cumtp=cumsum(tp), cumtpdup=cumsum(tpdup), fn=sum(fn),cumncalls=cumsum(ncalls)) %>%
		mutate(
			recall=cumtp/(max(cumtp) + fn),
			precision=cumtp/cumncalls) %>%
		ungroup() %>%
		filter(QUAL >= 0)
	
	roc_gr$simpletype = simpleEventType(roc_gr)
	roc_by_type = roc_gr %>% as.data.frame() %>%
		replace_na(list(fn=FALSE, tpdup=FALSE)) %>%
		mutate(QUAL=round(QUAL)) %>%
		dplyr::select(QUAL, tp, tpdup, fn, caller, simpletype) %>%
		group_by(caller, simpletype, QUAL) %>%
		summarise(tp=sum(tp), tpdup=sum(tpdup), fn=sum(fn), ncalls=n()) %>%
		group_by(caller, simpletype) %>%
		arrange(desc(QUAL)) %>%
		mutate(cumtp=cumsum(tp), cumtpdup=cumsum(tpdup), fn=sum(fn),cumncalls=cumsum(ncalls)) %>%
		mutate(
			recall=cumtp/(max(cumtp) + fn),
			precision=cumtp/cumncalls) %>%
		ungroup() %>%
		filter(QUAL >= 0)
	return(list(roc=roc, roc_by_type=roc_by_type, gr=roc_gr, caller_gr=c(caller_bpgr, caller_begr)))
}