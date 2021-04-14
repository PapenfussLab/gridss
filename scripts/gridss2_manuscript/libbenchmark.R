library(StructuralVariantAnnotation)
library(rtracklayer)
library(tidyverse)
library(ggrepel)
library(cowplot)
library(BSgenome.Hsapiens.UCSC.hg19)
options(stringsAsFactors = FALSE)
datadir = "./publicdata/"
privatedatadir = "./protecteddata/"
figdir = "./figures/"
pon_dir = paste0(datadir, "pon3792v1/")

gridss_fig_tp_colours = c("#6baed6", "#3182bd", "#08519c")
gridss_fig_fp_colours = c("#fb6a4a", "#de2d26", "#a50f15")

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

#' Includes self
lastTrueOrdinal = function(boolvec, missingValue=-1) {
	cummax(ifelse(boolvec, seq_along(boolvec), missingValue))
}
#' Includes self so distance can be zero
distanceToLastTrueValue = function(boolvec, missingValue=-1) {
	lto = lastTrueOrdinal(boolvec, missingValue=missingValue)
	ifelse(lto == missingValue, lto, seq_along(boolvec) - lto)
}
secondLastTrueOrdinal = function(boolvec, missingValue=-1) {
	padded = c(TRUE, TRUE, boolvec)
	lto = pmax(1, lastTrueOrdinal(padded, missingValue))
	slto = lag(lto)[lto]
	ifelse(is.na(slto) | slto - 2 <= 0, missingValue, slto - 2)[c(-1, -2)]
}
testthat::test_that("lastTrueOrdinal", {
	assertthat::assert_that(assertthat::are_equal(lastTrueOrdinal(c(FALSE, TRUE, FALSE, TRUE, TRUE, FALSE)), c(-1, 2, 2, 4, 5, 5)))
	assertthat::assert_that(assertthat::are_equal(secondLastTrueOrdinal(c(FALSE, TRUE, FALSE, TRUE, TRUE, FALSE)), c(-1, -1, -1, 2, 4, 4)))
	assertthat::assert_that(assertthat::are_equal(distanceToLastTrueValue(c(F, T, T, F, F, F, T, F, T)), c(-1, 0, 0, 1, 2, 3, 0, 1, 0)))

})

if (!exists("pon_bedpe_gr")) {
	pon_bedpe_gr=cached_read_file(paste0(pon_dir, "gridss_pon_breakpoint.bedpe"), read_gridss_breakpoint_pon)
}
#pon_bed_gr=cached_read_file(paste(pon_dir, "gridss_pon_single_breakend.bed", sep="/"), import), # not needed - GRIDSS2 has PON FILTER and nobody else calls single breakends
pon_filter = function(gr) {
	return (gr[!gridss_overlaps_breakpoint_pon(gr, pongr=pon_bedpe_gr)])
}
calc_roc_pass_all = function(truth_gr, caller_gr, sample_name, ...) {
	all_calls=calc_roc(truth_gr, caller_gr, ...)
	all_calls$roc$subset = "All calls"
	all_calls$roc_by_type$subset = "All calls"
	all_calls$gr$subset = "All calls"
	all_calls$gr$sample_name = sample_name
	all_calls$roc$sample_name = sample_name
	all_calls$roc_by_type$sample_name = sample_name
	pass_calls=calc_roc(truth_gr, caller_gr[is.na(caller_gr$FILTER) | caller_gr$FILTER %in% c("PASS", "", ".")], ...)
	pass_calls$roc$subset = "PASS only"
	pass_calls$roc_by_type$subset = "PASS only"
	pass_calls$gr$subset = "PASS only"
	pass_calls$gr$sample_name = sample_name
	pass_calls$roc$sample_name = sample_name
	pass_calls$roc_by_type$sample_name = sample_name
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
		caller_gr = caller_gr[caller_gr$breakendInRegion & caller_gr$isIntrachromosomal ]#& caller_gr$fullyContainedInRegion]
		caller_gr = caller_gr[is.na(caller_gr$partner) | caller_gr$partner %in% names(caller_gr)]
	}
	if (!is.null(additional_filter)) {
		truth_gr = additional_filter(truth_gr)
		caller_gr = additional_filter(caller_gr)
	}
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
# Consistency with SV calls:
# segments
# segments supported by SV
# distance to SV
evaluate_cn_transitions = function (cngr, svgr, margin=100000, expectConsistent=TRUE) {
	cn_transitions = with(cngr %>% as.data.frame(), IRanges::reduce(c(
		GRanges(seqnames=seqnames, ranges=IRanges(start=start, width=1)),
		GRanges(seqnames=seqnames, ranges=IRanges(start=end + 1, width=1)))))
	cn_transitions$caller = unique(cngr$caller)
	cn_transitions$cn_left = NA
	cn_transitions$cn_right = NA
	hits = findOverlaps(cn_transitions, with(cngr %>% as.data.frame(), GRanges(seqnames=seqnames, ranges=IRanges(start=start, width=1))))
	cn_transitions$cn_right[queryHits(hits)] = cngr[subjectHits(hits)]$cn
	hits = findOverlaps(cn_transitions, with(cngr %>% as.data.frame(), GRanges(seqnames=seqnames, ranges=IRanges(start=end + 1, width=1))))
	cn_transitions$cn_left[queryHits(hits)] = cngr[subjectHits(hits)]$cn
	cn_transitions$cn_delta = cn_transitions$cn_right - cn_transitions$cn_left
	
	hits = findOverlaps(svgr, cn_transitions, maxgap=margin, ignore.strand=TRUE) %>% as.data.frame() %>%
		mutate(distance=pmax(1, abs(start(svgr)[queryHits] - start(cn_transitions)[subjectHits])))
	best_sv_hit = hits %>% group_by(queryHits) %>% filter(distance==min(distance)) %>% ungroup()
	best_cn_hit = hits %>% group_by(subjectHits) %>% filter(distance==min(distance)) %>% ungroup()
	cn_transitions$distance = NA
	cn_transitions[best_cn_hit$subjectHits]$distance = best_cn_hit$distance
	svgr$distance = rep(NA, length(svgr))
	svgr[best_sv_hit$queryHits]$distance = best_sv_hit$distance
	
	exact_bp_hit = best_sv_hit %>%
		mutate(queryHits.p=match(svgr[queryHits]$partner, names(svgr))) %>%
		filter(!is.na(queryHits.p)) %>%
		inner_join(best_sv_hit, by=c("queryHits.p"="queryHits"), suffix=c("", ".p")) %>%
		filter(distance == 1 & distance.p == 1) %>%
		mutate(
			ori_local = as.character(strand(svgr))[queryHits],
			ori_remote = as.character(strand(svgr))[queryHits.p],
			sv_cn = cn_transitions$cn_delta[subjectHits.p] * ifelse(ori_remote == "+", -1, 1),
			cn_left = cn_transitions$cn_left[subjectHits] + ifelse(ori_local == "+", 0, sv_cn),
			cn_right = cn_transitions$cn_right[subjectHits] + ifelse(ori_local == "+", sv_cn, 0),
			cn_error=abs(cn_left - cn_right))
	sv_hit_count = best_sv_hit %>% group_by(subjectHits) %>% summarise(n=n())
	cn_transitions$sv_matches = 0
	cn_transitions$sv_matches[sv_hit_count$subjectHits] = sv_hit_count$n
	cn_transitions$sv_match_classification = "Missing SV"
	cn_transitions$sv_match_classification[best_sv_hit %>% filter(is.na(svgr$partner[queryHits])) %>% pull(subjectHits)] = "Single Breakend"
	cn_transitions$sv_match_classification[best_sv_hit %>% filter(!is.na(svgr$partner[queryHits])) %>% pull(subjectHits)] = "Breakpoint"
	cn_transitions$sv_match_rescued = FALSE
	if (!is.null(svgr$Recovered)) {
		cn_transitions$sv_match_rescued[best_sv_hit %>% filter(svgr$Recovered[queryHits]) %>% pull(subjectHits)] = TRUE
	}
	cn_transitions$sv_partner_matches = 0
	cn_transitions$sv_partner_matches[exact_bp_hit$subjectHits] = cn_transitions$sv_matches[exact_bp_hit$subjectHits.p]
	cn_transitions$can_evaluate_cn_error = cn_transitions$sv_matches == 1 & cn_transitions$sv_partner_matches == 1
	cn_transitions$cn_error = NA
	cn_transitions$cn_error[exact_bp_hit$subjectHits] = exact_bp_hit$cn_error
	cn_transitions$cn_error[!cn_transitions$can_evaluate_cn_error] = NA
	
	exact_sv_hit = findOverlaps(svgr, cn_transitions, maxgap=1, ignore.strand=TRUE) %>%
		as.data.frame() %>%
		mutate(distance=abs(start(svgr)[queryHits] + ifelse(as.character(strand(svgr)[queryHits]=="-"), -1, 0)) - start(cn_transitions)[subjectHits]) %>%
		group_by(queryHits) %>%
		filter(distance==min(distance))
	svgr$cn_left = rep(NA, length(svgr))
	svgr$cn_right = rep(NA, length(svgr))
	svgr$orphaned = rep(FALSE, length(svgr))
	svgr$cn_left[exact_sv_hit$queryHits] = cn_transitions$cn_left[exact_sv_hit$subjectHits]
	svgr$cn_right[exact_sv_hit$queryHits] = cn_transitions$cn_right[exact_sv_hit$subjectHits]
	svgr$orphaned[exact_sv_hit$queryHits] = FALSE
	if (expectConsistent & any(svgr$orphaned)) {
		browser()
	} else {
		orphan_sv_hits = findOverlaps(svgr, cngr, ignore.strand=TRUE) %>% as.data.frame()
		svgr$cn_left[orphan_sv_hits$queryHits] = ifelse(!svgr$orphaned, svgr$cn_left[orphan_sv_hits$queryHits], cngr$cn[orphan_sv_hits$subjectHits])
		svgr$cn_right = ifelse(!svgr$orphaned, svgr$cn_right, svgr$cn_left)
	}
	return(list(cn_transitions=cn_transitions, sv=svgr))
}
lnx_to_gr <- function(lnx_svs) {
	lnx_svs = lnx_svs %>% replace_na(list(InsertSeq=""))
	grs = GRanges(
		seqnames=lnx_svs$ChrStart,
		ranges=IRanges(start=lnx_svs$PosStart, width=1),
		strand=ifelse(lnx_svs$OrientStart == 1, "+", "-"),
		InsertSeq=lnx_svs$InsertSeq,
		partner=ifelse(lnx_svs$ChrEnd == 0, NA_character_, paste0(lnx_svs$Id, "h")),
		Id=lnx_svs$Id,
		SampleId=lnx_svs$SampleId,
		beid=paste0(lnx_svs$Id, ifelse(is.na(lnx_svs$ChrEnd), "b",  "o")))
	mcols(grs) = bind_cols(as.data.frame(mcols(grs)), lnx_svs %>% dplyr::select(-SampleId, Id))
	names(grs) = grs$beid
	lnx_svs = lnx_svs %>% filter(ChrEnd != 0)
	rc_insert_sequence = lnx_svs$InsertSeq
	rc_insert_sequence[!str_detect(rc_insert_sequence, "[^ACGTN]")] = as.character(reverseComplement(DNAStringSet(rc_insert_sequence[!str_detect(rc_insert_sequence, "[^ACGTN]")])))
	grh = GRanges(
		seqnames=lnx_svs$ChrEnd,
		ranges=IRanges(start=lnx_svs$PosEnd, width=1),
		strand=ifelse(lnx_svs$OrientEnd == 1, "+", "-"),
		insertSequence=ifelse(lnx_svs$OrientStart != lnx_svs$OrientEnd, lnx_svs$InsertSeq, rc_insert_sequence),
		partner=paste0(lnx_svs$Id, "o"),
		Id=lnx_svs$Id,
		SampleId=lnx_svs$SampleId,
		beid=paste0(lnx_svs$Id, "h"))
	mcols(grh) = bind_cols(as.data.frame(mcols(grh)), lnx_svs %>% dplyr::select(-SampleId, Id))
	names(grh)=grh$beid
	return(c(grs, grh))
}

# UCSC table export
hg19_gaps = with(read_tsv(paste0(datadir, "hg19_gap")),
								 GRanges(seqnames=str_replace(chrom, "chr", ""), ranges=IRanges(start=chromStart, end=chromEnd), type=type))
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
hg19_cytobands =  with(read_tsv(
	file=paste0(datadir, "cytoband.txt"),
	col_names=c("chr", "start", "end", "band", "type"),
	col_type= "ciicc"),
	GRanges(seqnames=chr, ranges=IRanges(start=start+1, end=end), band=band, type=type))
seqlevelsStyle(hg19_cytobands) = "NCBI"
hg19_centromeres = hg19_cytobands[hg19_cytobands$type == "acen"]

#### Load PCAWG data ####
if (!exists(".gr_sv_gridss")) {
	.gr_sv_gridss = list()
	.gr_cn_purple = list()
	.gr_sv_pcawg = list()
	.gr_cn_pcawg = list()
}
load_gr_sv_pcawg = function(sampleId, donorId) {
	if (is.null(.gr_sv_pcawg$sampleId)) {
		svfile=paste0(pcawg_dir, "/", sampleId, ".pcawg_consensus_1.6.161116.somatic.sv.bedpe")
		if (!file.exists(svfile)) {
			write(paste("Missing ", svfile), stderr())
			return(NULL)
		}
		sv_bedpe = read_delim(svfile, delim="\t", col_names=TRUE, col_types="cnncnncncccc")
		svgr = with(sv_bedpe, GRanges(
			seqnames=c(chrom1, chrom2),
			ranges=IRanges(start=c(start1 + 1, start2 + 1), end=c(end1, end2)),
			strand=c(strand1, strand2)
		))
		if (length(svgr) == 0) {
			svgr$sampleId=character(0)
			svgr$sourceId=character(0)
			svgr$partner=character(0)
		} else {
			svgr$sampleId=sampleId
			svgr$sourceId=c(paste0(sampleId, sv_bedpe$sv_id, "_o"), paste0(sampleId, sv_bedpe$sv_id, "_h"))
			svgr$partner=c(paste0(sampleId, sv_bedpe$sv_id, "_h"), paste0(sampleId, sv_bedpe$sv_id, "_o"))
		}
		names(svgr) = svgr$sourceId
		.gr_sv_pcawg$sampleId = svgr
	}
	return(.gr_sv_pcawg$sampleId)
}
load_gr_cn_pcawg = function(sampleId, donorId) {
	if (is.null(.gr_cn_pcawg$sampleId)) {
		cnfile=paste0(pcawg_dir, "/", sampleId, ".consensus.20170119.somatic.cna.txt")
		if (!file.exists(cnfile)) {
			write(paste("Missing ", cnfile), stderr())
			return(NULL)
		}
		cndf = read_delim(cnfile, delim="\t", col_names=TRUE, col_types="cnnnnnn", na="NA")
		cngr = with(cndf, GRanges(
			seqnames=chromosome,
			ranges=IRanges(start=start, end=end),
			sampleId=sampleId,
			cn=total_cn,
			cn_major=major_cn,
			cn_minor=minor_cn,
			star=star))
		.gr_cn_pcawg$sampleId = cngr
	}
	return(.gr_cn_pcawg$sampleId)
}
load_gr_sv_gridss = function(sampleId, donorId) {
	if (is.null(.gr_sv_gridss$sampleId)) {
		gridss_vcf = readVcf(paste0(privatedatadir, "/pcawg/", donorId, "T.purple.sv.vcf.gz"))
		gridss_vcf = gridss_vcf[rowRanges(gridss_vcf)$FILTER == "PASS",]
		gridss_gr = c(
			breakpointRanges(gridss_vcf, nominalPosition=TRUE),
			breakendRanges(gridss_vcf, nominalPosition=TRUE))
		if (length(gridss_gr) > 0) {
			names(gridss_gr) = paste0(sampleId, "_", names(gridss_gr))
			gridss_gr$partner = ifelse(is.na(gridss_gr$partner), NA, paste0(sampleId, "_", gridss_gr$partner))
		}
		.gr_sv_gridss$sampleId = gridss_gr
	}
	return(.gr_sv_gridss$sampleId)
}
load_gr_cn_purple = function(sampleId, donorId) {
	if (is.null(.gr_cn_purple$sampleId)) {
		purple_df = read_table2(paste0(privatedatadir, "/pcawg/", donorId, "T.purple.cnv.somatic.tsv"), col_types="ciiddddcccidiidd")
		purple_gr = with(purple_df, GRanges(
			seqnames=chromosome,
			ranges=IRanges(start=start, end=end),
			cn=minorAlleleCopyNumber + majorAlleleCopyNumber,
			cn_major=minorAlleleCopyNumber,
			cn_minor=majorAlleleCopyNumber,
			sampleId=sampleId))
		.gr_cn_purple$sampleId = purple_gr
	}
	return(.gr_cn_purple$sampleId)
}
#### PCAWG ####
pcawg_dir=paste0(datadir, "pcawgcnsv/")
pcawg_evaluate_cn_transitions = function(sampleId) {
	write(paste("Processing ", sampleId), stderr())
	cngr = load_gr_cn_pcawg(sampleId, NULL)
	svgr = load_gr_sv_pcawg(sampleId, NULL)
	# TODO: find sample pairing
	result = evaluate_cn_transitions(cngr, svgr, expectConsistent=FALSE)
	result$sv$sampleId=rep(sampleId, length(result$sv))
	result$cn_transition$sampleId=rep(sampleId, length(result$cn_transition))
	return(result)
}

# from http://github.com/PapenfussLab/sv_benchmark
import.repeatmasker.fa.out <- function(repeatmasker.fa.out) {
	rmdt <- read_table2(repeatmasker.fa.out, col_names=FALSE, skip=3)
	grrm <- GRanges(
		seqnames=rmdt$X5,
		ranges=IRanges(start=rmdt$X6 + 1, end=rmdt$X7),
		strand=ifelse(rmdt$X9=="C", "-", "+"),
		repeatType=rmdt$X10,
		repeatClass=rmdt$X11)
	return(grrm)
}
