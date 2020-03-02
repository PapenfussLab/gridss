# quieten all the dependency package conflict spam
library(BiocGenerics, quietly=TRUE, warn.conflicts=FALSE)
library(S4Vectors, quietly=TRUE, warn.conflicts=FALSE)
library(IRanges, quietly=TRUE, warn.conflicts=FALSE)
library(matrixStats, quietly=TRUE, warn.conflicts=FALSE)
library(DelayedArray, quietly=TRUE, warn.conflicts=FALSE)
library(XVector, quietly=TRUE, warn.conflicts=FALSE)
library(Biostrings, quietly=TRUE, warn.conflicts=FALSE)
suppressPackageStartupMessages(library(Biobase))
# Ok, now load direct package dependencies. It's a pity warn.conflicts does not get passed through
library(VariantAnnotation, quietly=TRUE, warn.conflicts=FALSE)
library(rtracklayer, quietly=TRUE, warn.conflicts=FALSE)
library(StructuralVariantAnnotation, quietly=TRUE, warn.conflicts=FALSE)
options(tidyverse.quiet = TRUE)
library(tidyverse, quietly=TRUE, warn.conflicts=FALSE)
library(stringr, quietly=TRUE, warn.conflicts=FALSE)
library(testthat, quietly=TRUE, warn.conflicts=FALSE)
library(stringdist, quietly=TRUE, warn.conflicts=FALSE)

options(stringsAsFactors=FALSE)
source("gridss.config.R")

#' Replaces the NA values in a with corresponding values in b
#' @param a,b objects to be tested or coerced.
#' @return The altered object.
'%na%' <- function(a, b) {
  if (is.null(a) || length(a) == 0) return(b)
  if (is.null(b) || length(b) == 0) return(a)
  return(ifelse(is.na(a), b, a))
}

#' sum of genotype fields
.genosum <- function(genotypeField, columns) {
	rowSums(genotypeField[,columns, drop=FALSE])
}
.addFilter = function(existing, filterName, appliesTo) {
  existing[appliesTo] = paste(existing[appliesTo], filterName, sep=";")
  return(existing)
}
addVCFHeaders = function(vcf) {
  info(header(vcf)) = unique(as(rbind(as.data.frame(info(header(vcf))), data.frame(
    row.names=c("TAF", "LOCAL_LINKED_BY", "REMOTE_LINKED_BY"),
    Number=c(".", ".", "."),
    Type=c("Float", "String", "String"),
    Description=c(
      "Overall unadjusted allele fraction of tumour samples (not weighted by sample depth, nor purity adjusted)",
      "Breakend linking information",
      "Partner breakend linking information"))),
    "DataFrame"))
  VariantAnnotation::fixed(header(vcf))$FILTER = unique(as(rbind(as.data.frame(VariantAnnotation::fixed(header(vcf))$FILTER), data.frame(
    row.names=c(
      "PON",
      "imprecise",
      "strand_bias",
      "homlen",
      "ihomlen",
      "largeNoRP",
      "smallNoSR",
      "small.del.ligation.fp",
      "small.inv.hom.fp",
      "small.replacement.fp",
      "normalSupport",
      "SRNormalSupport",
      "normalCoverage",
      "af",
      "NoAssembledRP",
      "LongPolyC",
      "minRead",
      "noAssembly",
      "cohortMinSize"),
    Description=c(
      "Found in panel of normals",
      "Imprecise variant",
      "Short event with excessive strand bias in split reads/soft clipped reads overlapping breakpoint",
      "Breakpoint homology length too long (gridss.max_homology_length)",
      "Inexact breakpoint homology length too long (gridss.max_inexact_homology_length)",
      "Large event not supported by any read pairs either directly or via assembly",
      "Short event not supported by any split reads either directly or via assembly",
      "Short deletion that appears to be a ligation artefact",
      "Short inversion with significant sequence homology",
      "Deletion with insertion of the same length that is not a simple inversion.",
      "Too many supporting reads from the normal sample(s) (gridss.allowable_normal_contamination)",
      "Short event with split reads support in the normal sample",
      "Insufficient normal coverage to determine somatic status (gridss.min_normal_depth)",
      "Variant allele fraction too low (gridss.min_af)",
      "Single breakend with no assembled read pairs",
      "Single breakend containing long polyC or polyG run. Likely to be an NovaSeq artefact.",
      "Too few reads directly support the variant (gridss.min_direct_read_support)",
      "no assembly support",
      "Variant is smaller than the minimum event size considered for this cohort"))), "DataFrame"))
  return(vcf)
}


gridss_overlaps_breakpoint_pon = function(gr,
    pon_dir=NULL,
    pongr=cached_read_file(paste(pon_dir, "gridss_pon_breakpoint.bedpe", sep="/"), read_gridss_breakpoint_pon),
    ...) {
  hasHit = rep(FALSE, length(gr))
  if (!is.null(pongr)) {
    hasHit[as.data.frame(findBreakpointOverlaps(gr, pongr[pongr$score >= gridss.pon.min_samples], sizemargin=NULL, restrictMarginToSizeMultiple=NULL, ...))$queryHits] = TRUE
  }
  return(hasHit)
}
gridss_overlaps_breakend_pon = function(gr,
    pon_dir=NULL,
    pongr=cached_read_file(paste(pon_dir, "gridss_pon_single_breakend.bed", sep="/"), import),
    ...) {
  hasHit = rep(FALSE, length(gr))
  if (!is.null(pongr)) {
    hasHit[queryHits(findOverlaps(gr, pongr[pongr$score >= gridss.pon.min_samples], ...))] = TRUE
  }
  return(hasHit)
}
#' For each GRanges breakend, indicates whether the variant
#' should be filtered
#' @param somatic_filters apply somatic filters.
#' Assumes the normal and tumour samples are the first and second respectively
gridss_breakpoint_filter = function(gr, vcf, bsgenome, min_support_filters=TRUE, somatic_filters=TRUE, support_quality_filters=TRUE, normalOrdinal, tumourOrdinal, pon_dir=NULL, pon_maxgap=4) {
	vcf = vcf[names(gr)]
	i = info(vcf)
	g = geno(vcf)
	isShort = is_short_deldup(gr)
	ihomlen = gridss_inexact_homology_length(gr, vcf)
	homlen = elementExtract(info(vcf)$HOMLEN, 1)
	homlen[is.na(homlen)] = 0
	filtered = rep("", length(gr))

	if (!is.null(pon_dir)) {
	  filtered = .addFilter(filtered, "PON", gridss_overlaps_breakpoint_pon(gr, pon_dir, maxgap=pon_maxgap))
	}

	if (support_quality_filters) {
	  # TODO: update this to a binomial test so we don't filter low confidence
	  # variants that are strand biased by chance
	  # strand bias can be NA as it is calculated from SR/SC support only
	  strandbias = pmax(i$SB, 1 - i$SB)
	  filtered = .addFilter(filtered, "strand_bias", !is.na(strandbias) & isShort & strandbias > gridss.max_allowable_short_event_strand_bias)

	  #filtered = .addFilter(filtered, "FlankingHighQualIndel", str_detect(gridss_gr$FILTER, "SINGLE_ASSEMBLY") & (i$RP + i$SR) / i%VF < 0.1 & i$RP >= 2 & i$SR >= 2 # fixed in GRIDSSv1.8.0
	}
	if (min_support_filters) {
	  filtered = .addFilter(filtered, "af", gridss_bp_af(gr, vcf, tumourOrdinal) < gridss.min_af)
	  # Multiple biopsy concordance indicates that assemblies with very few supporting reads are sus
	  filtered = .addFilter(filtered, "minRead", .genosum(g$RP,c(normalOrdinal, tumourOrdinal)) + .genosum(g$SR, c(normalOrdinal, tumourOrdinal)) < gridss.min_direct_read_support)
	  #filtered = .addFilter(filtered, "unanchored", i$ASSR + i$SR + i$IC == 0) # replaced by 'imprecise'
	  # very high coverage hits assembly threshold; we also need to keep transitive calls so we reallocate them to get the correct VF
	  filtered = .addFilter(filtered, "noAssembly", str_detect(gr$FILTER, "NO_ASSEMBLY") & i$VF < 100)
	  # homology FPs handled by normal and/or PON
	  filtered = .addFilter(filtered, "homlen", homlen > gridss.max_homology_length)
	  filtered = .addFilter(filtered, "ihomlen", ihomlen > gridss.max_inexact_homology_length & !(is_short_dup(gr)))
		# Added ASRP into filter otherwise breakpoint chains don't get called
	  filtered = .addFilter(filtered, "largeNoRP", !isShort & .genosum(g$RP,c(normalOrdinal, tumourOrdinal)) + .genosum(g$ASRP,c(normalOrdinal, tumourOrdinal)) == 0)
	  filtered = .addFilter(filtered, "smallNoSR", isShort & .genosum(g$SR,c(normalOrdinal, tumourOrdinal)) == 0)

	  filtered = .addFilter(filtered, "small.del.ligation.fp", is_likely_library_prep_fragment_ligation_artefact(gr, vcf))
	  filtered = .addFilter(filtered, "small.inv.hom.fp", is_small_inversion_with_homology(gr, vcf))
	  if (!is.null(bsgenome)) {
	    filtered = .addFilter(filtered, "small.replacement.fp", is_indel_artefact(gr, bsgenome))
	  }
	  filtered = .addFilter(filtered, "cohortMinSize", is_too_small_event(gr))
	}
	if (somatic_filters) {
		#normalaf <- gridss_af(gr, vcf, normalOrdinal)
	  filtered = .addFilter(filtered, "normalSupport", .genosum(g$VF,normalOrdinal) > gridss.allowable_normal_contamination * .genosum(g$VF,tumourOrdinal))
	  filtered = .addFilter(filtered, "SRNormalSupport", isShort & .genosum(g$SR, normalOrdinal) != 0)
	  filtered = .addFilter(filtered, "normalCoverage", .genosum(g$REF, normalOrdinal) + .genosum(g$REFPAIR, normalOrdinal) + .genosum(g$VF, normalOrdinal) < gridss.min_normal_depth)
	}
	return(filtered)
}
#' For each GRanges breakend, indicates whether the variant
#' should be filtered
#' @param somatic_filters apply somatic filters.
#' Assumes the normal and tumour samples are the first and second respectively
gridss_breakend_filter = function(gr, vcf, min_support_filters=TRUE, somatic_filters=TRUE, normalOrdinal, tumourOrdinal, pon_dir=NULL) {
  vcf = vcf[names(gr)]
  i = info(vcf)
  g = geno(vcf)
  filtered = rep("", length(gr))
  if (!is.null(pon_dir)) {
    filtered = .addFilter(filtered, "PON", gridss_overlaps_breakend_pon(gr, pon_dir))
  }
  if (min_support_filters) {
    filtered = .addFilter(filtered, "af", gridss_be_af(gr, vcf, tumourOrdinal) < gridss.min_af)
    filtered = .addFilter(filtered, "imprecise", i$IMPRECISE)
    filtered = .addFilter(filtered, "NO_ASSEMBLY", str_detect(gr$FILTER, "NO_ASSEMBLY"))
    # require direct read support as well
    #filtered = .addFilter(filtered, "minRead", .genosum(g$BSC,c(normalOrdinal, tumourOrdinal)) + .genosum(g$BUM, c(normalOrdinal, tumourOrdinal)) < gridss.min_direct_read_support * gridss.single_breakend_multiplier)
    # require at least one read pair included in the assembly
    # this is a relatively strict filter but does filter out most of the
    # noise from microsatellite sequences
    filtered = .addFilter(filtered, "NoAssembledRP", i$BASRP == 0)
    filtered = .addFilter(filtered, "LongPolyC", str_detect(gr$insSeq, "CCCCCCCCCCCCCCCC") | str_detect(gr$insSeq, "GGGGGGGGGGGGGGGG"))
  }
  if (somatic_filters) {
    filtered = .addFilter(filtered, "normalSupport", .genosum(g$BVF,normalOrdinal) > gridss.allowable_normal_contamination * .genosum(g$BVF,tumourOrdinal))
    filtered = .addFilter(filtered, "normalCoverage", .genosum(g$REF, normalOrdinal) + .genosum(g$REFPAIR, normalOrdinal) + .genosum(g$BVF, normalOrdinal) < gridss.min_normal_depth)
  }
  return(filtered)
}
gridss_inexact_homology_length = function(gr, vcf) {
  if (length(gr) == 0) {
    return(integer(0))
  }
  ihommat = as.matrix(info(vcf[names(gr)])$IHOMPOS)
  ihomlen = pmax(0, ihommat[,2] - ihommat[,1])
  ihomlen[is.na(ihomlen)] = 0
  return(ihomlen)
}
is_short_deldup = function(gr) {
  is_deldup = rep(FALSE, length(gr))
  if (!is.null(gr$partner)) {
    isbp <- gr$partner %in% names(gr)
    bpgr <- gr[isbp]
    bp_short_deldup = strand(bpgr) != strand(partner(bpgr)) &
      seqnames(bpgr) == seqnames(partner(bpgr)) &
      abs(start(bpgr)-start(partner(bpgr))) < gridss.short_event_size_threshold
    is_deldup[isbp] = bp_short_deldup
  }
  return(is_deldup)
}
is_short_dup = function(gr) {
  is_dup = rep(FALSE, length(gr))
  if (!is.null(gr$partner)) {
    isbp <- gr$partner %in% names(gr)
    bpgr <- gr[isbp]
    bp_short_dup = strand(bpgr) != strand(partner(bpgr)) &
      seqnames(bpgr) == seqnames(partner(bpgr)) &
      abs(start(bpgr) - start(partner(bpgr))) < gridss.short_event_size_threshold &
      ((start(bpgr) <= start(partner(bpgr)) & strand(bpgr) == "-") |
        (start(bpgr) >= start(partner(bpgr)) & strand(bpgr) == "+"))
    is_dup[isbp] = bp_short_dup
  }
  return(is_dup)
}
is_short_del = function(gr) {
  return(is_short_deldup(gr) & !is_short_dup(gr))
}
is_likely_library_prep_fragment_ligation_artefact = function(gr, vcf, minsize=100, maxsize=800, minihomlen=6) {
  result = rep(FALSE, length(gr))
  if (!is.null(gr$partner)) {
    isbp <- gr$partner %in% names(gr)
    bpgr = gr[isbp]
    svlen = abs(start(bpgr) - start(partner(bpgr)))
    ihomlen = gridss_inexact_homology_length(gr, vcf)
    result[isbp] = simpleEventType(bpgr) == "DEL" &
      svlen >= minsize & svlen <= maxsize &
      #maxaf <= 0.25
      ihomlen >= minihomlen
  }
  return(result)
}
is_indel_artefact = function(gr, bsgenome, minsizedelta=5, minEditDistancePerBase=0.5, maxEditDistancePerInversionBase=0.2) {
  result = rep(FALSE, length(gr))
  gr$isOfInterest = is_short_del(gr) & abs(abs(start(gr) - start(gr[ifelse(is.na(gr$partner), names(gr), gr$partner)])) - gr$insLen) < minsizedelta
  gr$isOfInterest = gr$isOfInterest & !is.na(gr$partner) & gr[ifelse(is.na(gr$partner), names(gr), gr$partner)]$isOfInterest
  ucscgr = gr
  seqlevelsStyle(ucscgr) = "UCSC"
  gr$isOfInterest = gr$isOfInterest & as.logical(seqnames(ucscgr) %in% seqnames(bsgenome))
  igr = gr[gr$isOfInterest]
  seqlevelsStyle(igr) = "UCSC"
  inseq = igr$insSeq
  igr=GRanges(seqnames=seqnames(igr), ranges=IRanges(start=pmin(start(igr), start(partner(igr))), end=pmax(end(igr), end(partner(igr)))))
  refseq = getSeq(bsgenome, names=igr, as.character=TRUE)
  revSeq = as.character(reverseComplement(DNAStringSet(refseq)))
  fwdEditDistance = stringdist(inseq, refseq, method="lv")
  invEditDistance = stringdist(inseq, revSeq, method="lv")
  fwdEditDistancePerBase = fwdEditDistance / ifelse(nchar(inseq) == 0, 1, nchar(inseq))
  invEditDistancePerBase = invEditDistance / ifelse(nchar(inseq) == 0, 1, nchar(inseq))
  isActualInversion = fwdEditDistancePerBase > minEditDistancePerBase & invEditDistancePerBase < maxEditDistancePerInversionBase
  result[gr$isOfInterest] = !isActualInversion
  return(result)
}
is_small_inversion_with_homology = function(gr, vcf, minhomlen=6, maxsize=40) {
  result = rep(FALSE, length(gr))
  if (!is.null(gr$partner)) {
    isbp <- gr$partner %in% names(gr)
    bpgr = gr[isbp]
    homlen = (as.integer(info(vcf[bpgr$sourceId])$HOMLEN) %na% 0)[isbp]
    #ihomlen = gridss_inexact_homology_length(gr, vcf)
    svlen = abs(start(bpgr) - start(partner(bpgr)))
    result[isbp] = simpleEventType(bpgr) == "INV" &
      svlen <= maxsize &
      homlen >= minhomlen
  }
  return(result)
}
is_too_small_event = function(gr, minSize=gridss.min_event_size) {
  result = rep(FALSE, length(gr))
  if (!is.null(gr$partner)) {
    isbp <- gr$partner %in% names(gr)
    bpgr = gr[isbp]
    svlen = abs(start(bpgr) - start(partner(bpgr))) + str_length(bpgr$insSeq) + ifelse(simpleEventType(bpgr) %in% c("DEL", "INS"), -1, 1)
    result[isbp] = simpleEventType(bpgr) %in% c("DEL", "DUP", "INS") & svlen < minSize
  }
  return(result)
}
#' @description filter out 'shadow' calls of strong multi-mapping calls
#' bwa overestimates the MAPQ of some multimapping reads
#' and makes FP calls to alternate locations.
is_shadow_breakpoint = function(bpgr, begr, vcf, breakendQualMultiple=3) {
  combinedgr <- bpgr
  combinedgr$partner <- NULL
  combinedgr <- c(combinedgr, begr)
  # require exact overlaps
  bestOverlap = findOverlaps(bpgr, combinedgr) %>%
    as.data.frame() %>%
    # bpgr offsets are the same as combinedgr ofsets
    filter(queryHits != subjectHits) %>%
    mutate(beQUAL = combinedgr$QUAL[subjectHits]) %>%
    group_by(queryHits) %>%
    summarise(beQUAL = max(beQUAL))
  bpgr$overlapQUAL = 0
  bpgr$overlapQUAL[bestOverlap$queryHits] = bestOverlap$beQUAL
  i <- info(vcf[bpgr$sourceId])
  better_call_filter = !is_short_event(bpgr) &
    bpgr$overlapQUAL > gridss.shadow_breakend_multiple * bpgr$QUAL &
    i$BANSRQ + i$BANRPQ > i$ASQ
  return(as.logical(better_call_filter))
}

is_short_event = function(gr) {
  seqnames(gr) == seqnames(partner(gr)) &
    abs(start(gr)-start(partner(gr))) < gridss.short_event_size_threshold
}
gridss_bp_af = function(gr, vcf, ordinal) {
  return(.gridss_af(gr, vcf, ordinal, !is_short_deldup(gr), includeBreakpointSupport=TRUE, includeBreakendSupport=FALSE))
}
gridss_be_af = function(gr, vcf, ordinal) {
  return(.gridss_af(gr, vcf, ordinal, includeRefPair=rep(TRUE, length(gr)), includeBreakpointSupport=FALSE, includeBreakendSupport=TRUE))
}
.gridss_af = function(gr, vcf, ordinal, includeRefPair, no_coverage_af=0, includeBreakpointSupport=TRUE, includeBreakendSupport=FALSE) {
  assertthat::are_equal(length(gr), length(includeRefPair))
  genotype = geno(vcf[names(gr)])
  g = lapply(names(genotype), function(field) { if (is.numeric(genotype[[field]])) { .genosum(genotype[[field]], ordinal) } else { genotype[[field]] } })
  names(g) <- names(genotype)
  support = rep(0, length(gr))
  if (includeBreakpointSupport) {
    support = support + g$VF
  }
  if (includeBreakendSupport) {
    support = support + g$BVF
  }
  vf_af = support / (support + g$REF + ifelse(!includeRefPair, 0, g$REFPAIR))
  vf_af[is.nan(vf_af)] = no_coverage_af
  return(vf_af)
}
annotate_overlaps = function(query, subject, ..., group_by_col_name="sampleId") {
	hits = as.data.frame(findBreakpointOverlaps(query, subject, ...))
	if (group_by_col_name %in% names(mcols(query)) & group_by_col_name %in% names(mcols(subject))) {
		hits = hits %>% filter(mcols(query)[[group_by_col_name]][queryHits] == mcols(subject)[[group_by_col_name]][subjectHits])
	}
	hits = hits %>%
		mutate(queryQUAL = query$QUAL[queryHits]) %>%
		mutate(subjectQUAL = subject$QUAL[subjectHits]) %>%
		mutate(QUAL=ifelse(is.na(queryQUAL), subjectQUAL, ifelse(is.na(subjectQUAL), queryQUAL, pmax(queryQUAL, subjectQUAL)))) %>%
		replace_na(list(QUAL=-1)) %>%
		arrange(desc(QUAL)) %>%
		# only allow 1-1 mappings
		filter(!duplicated(queryHits)) %>%
		filter(!duplicated(subjectHits))
	result = rep(NA_character_, length(query))
	result[hits$queryHits] = names(subject)[hits$subjectHits]
	return(result)
}
gr_join_to_df = function(gr1, vcf1, gr2, vcf2, ..., suffix=c(".1", ".2"), group_by_col_name="sampleId") {
	gr1$key = paste0(names(gr1), "_", annotate_overlaps(gr1, gr2, ..., group_by_col_name=group_by_col_name))
	gr2$key = paste0(annotate_overlaps(gr2, gr1, ..., group_by_col_name=group_by_col_name), "_", names(gr2))
	if (group_by_col_name %in% names(mcols(gr1)) & group_by_col_name %in% names(mcols(gr2))) {
		gr1$key = paste(mcols(gr1)[[group_by_col_name]], gr1$key)
		gr2$key = paste(mcols(gr2)[[group_by_col_name]], gr2$key)
	}
	df1 = as.data.frame(gr1)
	df2 = as.data.frame(gr2)
	df1$orientation = paste0(strand(gr1), strand(partner(gr1)))
	df2$orientation = paste0(strand(gr2), strand(partner(gr2)))
	if (!is.null(vcf1)) {
		df1 <- df1 %>% bind_cols(as.data.frame(info(vcf1[gr1$sourceId])))
	}
	if (!is.null(vcf2)) {
		df2 <- df2 %>% bind_cols(as.data.frame(info(vcf2[gr2$sourceId])))
	}
	merged_df = full_join(df1, df2, by="key", suffix=suffix)
	return(merged_df)
}
query_structural_variants_samples = function(dbConnect) {
	query = paste(
		"SELECT DISTINCT sampleId",
		"FROM structuralVariant",
		sep = " ")
	df = dbGetQuery(dbConnect, query)
	return(df$sampleId)
}
simpleEventType <- function(gr) {
  isbp = !is.na(gr$partner)
  result = rep("BE", length(gr))
  bpgr = gr[isbp]
  same_chr = as.character(seqnames(bpgr)) == as.character(seqnames(partner(bpgr)))
  insSeq = rep("", length(bpgr))
  if (!is.null(bpgr$insSeq)) {
    insSeq = bpgr$insSeq
  }
  if (!is.null(bpgr$insertSequence)) {
    insSeq = bpgr$insertSequence
  }
  more_ins_than_length = str_length(insSeq) >= abs(start(bpgr)-start(partner(bpgr))) * 0.7
  same_strand = strand(bpgr) == strand(partner(bpgr))
  result[isbp] = ifelse(!same_chr, "ITX", # inter-chromosomosal
          ifelse(more_ins_than_length, "INS", # TODO: improve classification of complex events
            ifelse(same_strand, "INV",
              ifelse(xor(start(bpgr) < start(partner(bpgr)), strand(bpgr) == "-"), "DEL", "DUP"))))
  return(result)
}

load_full_gridss_gr = function(filename, directory="", filter=c("somatic", "qual")) {
	vcf = readVcf(paste0(directory, filename), "hg19")
	if ("somatic" %in% filter) {
		vcf = gridss_somatic_filter(vcf)
	}
	if ("qual" %in% filter) {
		vcf = gridss_qual_filter(vcf)
	}
	gr = breakpointRanges(vcf)
	gr$filename = filename
	gr$af = gridss_somatic_af(gr, vcf)
	info_df = info(vcf[gr$sourceId,])
	info_df$filtered = should_filter(vcf[gr$sourceId,])
	mcols(gr) = bind_cols(as.data.frame(mcols(gr)), as.data.frame(info_df))
	return(gr)
}
full_gridss_annotate_gr = function(gr, vcf, geno_suffix=c(".normal", ".tumour")) {
	vcf = vcf[gr$sourceId,]
	g = geno(vcf)
	extra_df <- as.data.frame(as.data.frame(info(vcf)))
	for (i in seq_along(geno_suffix)) {
		if (!is.na(geno_suffix[i])) {
			gdf = data.frame(QUAL=g$QUAL[,i])
			for (col in names(g)) {
				gdf[[paste0(col, geno_suffix[i])]] = g[[col]][,i]
			}
			extra_df <- bind_cols(extra_df, gdf)
		}
	}
	mcols(gr) = bind_cols(as.data.frame(mcols(gr)), extra_df)
}

transitive_of = function(gr, max_traversal_length, max_positional_error, ...) {
	paths = transitive_paths(gr, max_traversal_length, max_positional_error, ...)
	result = rep(NA_character_, length(gr))
	result[paths$ordinal] = paths$bp_path
	return(result)
}
#' Calculates all transitive paths for the given set of breakpoints
#' @param gr breakpoint GRanges
#' @param max_traversal_length max length of sequence to traverse
#' @param max_positional_error maximum position error in breakpoint calls
#' @param max_depth maximum number of breakpoints to traverse.
#' Highly connected breakpoint sets are O(e^max_depth)
#' @param group_by_col_name column name to separate GRanges by.
#' The default of sampleId causes each sample to calculated independently
transitive_paths = function(gr, max_traversal_length, max_positional_error, max_depth=5, group_by_col_name="sampleId") {
	gr$ordinal = seq(length(gr))
	hits = as.data.frame(findOverlaps(gr, gr, maxgap=max_traversal_length, ignore.strand=TRUE))
	hits = hits %>% filter(hits$subjectHits != hits$queryHits)
	if (group_by_col_name %in% names(mcols(gr))) {
		hits = hits %>% filter(mcols(gr)[[group_by_col_name]][queryHits] == mcols(gr)[[group_by_col_name]][subjectHits])
	}
	hits = hits %>% mutate(
			query_strand = as.character(strand(gr)[queryHits]),
			subject_strand = as.character(strand(gr)[subjectHits])) %>%
		mutate(
			upstream_distance=(start(gr)[queryHits] - start(gr)[subjectHits]) * ifelse(query_strand =="+", 1, -1))

	terminal_hits = hits %>%
		filter(query_strand == subject_strand) %>%
		filter(abs(upstream_distance) < max_positional_error) %>%
		mutate(current=queryHits, terminal_ordinal=subjectHits, terminal_distance=upstream_distance) %>%
		dplyr::select(current, terminal_ordinal, terminal_distance)
	transitive_hits = hits %>%
		filter(query_strand != subject_strand) %>%
		filter(upstream_distance > -max_positional_error) %>%
		mutate(
			transitive_start=queryHits,
			transitive_end=partner(gr)$ordinal[subjectHits],
			transitive_distance=upstream_distance,
			path_name=names(gr)[subjectHits]) %>%
		dplyr::select(transitive_start, transitive_end, transitive_distance, path_name)

	activedf = data.frame(
			ordinal=gr$ordinal,
			terminal_ordinal=partner(gr)$ordinal) %>%
		dplyr::inner_join(terminal_hits, by=c("ordinal"="current"), suffix=c("", ".first")) %>%
		mutate(
			distance=-terminal_distance,
			bp_path=names(gr)[terminal_ordinal.first],
			current=partner(gr)$ordinal[terminal_ordinal.first],
			path_length=1) %>%
		dplyr::select(ordinal, current, bp_path, path_length, terminal_ordinal, distance)

	resultdf = data.frame(ordinal=integer(), name=character(), bp_path=character(), path_length=integer(), distance=integer(), stringsAsFactors=FALSE)
	while(nrow(activedf) > 0 & max_depth > 0) {
		activedf = activedf %>% dplyr::inner_join(transitive_hits, by=c("current"="transitive_start")) %>%
			mutate(
				distance = distance + transitive_distance,
				current=transitive_end,
				bp_path=paste(bp_path, path_name),
				path_length=path_length + 1) %>%
			dplyr::select(ordinal, current, bp_path, path_length, terminal_ordinal, distance) %>%
			filter(distance <= max_traversal_length + max_positional_error & distance >= -max_positional_error) %>%
			filter(current != ordinal & current != terminal_ordinal) # don't follow loops
		current_terminal = activedf %>%
			dplyr::inner_join(terminal_hits, by=c("current"="current", "terminal_ordinal"="terminal_ordinal")) %>%
			mutate(distance = distance + terminal_distance,
						 name=names(gr)[ordinal]) %>%
			dplyr::select(ordinal, name, bp_path, path_length, distance)
		resultdf = bind_rows(resultdf, current_terminal)
		max_depth = max_depth - 1
	}
	return(resultdf)
}
test_that("transitive_paths", {
	#debugonce(transitive_paths)
	gr = GRanges(seqnames=c("A", "C", "A", "B", "B", "C", "B", "C"),
							 ranges=IRanges(start=c(1, 5000, 1, 1000, 2000, 5100, 2000, 5100), width=1),
							 strand=c("+", "-", "+", "-", "+", "-", "+", "-"))
	names(gr) = c("AC1", "AC2", "AB1", "AB2", "BC1", "BC2", "diffSampleBC1", "diffSampleBC2")
	gr$partner = c("AC2", "AC1", "AB2", "AB1", "BC2", "BC1", "diffSampleBC2", "diffSampleBC1")
	gr$sampleId = c(1, 1, 1, 1, 1, 1, 2, 2)
	paths = transitive_paths(gr, 2000, 200)
	expect_equal(2, nrow(paths))
	expect_equal(rep(1000-100, 2), paths$distance)
	expect_equal(c("AC1", "AC2"), paths$name)
})
get_db_comparision_df = function(dbExisting, dbNew, suffix=c(".old", ".new"), sampleIds=NULL, line_annotation_bed=NULL) {
  if (is.null(sampleIds)) {
    common_sample_ids = query_structural_variants_samples(dbExisting)
    common_sample_ids <- common_sample_ids[common_sample_ids %in% query_structural_variants_samples(dbNew)]
  } else {
    common_sample_ids = sampleIds
  }
  grex <- query_structural_variants_for_sample_as_granges(dbExisting, common_sample_ids)
  grnew <- query_structural_variants_for_sample_as_granges(dbNew, common_sample_ids)
  grex$transitive=transitive_of(grex, 2000, 300)
  grnew$transitive=transitive_of(grnew, 2000, 300)
  # make a proxy QUAL since gr_join_to_df needs it to resolve matches in favour of the 'better' one
  grex$QUAL <- ifelse(is.na(grex$ploidy), grex$af, grex$ploidy)
  grnew$QUAL <- ifelse(is.na(grnew$ploidy), grnew$af, grnew$ploidy)
  grex$type = simpleEventType(grex)
  grnew$type = simpleEventType(grnew)

  if (!is.null(line_annotation_bed)) {
    grex$isline <- overlapsAny(grex, line_annotation_bed, maxgap=10000)
    grnew$isline <- overlapsAny(grnew, line_annotation_bed, maxgap=10000)
  }
  fullmatchdf <- gr_join_to_df(grex, NULL, grnew, NULL, maxgap=500, sizemargin=2, suffix=suffix)
  fullmatchdf = fullmatchdf %>%
    mutate(called=c("Neither", "Existing", "New", "Both")
           [1+ifelse(!is.na(fullmatchdf[[paste0("start", suffix[1])]]), 1, 0) + ifelse(!is.na(fullmatchdf[[paste0("start", suffix[2])]]), 2, 0)]) %>%
    mutate(transitive=c("Neither", "Existing", "New", "Both")
           [1+ifelse(!is.na(fullmatchdf[[paste0("transitive", suffix[1])]]), 1, 0) + ifelse(!is.na(fullmatchdf[[paste0("transitive", suffix[2])]]), 2, 0)]) %>%
    mutate(type=ifelse(!is.na(fullmatchdf[[paste0("type", suffix[1])]]), fullmatchdf[[paste0("type", suffix[1])]], fullmatchdf[[paste0("type", suffix[2])]]),
           sampleId=ifelse(!is.na(fullmatchdf[[paste0("sampleId", suffix[1])]]), fullmatchdf[[paste0("sampleId", suffix[1])]], fullmatchdf[[paste0("sampleId", suffix[2])]]),
           orientation=ifelse(!is.na(fullmatchdf[[paste0("orientation", suffix[1])]]), fullmatchdf[[paste0("orientation", suffix[1])]], fullmatchdf[[paste0("orientation", suffix[2])]])) %>%
    mutate(svlen=ifelse(!is.na(fullmatchdf[[paste0("svlen", suffix[2])]]), fullmatchdf[[paste0("svlen", suffix[2])]], fullmatchdf[[paste0("svlen", suffix[1])]]))
  if (!is.null(line_annotation_bed)) {
    line1 <- fullmatchdf[[paste0("isline", suffix[1])]]
    line2 <- fullmatchdf[[paste0("isline", suffix[2])]]
    fullmatchdf = fullmatchdf %>% mutate(isline=(!is.na(line1) & line1) | (!is.na(line2) & line2))
  }
  return(fullmatchdf)
}

#' @description Determines which A-C transitives calls can be explained by an A-B-C.
#' Transitive calls are caused by fragments that completely span the B fragment.
#' For paired end sequencing, transitive calls will be IMPRECISE as all split reads
#' support the A-B and B-C breakpoints.
#' @param gr input breakpoint GRanges object
#' @param max_traversed_length maximum length of sequence in the traversed breakpoint chain.
#' @param min_segment_length minimum length of traversed sequence. Aligners typically cannot align less than around 20 bases.
#' @param allow_loops allow breakpoints to be traversed multiple times
#' @param max_hops maximum number of breakends to traverse
#' @param report When report is "shortest", only the shortest transitive breakpoint chain is reported,
#' "all" reports all chains,
#' and "max2" reports the two chains with the least number of hops
#' @return id of the
transitive_breakpoints <- function(
    gr,
    find_transitive_for=rep(TRUE, length(gr)),
    can_traverse_through=rep(TRUE, length(gr)),
    max_traversed_length=1000,
    min_segment_length=0,
    transitive_call_slop=100,
    allow_loops=FALSE,
    max_hops=4,
    max_intermediate_paths=1000,
    max_active_paths=100000,
    report=c("shortest", "max2", "all")) {
  ordinal_lookup = seq_len(length(gr))
  names(ordinal_lookup) = names(gr)
  partner_lookup = ordinal_lookup[gr$partner]
  next_df = .traversal_next_breakpoint(gr, max_traversed_length, min_segment_length) %>%
    # prevent self-intersections
    filter(queryHits != subjectHits & queryHits != partner_lookup[subjectHits]) %>%
    # incorporate breakpoint inserted sequence into the traversal length
    mutate(
      min_traversed=min_traversed + gr$insLen[subjectHits],
      max_traversed=max_traversed + gr$insLen[subjectHits],
      source_from=partner_lookup[queryHits],
      source_to=queryHits,
      dest_from=subjectHits,
      dest_to=partner_lookup[subjectHits]) %>%
    dplyr::select(-queryHits, -subjectHits) %>%
    filter(
      can_traverse_through[source_from] &
      can_traverse_through[source_to] &
      can_traverse_through[dest_from] &
      can_traverse_through[dest_to])

  terminal_df = .adjacent_breakends(gr, gr, maxgap=transitive_call_slop, allowed_orientation=c("--", "++")) %>%
    filter(
      queryHits != subjectHits & queryHits != partner_lookup[subjectHits] &
      min_traversed < transitive_call_slop & max_traversed > -transitive_call_slop &
      find_transitive_for[queryHits] & can_traverse_through[subjectHits])

  result_df = data.frame(transitive_start=character(), transitive_end=character(), bp_path=character(), min_length=integer(), max_length=integer())
  # start with the first traversal away from the putative transitive call
  active_df = terminal_df %>%
    mutate(
      terminal_start=queryHits,
      terminal_end=partner_lookup[queryHits],
      current_to=partner_lookup[subjectHits],
      bp_path=names(gr)[subjectHits],
      min_length=min_traversed + gr$insLen[subjectHits],
      max_length=max_traversed + gr$insLen[subjectHits]) %>%
    dplyr::select(terminal_start, terminal_end, current_to, bp_path, min_length, max_length) %>%
    filter(find_transitive_for[terminal_start])
  path_key=integer()
  path_value1=integer()
  path_value2=integer()
  i = 0
  while (nrow(active_df) > 0 & i < max_hops) {
    # continue traversing
    active_df = active_df %>%
      dplyr::select(terminal_start, terminal_end, current_to, bp_path, min_length, max_length) %>%
      inner_join(next_df, by=c("current_to"="source_to")) %>%
      filter(allow_loops | !str_detect(bp_path, stringr::fixed(names(gr)[dest_from]))) %>%
      mutate(
        bp_path=paste0(bp_path, ";", names(gr)[dest_from]),
        current_to=dest_to,
        min_length=min_length + min_traversed,
        max_length=max_length + max_traversed) %>%
      dplyr::select(terminal_start, terminal_end, current_to, bp_path, min_length, max_length) %>%
      filter(min_length < max_traversed_length)
    if (FALSE && report %in% c("shortest", "max2")) {
      # Only follow the two shortest paths to any given intermediate node
      # technically we're going to follow all paths of the shortest two lengths
      # but as this is just an optimisation to prevent exponential expansion of
      # highly connected path graphs it's good enough
      # NB: only works on up to 2**16 nodes
      key = active_df$current_to * length(gr) + active_df$terminal_start
      key_offset = match(key, path_key)
      value = active_df$min_length
      is_only_path = is.na(key_offset)
      is_best_path = !is_only_path & value < path_value1[key_offset]
      is_second_best_path = !is_only_path & !is_best_path & value < path_value2[key_offset]
      # Update second value
      path_value2[key_offset[is_second_best_path]] = value[is_second_best_path]
      # Update first value
      path_value2[key_offset[is_best_path]] = path_value1[key_offset[is_best_path]]
      path_value1[key_offset[is_best_path]] = value[is_best_path]
      # Add new best values
      path_key = c(path_key, key[is_only_path])
      path_value1 = c(path_value1, value[is_only_path])
      path_value2 = c(path_value2, rep(2 * max_traversed_length, sum(is_only_path))) # placeholder
      active_df = active_df %>% filter(is_only_path | is_best_path | is_second_best_path)
    }
    active_df = active_df %>%
      group_by(terminal_start) %>%
      top_n(as.integer(min(max_intermediate_paths, 1 + max_active_paths / (length(unique(active_df$terminal_start)) + 1))), wt=min_length) %>%
      ungroup()
    # check for terminal completion
    active_df = active_df %>% left_join(terminal_df, by=c("current_to"="queryHits", "terminal_end"="subjectHits"))
    result_df = active_df %>%
      filter(!is.na(min_traversed)) %>%
      mutate(min_length=min_length + min_traversed,
             max_length=max_length + max_traversed,
             transitive_start=names(gr)[terminal_start],
             transitive_end=names(gr)[terminal_end]) %>%
      dplyr::select(transitive_start, transitive_end, bp_path, min_length, max_length) %>%
      bind_rows(result_df)
    if (report == "shortest") {
      # any further solutions will be longer so we don't need to consider them
      active_df = active_df %>% filter(is.na(min_traversed))
      best_min = result_df
    }
    if (report == "max2") {
      found2 = result_df %>%
        group_by(transitive_start) %>%
        summarise(n=n()) %>%
        filter(n >= 2) %>%
        pull(transitive_start)
      active_df = active_df %>% filter(!(names(gr)[terminal_start] %in% found2))
    }
    i = i + 1
  }
  if (report == "shortest") {
    # Just the shortest result for each transitive call
    result_df = result_df %>%
      group_by(transitive_start) %>%
      top_n(1, min_length) %>%
      ungroup()
  }
  if (report == "max2") {
    result_df = result_df %>%
      group_by(transitive_start) %>%
      top_n(1, min_length) %>%
      ungroup()
  }
  return(result_df)
}
.traversal_next_breakpoint <- function(gr, max_traversed_length, min_segment_length) {
  .adjacent_breakends(gr, gr, max_traversed_length, allowed_orientation=c("-+", "+-")) %>%
    filter(max_traversed > min_segment_length) %>%
    mutate(min_traversed=pmax(min_segment_length, min_traversed))
}
#' @description determines which breakends are near the given breakend
#' @param maxgap maximum distance between adjacent breakends
#' @param allowed_orientation allowed breakend orientation of the query then subject breakends
#' @return adjacent breakends of the given orientations.
#' Traversal distances are defined as the number of number of bases traversed into the DNA segment
#' arriving from the query breakend until the subject breakend position is reached.
#' Note that the subject breakend orientation does not affect the traversal distance.
.adjacent_breakends <- function(query, subject, maxgap, allowed_orientation=c("--", "-+", "+-", "++")) {
  allowed_orientation = match.arg(allowed_orientation, several.ok=TRUE)
  findOverlaps(query, subject, maxgap=maxgap, ignore.strand=TRUE) %>%
    as.data.frame() %>%
    dplyr::filter(paste0(strand(query[queryHits]), strand(subject[subjectHits])) %in% allowed_orientation) %>%
    dplyr::mutate(
      start_end_traversed=(end(subject[subjectHits]) - start(query[queryHits])) * ifelse(as.logical(strand(query[queryHits])=="-"), 1, -1),
      end_start_traversed=(start(subject[subjectHits]) - end(query[queryHits])) * ifelse(as.logical(strand(query[queryHits])=="-"), 1, -1)) %>%
    dplyr::mutate(min_traversed=pmin(start_end_traversed, end_start_traversed)) %>%
    dplyr::mutate(max_traversed=pmax(start_end_traversed, end_start_traversed)) %>%
    dplyr::select(-start_end_traversed, -end_start_traversed)
}
#' Adjusts the nominal breakpoint position towards
#' @param  align alignment position within any interval of uncertainty.
#' In the case of an even width interval, centre alignment adjusts to the centre
#' position closest to the initial location.
#' Adjusts the nominal position of a breakpoints
align_breakpoints <- function(vcf, align=c("centre"), is_higher_breakend=names(vcf) < info(vcf)$MATEID) {
  if (length(vcf) == 0) {
    return(vcf)
  }
  align = match.arg(align)
  if (!all(elementNROWS(info(vcf)$CIPOS) == 2)) {
    stop("CIPOS not specified for all variants.")
  }
  is_higher_breakend[is.na(is_higher_breakend)] = FALSE
  nominal_start = start(rowRanges(vcf))
  cipos = t(matrix(unlist(info(vcf)$CIPOS), nrow=2))
  ciwdith = cipos[,2] - cipos[,1]
  orientations = .vcfAltToStrandPair(rowRanges(vcf)$ALT)
  if (align == "centre") {
    citargetpos = nominal_start + cipos[,1] + ciwdith / 2.0
    adjust_by = citargetpos - nominal_start
    adjust_in_opposite_direction_to_partner = orientations %in% c("--", "++")
    adjust_by = ifelse(is_higher_breakend & adjust_in_opposite_direction_to_partner, ceiling(adjust_by), floor(adjust_by))
  } else {
    stop("Only centre alignment is currently implemented.")
  }
  isbp = str_detect(VariantAnnotation::fixed(vcf)$ALT, "[\\]\\[]")
  is_adjusted_bp =  isbp & !is.na(adjust_by) & adjust_by != 0
  rowRanges(vcf) = shift(rowRanges(vcf), ifelse(!is_adjusted_bp, 0, adjust_by))
  info(vcf)$CIPOS = info(vcf)$CIPOS - adjust_by
  if (!is.null(info(vcf)$CIEND)) {
    info(vcf)$CIEND = info(vcf)$CIEND - adjust_by
  }
  if (!is.null(info(vcf)$IHOMPOS)) {
    info(vcf)$IHOMPOS = info(vcf)$IHOMPOS - adjust_by
  }
  alt = unlist(rowRanges(vcf)$ALT)
  partner_alt = stringr::str_match(alt, "^([^\\]\\[]*)[\\]\\[]([^:]+):([0-9]+)([\\]\\[])([^\\]\\[]*)$")
  # [,2] anchoring bases
  # [,3] partner chr
  # [,4] old partner position
  partner_pos = ifelse(is.na(partner_alt[,4]), NA_integer_, as.integer(partner_alt[,4])) + ifelse(adjust_in_opposite_direction_to_partner, -adjust_by, adjust_by)
  # [,5] partner orientation
  # [,6] anchoring bases
  # adjust ALT for breakpoints. anchoring bases get replaced with N since we don't know
  VariantAnnotation::fixed(vcf)$ALT = as(ifelse(!is_adjusted_bp, alt,
                                                paste0(
                                                  str_pad("", stringr::str_length(partner_alt[,2]), pad="N"),
                                                  partner_alt[,5],
                                                  partner_alt[,3],
                                                  ":",
                                                  format(partner_pos, trim=TRUE, scientific=FALSE),
                                                  partner_alt[,5],
                                                  str_pad("", stringr::str_length(partner_alt[,6]), pad="N"))), "CharacterList")
  info(vcf)$CIRPOS = NULL # TODO: remove CIRPOS from GRIDSS entirely
  return(vcf)
}
.vcfAltToStrandPair = function(alt) {
  chralt = unlist(alt)
  ifelse(startsWith(chralt, "."), "-",
         ifelse(endsWith(chralt, "."), "+",
                ifelse(startsWith(chralt, "]"), "-+",
                       ifelse(startsWith(chralt, "["), "--",
                              ifelse(endsWith(chralt, "]"), "++",
                                     ifelse(endsWith(chralt, "["), "+-", ""))))))

}

readVcf = function(file, ...) {
  raw_vcf = VariantAnnotation::readVcf(file=file, ...)
  #id = read_tsv(file, comment="#", col_names=c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", seq_len(ncol(geno(raw_vcf)[[1]]))), cols_only(ID=col_character()))$ID
  #assertthat::assert_that(all(alt(raw_vcf) != ""), "VariantAnnotation 1.29.11 or later is required")
  if (!all(unlist(alt(raw_vcf)) != "")) {
    write("Performing work-around for https://github.com/Bioconductor/VariantAnnotation/issues/8", stderr())
    alt = read_tsv(file, comment="#", col_names=c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", seq_len(ncol(geno(raw_vcf)[[1]]))), cols_only(ALT=col_character()))$ALT
    VariantAnnotation::fixed(raw_vcf)$ALT = CharacterList(lapply(as.character(alt), function(x) x))
  }
  # Work-around for https://github.com/PapenfussLab/gridss/issues/156
  # since we don't have all the info, we'll just pro-rata
  # is.nanan = function(x) is.na(x) | is.nan(x)
  # bp_pro_rata = (geno(raw_vcf)$ASSR + geno(raw_vcf)$ASRP) / rowSums(geno(raw_vcf)$ASSR + geno(raw_vcf)$ASRP)
  # be_pro_rata = (geno(raw_vcf)$BASSR + geno(raw_vcf)$BASRP) / rowSums(geno(raw_vcf)$BASSR + geno(raw_vcf)$BASRP)
  # bp_pro_rata[is.nanan(bp_pro_rata)] = 0
  # be_pro_rata[is.nanan(be_pro_rata)] = 0
  # geno(raw_vcf)$ASQ[is.nanan(geno(raw_vcf)$ASQ)] = (info(raw_vcf)$ASQ * bp_pro_rata)[is.nanan(geno(raw_vcf)$ASQ)]
  # geno(raw_vcf)$RASQ[is.nanan(geno(raw_vcf)$RASQ)] = (info(raw_vcf)$RASQ * bp_pro_rata)[is.nanan(geno(raw_vcf)$RASQ)]
  # geno(raw_vcf)$CASQ[is.nanan(geno(raw_vcf)$CASQ)] = (info(raw_vcf)$CASQ * bp_pro_rata)[is.nanan(geno(raw_vcf)$CASQ)]
  # geno(raw_vcf)$BAQ[is.nanan(geno(raw_vcf)$BAQ)] = (info(raw_vcf)$BAQ * be_pro_rata)[is.nanan(geno(raw_vcf)$BAQ)]
  # geno(raw_vcf)$QUAL[is.nanan(geno(raw_vcf)$QUAL)] = (
  #   geno(raw_vcf)$ASQ +
  #   geno(raw_vcf)$RASQ +
  #   geno(raw_vcf)$CASQ +
  #   geno(raw_vcf)$RPQ +
  #   geno(raw_vcf)$SRQ +
  #   geno(raw_vcf)$IQ
  # )[is.nanan(geno(raw_vcf)$QUAL)]
  # geno(raw_vcf)$BQ[is.nanan(geno(raw_vcf)$BQ)] = (
  #   geno(raw_vcf)$BAQ +
  #   geno(raw_vcf)$BUMQ +
  #   geno(raw_vcf)$BSCQ
  # )[is.nanan(geno(raw_vcf)$BQ)]
  raw_vcf = fix_parid(raw_vcf)
  raw_vcf = flatten_mateid(raw_vcf)
  return(raw_vcf)
}
cached_read_file = function(file, read_function) {
  cache_filename = paste0(file, ".cached.parsed.rds")
  if (file.exists(cache_filename)) {
    write(paste(Sys.time(), "Loading from cache ", cache_filename), stderr())
    result = readRDS(cache_filename)
  } else {
    result = read_function(file)
    saveRDS(result, file=cache_filename)
  }
}

read_gridss_breakpoint_pon = function(file) {
  df = read_tsv(file,
                col_names=c("chr1", "start1", "end1", "chr2", "start2", "end2", "name", "score", "strand1", "strand2"),
                col_types="ciiciicccc")
  gro = GRanges(
    seqnames=df$chr1,
    ranges=IRanges(
      start=df$start1 + 1,
      end=df$end1),
    strand=df$strand1,
    partner=paste0(seq_len(nrow(df)), "h"),
    score=df$score)
  names(gro) = paste0(seq_len(nrow(df)), "o")
  grh = GRanges(
    seqnames=df$chr2,
    ranges=IRanges(
      start=df$start2 + 1,
      end=df$end2),
    strand=df$strand2,
    partner=paste0(seq_len(nrow(df)), "o"),
    score=df$score)
  names(grh) = paste0(seq_len(nrow(df)), "h")
  return(c(gro, grh))
}
#' takes two lists of lists of the same shape
#' zips each corresponding element using ZIP_FUN
#' Then for each original list, applies AGGREGATE_FUN to the result of the zip
#' ZIP_FUN must be a vectorised function and is applied once
#' AGGREGATE_FUN is applied once for each non-empty nested list
#' @return a vector of the same length as the input lists.
#' If the input list has no elements,
#' the resultant vector returns NA without applying AGGREGATE_FUN to the empty list
zip_aggregate_pair_of_list_of_list = function(list1, list2, ZIP_FUN, AGGREGATE_FUN) {
  if (any(elementNROWS(list1) != elementNROWS(list2))) {
    stop("Lists are not of the same shape")
  }
  vec1 = unlist(list1)
  vec2 = unlist(list2)
  zipped = ZIP_FUN(vec1, vec2)
  result = rep(NA, length(list1))
  aggregatedf = data.frame(
      ordinal = rep(seq_len(length(list1)), times=elementNROWS(list1)),
      value = zipped) %>%
    group_by(ordinal) %>%
    summarise(value = AGGREGATE_FUN(value))
  result[aggregatedf$ordinal] = aggregatedf$value
  return(result)
}
linked_assemblies <- function(vcf, exclude_breakends_without_local_assembly=TRUE) {
  asm_linked_df <- data.frame(
    sourceId=names(vcf),
    mateid=info(vcf)$MATEID,
    start=start(rowRanges(vcf)),
    asm=zip_aggregate_pair_of_list_of_list(
      info(vcf)$BEID,
      info(vcf)$BEIDL,
      ZIP_FUN = function(x, y) paste(x, y, sep="/"),
      AGGREGATE_FUN = function(x) paste(x, collapse=";")),
    stringsAsFactors = FALSE
  ) %>%
    separate_rows(asm, sep=";") %>%
    filter(asm != "")

  asm_paired_df = .linked_assemblies_fix_indels(asm_linked_df)
  if (exclude_breakends_without_local_assembly) {
    asm_paired_df = .linked_assemblies_fix_breakends(vcf, asm_linked_df)
  }
  asm_paired_df = asm_paired_df %>%
    dplyr::select(sourceId, linked_by=asm) %>%
    group_by(linked_by) %>%
    filter(n() == 2) %>%
    ungroup()
  asm_paired_df = exclude_self_links(asm_paired_df, vcf)
  return(asm_paired_df)
}
exclude_self_links = function(linked_df, vcf) {
  self_links = linked_df %>% dplyr::select(sourceId, linked_by) %>%
    mutate(mateid = info(vcf[sourceId])$MATEID) %>%
    filter(!is.na(mateid)) %>%
    group_by(linked_by) %>%
    filter(length(unique(c(sourceId, mateid))) != 4) %>%
    pull(linked_by)
  return(linked_df %>% filter(!(linked_by %in% self_links)))

}
#' assemblies spanning SV indels will have the same asm
#' linking id as the flanking variants as well as having
#' the same id on both sides of the indel
#' we need to split these up so asm identifiers are unique
.linked_assemblies_fix_indels = function(asm_linked_df) {
  asm_linked_df %>%
    group_by(asm) %>%
    arrange(start) %>%
    mutate(
      is_indel_start = (lead(mateid) == sourceId) %na% FALSE,
      is_indel_end = (lag(mateid) == sourceId) %na% FALSE,
      requires_asm_renaming=n() > 2,
      paired_with_next=row_number() == 1 | is_indel_end) %>%
    # trim starting indel bound
    filter(!(is_indel_start & is.na(lag(sourceId)))) %>%
    # trim ending indel bound
    filter(!(is_indel_end & is.na(lead(sourceId)))) %>%
    mutate(asm_new=
      ifelse(!requires_asm_renaming, asm,
      paste(asm, (row_number() + ifelse(paired_with_next, 1, 0)) / 2, sep="-"))) %>%
    ungroup() %>%
    dplyr::select(sourceId, asm=asm_new)
}
#' Linked breakends require an additional assembly to ensure
#' we actually have a local assembly at the breakend
.linked_assemblies_fix_breakends = function(vcf, asm_linked_df) {
  i = info(vcf)
  breakends_with_less_than_two_assemblies = names(vcf)[
    is.na(i$MATEID) & i$AS + i$RAS + i$CAS < 2]
  asm_linked_df %>% filter(!(sourceId %in% breakends_with_less_than_two_assemblies))
}
require(testthat)
test_that(".linked_assemblies_fix_indels corrects indels", {
  asm_linked_df = data.frame(
    sourceId=c("in", "indel1o", "indel1h","indel2o", "indel2h", "out", "indel3o", "indel3h", "other_left", "other_right"),
    mateid=c("a", "indel1h", "indel1o","indel2h", "indel2o", NA, "indel3h", "indel3o", "b", "c"),
    start=c(1, 2, 3, 4, 5, 6, 2, 4, 3, 4),
    asm=c(rep("asm1", 6), rep("asm2", 2), rep("asm3", 2)),
    stringsAsFactors=FALSE)
  result = .linked_assemblies_fix_indels(asm_linked_df)
  # should remove isolated indel
  expect_that(nrow(result), equals(8))
  # should now have 4 links due to:
  # - removing asm2
  # - splitting asm1
  expect_that(length(unique(result$asm)), equals(4))
})


transitive_calls = function(vcf, bpgr,
  report=c("shortest", "max2", "all")
  ) {
  is_imprecise = info(vcf[names(bpgr)])$IMPRECISE
  is_imprecise[is.na(is_imprecise)] = FALSE
  # transitive calling reduction
  transitive_df = transitive_breakpoints(
    # only find transitive paths for the imprecise calls
    find_transitive_for=is_imprecise,
    # only traverse through precise paths
    # TODO: should we relax this criterion?
    can_traverse_through=!is_imprecise,
    bpgr,
    min_segment_length=20,
    report=report)
  transitive_df = transitive_df %>%
    #filter(info(vcf[transitive_df$transitive])$IMPRECISE) %>%
    mutate(full_path=bp_path)
  if (nrow(transitive_df) != 0) {
    transitive_df = transitive_df %>%
      filter(str_detect(transitive_start, "o$")) %>%
      group_by(transitive_start) %>%
      mutate(has_multiple_paths=n() > 1) %>%
      mutate(pathid=row_number()) %>%
      separate_rows(bp_path, sep=";") %>%
      group_by(transitive_start, full_path) %>%
      mutate(ordinal=row_number()) %>%
      mutate(
        remote_sourceId=bpgr[bp_path]$partner,
        prev_sourceId=ifelse(is.na(lag(remote_sourceId)), transitive_start, lag(remote_sourceId)),
        next_sourceId=ifelse(is.na(lead(bp_path)), transitive_end, lead(bp_path))) %>%
      ungroup()
    transitive_df = bind_rows(
      transitive_df %>%
        mutate(linked_by=paste(transitive_start, ordinal, sep="/"), sourceId=bp_path) %>%
        dplyr::select(sourceId, linked_by, has_multiple_paths, pathid, transitive_start, transitive_end),
      transitive_df %>%
        mutate(linked_by=paste(transitive_start, ordinal, sep="/"), sourceId=prev_sourceId) %>%
        dplyr::select(sourceId, linked_by, has_multiple_paths, pathid, transitive_start, transitive_end),
      transitive_df %>%
        mutate(linked_by=paste(transitive_start, ordinal+1, sep="/"), sourceId=remote_sourceId) %>%
        dplyr::select(sourceId, linked_by, has_multiple_paths, pathid, transitive_start, transitive_end),
      transitive_df %>%
        mutate(linked_by=paste(transitive_start, ordinal+1, sep="/"), sourceId=next_sourceId) %>%
        dplyr::select(sourceId, linked_by, has_multiple_paths, pathid, transitive_start, transitive_end)) %>%
      distinct()
  } else {
    # work-around for https://github.com/tidyverse/tidyr/issues/470
    transitive_df = data.frame(linked_by=character(), sourceId=character(), has_multiple_paths=logical(), pathid=numeric(), stringsAsFactors=FALSE)
  }
}
#' @description Calculates a somatic score for the given variant.
#' The score is defined as the log likelihood ratio of the observed
#' variant read counts being caused by tumour contamination of the normal,
#' and the variant being heterozygous in the normal.
#' @param vcf GRIDSS VCF object to calculate somatic scores for
#' @param normalOrdinal Name or index of the normal genotype in the VCF
#' @param tumourOrdinal Name or index of the tumour genotype in the VCF
#' @param contamination_rate estimated rate of tumour contamination. Typically,
#' this should be set to the tumour variant allele frequency multiplied by the
#' copy number ratio of the tumour and normal.
#' @param normal_ploidy Expected normal ploidy.
gridss_breakpoint_somatic_llr = function(vcf, normalOrdinal, tumourOrdinal, contamination_rate, normal_ploidy = 2, bpgr=breakpointRanges(vcf)) {
  g = geno(vcf)
  df = data.frame(
    is_bp = names(vcf) %in% names(bpgr),
    is_short_bp = names(vcf) %in% names(bpgr)[is_short_deldup(bpgr)])
  row.names(df) = names(vcf)
  partner_lookup = match(as.character(bpgr$partner), row.names(df))
  names(partner_lookup) = names(bpgr)
  df$partner = NA_integer_
  df[df$is_bp, "partner"] = partner_lookup[row.names(df)[df$is_bp]]
  df = df %>% mutate(
    n_ref = .genosum(g$REF,c(normalOrdinal)) + ifelse(is_short_bp, 0, .genosum(g$REFPAIR,c(normalOrdinal))),
    n_var = ifelse(is_bp, .genosum(g$VF,c(normalOrdinal)), .genosum(g$BVF,c(normalOrdinal))))
  df = df %>% mutate(
    n_ref = ifelse(!is_bp, n_ref, n_ref + df$n_ref[partner]),
    # assume variant ploidy is 1, then:
    germline_het_af = ifelse(is_bp, 1 / (1 + 2*(normal_ploidy - 1)), 1 / normal_ploidy),
    contamination_af = contamination_rate
  )
  df = df %>% mutate(
    germline_het_log_p = dbinom(n_var, n_ref + n_var, germline_het_af, log=TRUE),
    contamination_log_p = dbinom(n_var, n_ref + n_var, contamination_af, log=TRUE)
  )
  return (df$contamination_log_p - df$germline_het_log_p)
}

passes_final_filters = function(vcf, include.existing.filters=TRUE) {
  return(
    (!include.existing.filters | rowRanges(vcf)$FILTER %in% c(".", "PASS")) &
    ifelse(is.na(info(vcf)$MATEID),
         rowRanges(vcf)$QUAL >= gridss.min_qual * gridss.single_breakend_multiplier,
         rowRanges(vcf)$QUAL >= gridss.min_qual))
}

linked_by_breakend_breakend_insertion_classification = function(begr, maxgap=gridss.insertion.maxgap) {
  if (is.null(begr$sampleId)) {
    begr$sampleId = rep("placeholder", length(begr))
  }
  hits = findOverlaps(begr, begr, maxgap=maxgap, ignore.strand=TRUE) %>%
    as.data.frame() %>%
    filter(begr$sampleId[queryHits] == begr$sampleId[subjectHits]) %>% # intra-sample
    filter(as.logical(strand(begr)[queryHits] != strand(begr)[subjectHits])) %>% # opposite strand
    # matching pairs with the best qual
    mutate(qqual=begr$QUAL[queryHits], squal=begr$QUAL[subjectHits]) %>%
    group_by(queryHits) %>%
    top_n(1, squal) %>%
    group_by(subjectHits) %>%
    top_n(1, qqual) %>%
    ungroup() %>%
    filter(queryHits %in% subjectHits & subjectHits %in% queryHits) %>%
    filter(queryHits < subjectHits) %>% # remove symmetry
    mutate(linked_by=paste0("bebeins", row_number()))
  return( bind_rows(
    hits %>% mutate(sourceId=names(begr)[queryHits]) %>% dplyr::select(sourceId, linked_by),
    hits %>% mutate(sourceId=names(begr)[subjectHits]) %>% dplyr::select(sourceId, linked_by)
  ) %>% distinct())
}
linked_by_breakpoint_breakend_insertion_classification = function(bpgr, begr, maxgap=gridss.insertion.maxgap) {
  if (is.null(begr$sampleId)) {
    begr$sampleId = rep("placeholder", length(begr))
  }
  if (is.null(bpgr$sampleId)) {
    bpgr$sampleId = rep("placeholder", length(bpgr))
  }
  hits = findOverlaps(bpgr, begr, maxgap=maxgap, ignore.strand=TRUE) %>%
    as.data.frame() %>%
    filter(bpgr$sampleId[queryHits] == begr$sampleId[subjectHits]) %>% # intra-sample
    filter(as.logical(strand(bpgr)[queryHits] != strand(begr)[subjectHits])) %>% # opposite strand
    # matching pairs with the best qual
    mutate(qqual=bpgr$QUAL[queryHits], squal=begr$QUAL[subjectHits]) %>%
    group_by(queryHits) %>%
    top_n(1, squal) %>%
    group_by(subjectHits) %>%
    top_n(1, qqual) %>%
    ungroup() %>%
    mutate(linked_by=paste0("bpbeins", row_number()))
  return( bind_rows(
    hits %>% mutate(sourceId=names(bpgr)[queryHits]) %>% dplyr::select(sourceId, linked_by),
    hits %>% mutate(sourceId=names(begr)[subjectHits]) %>% dplyr::select(sourceId, linked_by)
  ) %>% distinct())
}
linked_by_simple_inversion_classification = function(bpgr, maxgap=gridss.inversion.maxgap) {
  bpgr = bpgr[strand(bpgr) == strand(partner(bpgr))]
  bpgr = bpgr[seqnames(bpgr) == seqnames(partner(bpgr))]
  if (is.null(bpgr$sampleId)) {
    bpgr$sampleId = rep("placeholder", length(bpgr))
  }
  hits = as.data.frame(findBreakpointOverlaps(bpgr, bpgr, maxgap=maxgap, ignore.strand=TRUE, sizemargin=NULL, restrictMarginToSizeMultiple=NULL)) %>%
    filter(bpgr$sampleId[queryHits] == bpgr$sampleId[subjectHits]) %>% # intra-sample
    filter(as.logical(strand(bpgr)[queryHits] != strand(bpgr)[subjectHits])) %>% # opposite strand
    # matching pairs with the best qual
    mutate(qqual=bpgr$QUAL[queryHits], squal=bpgr$QUAL[subjectHits]) %>%
    group_by(queryHits) %>%
    top_n(1, squal) %>%
    group_by(subjectHits) %>%
    top_n(1, qqual) %>%
    ungroup() %>%
    filter(queryHits < subjectHits) %>%
    filter(queryHits < match(bpgr[queryHits]$partner, names(bpgr))) %>%
    mutate(linked_by=paste0("inv", row_number()))
    return( bind_rows(
      hits %>% mutate(sourceId=names(bpgr)[queryHits]) %>% dplyr::select(sourceId, linked_by),
      hits %>% mutate(sourceId=bpgr[queryHits]$partner) %>% dplyr::select(sourceId, linked_by),
      hits %>% mutate(sourceId=names(bpgr)[subjectHits]) %>% dplyr::select(sourceId, linked_by),
      hits %>% mutate(sourceId=bpgr[subjectHits]$partner) %>% dplyr::select(sourceId, linked_by)
    ) %>% distinct())
}
linked_by_dsb = function(bpgr, maxgap=gridss.dsb.maxgap, ...) {
  linked_by_adjacency(bpgr, maxgap=maxgap, select="unique", link_label="dsb")
}

linked_by_different_foldback_inversion_paths = function(gr, max_inversion_length = 1000) {
  if (is.null(gr$sampleId)) {
    gr$sampleId = rep("placeholder", length(gr))
  }
  gr_foldback = gr[!is.na(gr$partner)]
  gr_foldback = gr_foldback[abs(start(gr_foldback) - start(partner(gr_foldback))) < max_inversion_length]
  gr_foldback = gr_foldback[seqnames(gr_foldback) == seqnames(partner(gr_foldback))]
  gr_foldback = gr_foldback[strand(gr_foldback) == strand(partner(gr_foldback))]
  gr_foldback = gr_foldback[findOverlaps(gr_foldback, gr_foldback, maxgap=max_inversion_length, ignore.strand=FALSE) %>%
    as.data.frame() %>%
    filter(
      queryHits != subjectHits &
      gr_foldback$sampleId[queryHits] == gr_foldback$sampleId[subjectHits] &
      names(gr_foldback)[queryHits] == gr_foldback$partner[subjectHits]) %>%
    pull(queryHits)]
  insSeq = .insSeq(gr)

  hitdf = findOverlaps(gr, gr_foldback) %>%
    as.data.frame() %>%
    filter(
      #names(gr)[queryHits] != names(gr_foldback)[subjectHits]  &
        gr$sampleId[queryHits] == gr_foldback$sampleId[subjectHits] &
        str_length(insSeq[queryHits]) > 0) %>%
    mutate(
      spanningIns=insSeq[queryHits],
      spanningAnchor="", # ignore the sequence past the inversion for now
      foldbackIns=insSeq[names(gr_foldback)[subjectHits]],
      foldbackLength=abs(start(gr_foldback[subjectHits]) - start(partner(gr_foldback)[subjectHits])),
      foldbackAnchor=get_partner_anchor_sequence(gr_foldback[subjectHits], foldbackLength),
      spanningSeq=ifelse(strand(gr[queryHits]) =="+", paste0(spanningIns, spanningAnchor), paste0(spanningAnchor, spanningIns)),
      foldbackSeq = ifelse(strand(gr_foldback[subjectHits]) =="+", paste0(foldbackIns, foldbackAnchor), paste0(foldbackAnchor, foldbackIns)),
      foldbackInvertedSeq = ifelse(strand(gr_foldback[subjectHits]) =="+",
        paste0(foldbackIns, as.character(reverseComplement(DNAStringSet(foldbackAnchor)))),
        paste0(as.character(reverseComplement(DNAStringSet(foldbackAnchor))), foldbackIns)),
      targetLength=pmin(str_length(spanningIns), str_length(foldbackIns) + foldbackLength)) %>%
    mutate(
      # normalise seq to be traversing from the left
      spanningSeq = ifelse(strand(gr[queryHits])=="+", spanningSeq, as.character(reverseComplement(DNAStringSet((spanningSeq))))),
      foldbackSeq = ifelse(strand(gr_foldback[subjectHits]) =="+", foldbackSeq, as.character(reverseComplement(DNAStringSet((foldbackSeq))))),
      foldbackInvertedSeq = ifelse(strand(gr_foldback[subjectHits]) =="+", foldbackInvertedSeq, as.character(reverseComplement(DNAStringSet((foldbackInvertedSeq)))))) %>%
    mutate(
      spanningSeq = str_sub(spanningSeq, end=targetLength),
      foldbackSeq = str_sub(foldbackSeq, end=targetLength),
      foldbackInvertedSeq = str_sub(foldbackInvertedSeq, end=targetLength)) %>%
    mutate(
      editDistance=stringdist(spanningSeq, foldbackSeq, method="lv"),
      invEditDistance=stringdist(spanningSeq, foldbackInvertedSeq, method="lv"))
  return(hitdf)
}
sequence_common_prefix = function(gr, anchor_bases, bsgenome, ...) {
  if (is.null(gr$sampleId)) {
    gr$sampleId = rep("placeholder", length(gr))
  }
  insSeq = .insSeq(gr)
  anchor_sequence = get_partner_anchor_sequence(gr, anchor_bases, bsgenome)
  hitdf = findOverlaps(gr, gr, ...) %>%
    as.data.frame() %>%
    filter(
      queryHits < subjectHits &
        names(gr)[queryHits] != (gr$partner[subjectHits] %na% "NA_placeholder") &
        gr$sampleId[queryHits] == gr$sampleId[subjectHits]) %>%
    mutate(
      queryIns=insSeq[queryHits],
      subjectIns=insSeq[subjectHits],
      queryAnchor=anchor_sequence[names(gr)[queryHits]],
      subjectAnchor=anchor_sequence[names(gr)[subjectHits]],
      querySeq = ifelse(strand(gr[queryHits]) =="+", paste0(queryIns, queryAnchor), paste0(queryAnchor, queryIns)),
      subjectSeq = ifelse(strand(gr[subjectHits]) =="+", paste0(subjectIns, subjectAnchor), paste0(subjectAnchor, subjectIns))) %>%
    mutate(
      # normalise seq to be breakpoint on the left
      querySeq = ifelse(strand(gr[queryHits])=="+", querySeq, as.character(reverseComplement(DNAStringSet((querySeq))))),
      subjectSeq = ifelse(strand(gr[subjectHits])=="+", subjectSeq, as.character(reverseComplement(DNAStringSet((subjectSeq)))))) %>%
    mutate(
      breakendLength=pmin(str_length(querySeq), str_length(subjectSeq)),
      querySeq = str_sub(querySeq, end=breakendLength),
      subjectSeq = str_sub(subjectSeq, end=breakendLength)) %>%
    filter(
      breakendLength > 0) %>%
    mutate(
      edit_distance=stringdist(querySeq, subjectSeq, method="lv"),
      per_base_edit_distance=edit_distance/ breakendLength)
  if (!is.null(gr$id)) {
    hitdf = hitdf %>% mutate(
      qid=gr$id[queryHits],
      qbeid=gr$beid[queryHits],
      sid=gr$id[subjectHits],
      sbeid=gr$beid[subjectHits])
  }
  if (!is.null(gr$sourceId)) {
    hitdf = hitdf %>% mutate(
      qsourceId=gr$sourceId[queryHits],
      ssourceId=gr$sourceId[subjectHits])
  }
  return(hitdf)
}
.insSeq = function(gr) {
  insSeq = rep("", length(gr))
  if (!is.null(gr$insSeq)) {
    insSeq = gr$insSeq
  }
  if (!is.null(gr$insertSequence)) {
    insSeq = gr$insertSequence
  }
  names(insSeq) = names(gr)
  return(insSeq)
}
#' Returns the sequence encountered at the other side of a breakpoint when traversing
#' across it.
get_partner_anchor_sequence = function(gr, anchor_length, bsgenome) {
  if (length(gr) == 0) {
    return(character(0))
  }
  seq = rep("", length(gr))
  names(seq) = names(gr)
  seqlevelsStyle(gr) = "UCSC"
  in_ref = seqnames(gr) %in% seqnames(bsgenome)
  if(any(!in_ref)) {
    msg = paste("Ignoring contigs not in reference genome: ", paste0(unique(seqnames(gr)[!in_ref])), collapse = " ")
    warning(msg)
  }
  isbp = !is.na(gr$partner)
  partnergr = partner(gr[isbp])
  anchor_gr = GRanges(seqnames=seqnames(partnergr), ranges=IRanges(
    start=ifelse(strand(partnergr)=="+", start(partnergr) - anchor_length[isbp] + 1, start(partnergr)),
    end=ifelse(strand(partnergr)=="+", end(partnergr), end(partnergr) + anchor_length[isbp] - 1)),
    strand=ifelse(strand(partnergr)=="+", "-", "+"))
  seq[isbp & as.logical(in_ref)] = getSeq(bsgenome, anchor_gr[seqnames(anchor_gr) %in% seqnames(bsgenome)], as.character=TRUE)
  return(seq)
}
get_anchor_support_width = function(vcf, gr) {
  # Does not take into account breakpoint interval encoded by the XNX CIGAR elements
  return(GenomicAlignments::cigarWidthAlongReferenceSpace(info(vcf[names(gr)])$SC, N.regions.removed=TRUE) - 1)
}
get_partner_anchor_support_width = function(vcf, gr) {
  result = get_anchor_support_width(vcf, gr)
  names(result) = names(gr)
  isbp = !is.na(gr$partner)
  result[isbp] = result[gr[isbp]$partner]
  return(result)
}
#' Links breakends by their proximity
#'
#' @param bpgr breakpoint GRanges to link
#' @param maxgap maximum distance between breakends
#' @param strand strands to match
#' @param require_segment require that breakends can be directly joined by a DNA segment between the breakend positions
#' @param select which matching breakend pairs to report
#' @param allow_self_intersection report breakpoints where the two breakends satisfy the adjacency matching criteria
linked_by_adjacency = function(
    bpgr,
    maxgap,
   strand=c("opposite"),
   require_segment=FALSE,
   select=c("closest", "bestqual", "unique", "all"),
   link_label="adj",
   allow_self_intersection=FALSE) {
  strand = match.arg(strand)
  select = match.arg(select)
  if (is.null(bpgr$sampleId)) {
    bpgr$sampleId = rep("placeholder", length(bpgr))
  }
  if (is.null(bpgr$beid)) {
    bpgr$beid = rep(NA_character_, length(bpgr))
  }
  hitdf = findOverlaps(bpgr, bpgr, maxgap=maxgap, ignore.strand=TRUE) %>%
    as.data.frame() %>%
    # intra-sample
    filter(bpgr$sampleId[queryHits] == bpgr$sampleId[subjectHits])
  if (strand == "opposite") {
    hitdf = hitdf %>% filter(as.logical(strand(bpgr)[queryHits] != strand(bpgr)[subjectHits]))
  }
  if (require_segment) {
    if (strand != "opposite") {
      stop("NYI: segment between breakends not possible unless breakends have differen orientations")
    }
    # breakpoints point away from each other
    hitdf = hitdf %>%
      filter(as.logical((start(bpgr)[queryHits] <= start(bpgr)[subjectHits] & strand(bpgr)[queryHits] == "-") |
               (start(bpgr)[queryHits] >= start(bpgr)[subjectHits] & strand(bpgr)[queryHits] == "+")))
  }
  hitdf = hitdf %>%
    mutate(
      distance=abs(start(bpgr)[queryHits] - start(bpgr)[subjectHits]),
      qqual=bpgr$QUAL[queryHits],
      squal=bpgr$QUAL[subjectHits]) %>%
    group_by(queryHits) %>%
    mutate(
      is_closest_query = distance == min(distance),
      is_best_qual_query = qqual == min(qqual),
      partners_query = n()) %>%
    group_by(subjectHits) %>%
    mutate(
      is_closest_subject = distance == min(distance),
      is_best_qual_subject = qqual == min(squal),
      partners_subject = n()) %>%
    ungroup() %>%
    mutate(
      is_closest = is_closest_query & is_closest_subject,
      is_best_qual = is_best_qual_query & is_best_qual_subject,
      is_unique = partners_query == 1 & partners_subject == 1) %>%
    dplyr::select(
      -is_closest_query, -is_closest_subject,
      -is_best_qual_query, -is_best_qual_subject,
      -partners_query, -partners_subject)
  if (select == "closest") {
    hitdf = hitdf %>% filter(is_closest)
  } else if (select == "bestqual") {
    hitdf = hitdf %>% filter(is_best_qual_query)
  } else if (select == "unique") {
    hitdf = hitdf %>% filter(is_unique)
  } else if (select == "all") {
    hitdf = hitdf
  } else {
    stop(paste0("'", select, "' is not a valid value for select"))
  }
  if (!allow_self_intersection) {
    # remove self-intersection
    hitdf = hitdf %>% filter(bpgr$partner[queryHits] != names(bpgr)[subjectHits])
  }
  hitdf = hitdf %>%
    filter(queryHits < subjectHits) %>% # remove symmetry
    mutate(linked_by=paste0(link_label, row_number()))
  return(bind_rows(
    hitdf %>% mutate(sampleId=bpgr$sampleId[queryHits],   beid=bpgr$beid[queryHits],   sourceId=names(bpgr)[queryHits]  ) %>% dplyr::select(sampleId, beid, sourceId, linked_by),
    hitdf %>% mutate(sampleId=bpgr$sampleId[subjectHits], beid=bpgr$beid[subjectHits], sourceId=names(bpgr)[subjectHits]) %>% dplyr::select(sampleId, beid, sourceId, linked_by)) %>%
    distinct())
}
#'
#' @param bsgenome BSGenome reference genome used in the VCF
linked_by_equivalent_variants = function(vcf, gr, min_anchor_bases=20, max_per_base_edit_distance=0.1, bsgenome) {
  # go as far as we have support on the other side for when comparing
  anchor_bases = pmax(min_anchor_bases, get_partner_anchor_support_width(vcf, gr))
  similar_calls_df = sequence_common_prefix(gr, anchor_bases=anchor_bases, maxgap=5, bsgenome) %>%
    filter(per_base_edit_distance <= max_per_base_edit_distance) %>%
    dplyr::select(ssourceId, qsourceId) %>%
    mutate(linked_by=paste0("eqv", row_number())) %>%
    gather(sorq, sourceId, ssourceId, qsourceId) %>%
    dplyr::select(-sorq)
  return(similar_calls_df)
}

passes_very_hard_filters = function(filters) {
  fails = rep(FALSE, length(filters))
  for (vhf in gridss.very_hard_filters) {
    fails = fails | str_detect(filters, stringr::fixed(paste0(";", vhf)))
  }
  return(!fails)
}
passes_soft_filters = function(filters) {
  for (softFilter in gridss.soft_filters) {
    filters = str_replace(filters, stringr::fixed(paste0(";", softFilter)), "")
  }
  return(filters == "" | filters == ";" | filters == "PASS")
}


# SVA functions no longer exported by SVA
.elementExtract.List <- function(x, offset=1) {
  lengths <- S4Vectors::elementNROWS(x)
  flat <- BiocGenerics::unlist(x)
  hasValue <- lengths >= offset
  flatOffset <- head(c(1, 1 + cumsum(lengths)), -1) + offset - 1
  flatOffset[!hasValue] <- length(flat) + 1 # out of bounds
  # need to strip XStringSet since that throws an error
  # on out of bounds instead of returning a correctly typed NA
  return(.unXStringSet(flat)[flatOffset])
}
.elementExtract.ANY <- function(x, offset=1) {
  if (is.null(x)) return(x)
  if (is.vector(x)) {
    if (offset==1) return(x)
    return(x[rep(length(x) + 1, length(x))])
  }
  result <- sapply(x, function(r) r[offset], USE.NAMES=FALSE)
  return(result)
}
.elementExtract.XStringSet <- function(x, offset=1) {
  return(.elementExtract.ANY(as.character(x), offset))
}
#' Extracts the element of each element at the given position
#'
#' @param x list-like object
#' @param offset offset of list
#' @export
setGeneric("elementExtract", function(x, offset=1) standardGeneric("elementExtract"))
setMethod("elementExtract", "XStringSet", .elementExtract.XStringSet)
setMethod("elementExtract", "List", .elementExtract.List)
setMethod("elementExtract", "ANY", .elementExtract.ANY)

#' converts an XStringSet to a character
setGeneric(".unXStringSet", function(x) x)
setMethod(".unXStringSet", "XStringSet", function(x) as.character(x))

#' Uses b if a is NULL
#' @export
'%null%' <- function(a, b) {
  if (is.null(a)) return(b)
  return (a)
}

#' vectorised pairwise longest common prefix
#' Returns the length of the longest common prefix for
#' each string pair
#' @export
pairwiseLCPrefix <- function(s1, s2, ignore.case=FALSE) {
  s1 <- as.character(s1)
  s2 <- as.character(s2)
  if (ignore.case) {
    s1 <- toupper(s1)
    s2 <- toupper(s2)
  }
  prefixLength <- rep(0, max(length(s1), length(s2)))
  matchi <- TRUE
  i <- 1
  while (any(matchi)) {
    s1i <- substring(s1, i, i)
    s2i <- substring(s2, i, i)
    matchi <- s1i != "" & s1i == s2i
    prefixLength <- prefixLength + as.integer(matchi)
    i <- i + 1
  }
  return(prefixLength)
}

fix_parid = function(vcf) {
  if ("PARID" %in% row.names(info(header(vcf))) & !("MATEID" %in% row.names(info(header(vcf))))) {
    parid = info(vcf)$PARID
    info(vcf)$PARID = NULL
    parid_header_ordinal = which(row.names(info(header(vcf))) == "PARID")
    row.names(info(header(vcf)))[parid_header_ordinal] = "MATEID"
    info(header(vcf))$Description[parid_header_ordinal] = "ID of mate breakends"
    info(vcf)$MATEID = parid
    write("WARNING: MATEID header not found. Assuming VCF was generated prior to GRIDSS 2.8.0 and rewriting as MATEID.", stderr())
  }
  return(vcf)
}
flatten_mateid = function(vcf) {
  if (!("MATEID" %in% names(info(vcf)))) {
    stop("Missing MATEID")
  }
  mateid = info(vcf)$MATEID
  if (any(elementNROWS(mateid)) > 1) {
    stop("Multiple MATEID for a single record not supported.")
  }
  info(vcf)$MATEID = as.character(mateid)
  return(vcf)
}

