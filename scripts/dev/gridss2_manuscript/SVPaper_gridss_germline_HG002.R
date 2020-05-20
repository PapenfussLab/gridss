library(StructuralVariantAnnotation)
library(cowplot)
library(tidyverse)
library(rtracklayer)
theme_set(theme_bw())

#' Replaces the NA values in a with corresponding values in b
#' @param a,b objects to be tested or coerced.
#' @return The altered object.
'%na%' <- function(a, b) {
  if (is.null(a) || length(a) == 0) return(b)
  if (is.null(b) || length(b) == 0) return(a)
  return(ifelse(is.na(a), b, a))
}

datadir="S:/hartwig/svtoolkit/"


giab_tier1_regions = import(paste0(datadir, "germline/HG002_SVs_Tier1_v0.6.bed"))
giab_tier1_vcf = readVcf(paste0(datadir, "germline/HG002_SVs_Tier1_v0.6.vcf.gz"))
giab_tier1_gr = breakpointRanges(giab_tier1_vcf)
giab_tier1_gr$caller = "GIAB Tier 1 Truth"
giab_tier1_pass_gr = giab_tier1_gr[giab_tier1_gr$FILTER == "PASS"]

score_for_caller = function(caller_name, vcf) {
  qual = switch(
    caller_name,
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

caller_files = list(
  gridss2=paste0(datadir, "germline/gridss_2.9.1_mapq10.vcf"),
  #gridss2_map10=paste0(datadir, "vcfs/gridss_mapq10/HG002_60x/HG002_60x.gridss.vcf"), # turns out to be worse than mapq20
  #svaba=paste0(datadir, "vcfs/svaba/HG002_60x/HG002_60x.svaba.germline.sv.vcf"), # unfiltered handled by PASS logic
  svaba=paste0(datadir, "vcfs/svaba/HG002_60x/HG002_60x.svaba.unfiltered.germline.sv.vcf"),
  breakdancer=paste0(datadir, "germline/breakdancer_1.4.5.vcf"),
  crest=paste0(datadir, "germline/crest.vcf"),
  delly=paste0(datadir, "germline/delly_0.7.6.vcf"),
  gridss1=paste0(datadir, "germline/gridss_1.6.1.vcf"),
  hydra=paste0(datadir, "germline/hydra_master20160129.vcf"),
  manta=paste0(datadir, "germline/manta_1.1.1.vcf"), # TODO update
  pindel=paste0(datadir, "germline/pindel_0.2.5b6.vcf"))

cgrs = list()
for (caller_name in names(caller_files)) {
  if (is.null(cgrs[[caller_name]])) {
    caller_file = caller_files[[caller_name]]
    write(paste("Processing", caller_name), stderr())
    caller_vcf = readVcf(caller_file)
    VariantAnnotation::fixed(caller_vcf)$QUAL = score_for_caller(caller_name, caller_vcf)
    caller_bpgr = breakpointRanges(caller_vcf)
    caller_begr = breakendRanges(caller_vcf)
    caller_gr = c(caller_bpgr, caller_begr)
    caller_gr$caller = caller_name
    cgrs[[caller_name]] = caller_gr
  }
}

calc_roc_pass_all = function(truth_gr, caller_gr, ...) {
  all_calls=calc_roc(truth_gr, caller_gr, ...)
  all_calls$roc$subset = "All calls"
  all_calls$gr$subset = "All calls"
  pass_calls=calc_roc(truth_gr, caller_gr[is.na(caller_gr$FILTER) | caller_gr$FILTER %in% c("PASS", "", ".")], ...)
  pass_calls$roc$subset = "PASS only"
  pass_calls$gr$subset = "PASS only"
  return(list(roc=bind_rows(all_calls$roc, pass_calls$roc), gr=c(all_calls$gr, pass_calls$gr)))
}
calc_roc = function(truth_gr, caller_gr, filter_to_region_gr=NULL, bpmaxgap=100, bemaxgap=100, minsize=50, maxsizedifference=25, filter_region_margin=100, additional_filter=function(gr) { gr } ) {
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
  truth_gr = truth_gr[abs(truth_gr$svLen) >= minsize]
  caller_bpgr = caller_bpgr[!is.na(caller_bpgr$svLen) & abs(caller_bpgr$svLen) >= minsize]
  
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
  return(list(roc=roc, gr=roc_gr, caller_gr=c(caller_bpgr, caller_begr)))
}

# Add breakpoint only variants
bponly=cgrs[["gridss2"]]
bponly$caller=paste0(bponly$caller, " breakpoints only")
bponly = bponly[!is.na(bponly$partner)]
cgrs[[unique(bponly$caller)]] = bponly
# GRIDSS 1.6.1 was run with the then-beta non-standard breakend reporting capability - drop it as that's a GRIDSS2 feature
cgrs[["gridss1"]] = cgrs[["gridss1"]][!is.na(cgrs[["gridss1"]]$partner)]

roclist = lapply(cgrs, function(gr) {calc_roc_pass_all(giab_tier1_pass_gr, gr, giab_tier1_regions)})
combined_rocgr = unlist(GRangesList(lapply(roclist, function(x) x$gr)))
rocdf = bind_rows(lapply(roclist, function(x) x$roc))

# save.image(paste0(datadir, "workspace.RData"))

rocdf %>% group_by(caller) %>% summarise(recall=max(recall))

require("viridis")
ggplot(rocdf %>% mutate(caller=factor(caller, levels=c("gridss2", "gridss2 breakpoints only", "gridss1", "breakdancer", "crest", "delly", "hydra", "manta", "pindel", "svaba")))) +
  aes(x=recall, y=precision, colour=caller, linetype=subset) +
  #scale_shape_manual(values=c(3,4)) +
  geom_line() +
  scale_linetype_manual(values=c("dashed", "solid")) +
  #scale_color_viridis() +
  #scale_color_brewer(palette="Paired") +
  scale_color_manual(values = c("#000000", "#444444", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  #scale_color_brewer(palette = "Dark2") +
  labs(title="HG002 60x\nGIAB Tier 1 truth set")

ggplot(rocdf %>% filter(str_detect(caller, "gridss"))) +
  aes(x=recall, y=precision, colour=caller, shape=subset) +
  geom_point() +
  #scale_color_brewer(palette="Dark2") +
  labs(title="HG002 60x\nGIAB Tier 1 truth set")

ggplot(as.data.frame(combined_rocgr)) +
  aes(x=abs(svLen), fill=tp) +
  geom_histogram() +
  facet_grid(caller ~. , scales="free") +
  scale_x_log10()

ggplot(as.data.frame(combined_rocgr) %>% filter(tp | (fn %na% FALSE), str_detect(caller, "gridss") | str_detect(caller, "manta"))) +
  aes(x=abs(svLen), fill=ifelse(tp, "TP", ifelse(fn, "FN", "FP"))) +
  geom_histogram(bins=100) +
  scale_x_log10() +
  facet_grid(caller ~. , scales="free")


# What calls did we miss:
gridss_gr = bponly
gridss_gr = gridss_gr[gridss_gr$FILTER == "PASS"]
gridss_gr = gridss_gr[names(gridss_gr) %in% gridss_gr$partner]
good_gr = cgrs[["gridss1"]]
good_gr = good_gr[gridss_gr$FILTER == "PASS"]
good_gr = good_gr[names(good_gr) %in% good_gr$partner]

missed_by_gridss_gr = good_gr[countBreakpointOverlaps(good_gr, giab_tier1_pass_gr, maxgap=100) > 0 & countBreakpointOverlaps(good_gr, gridss_gr, maxgap=50) == 0]
missed_by_gridss_gr$behit = overlapsAny(missed_by_gridss_gr, cgrs[["gridss2"]][is.na(cgrs[["gridss2"]]$partner)])
missed_by_gridss_gr$behit350 = overlapsAny(missed_by_gridss_gr, cgrs[["gridss2"]][is.na(cgrs[["gridss2"]]$partner) & cgrs[["gridss2"]]$QUAL >= 350])


tgr = giab_tier1_pass_gr
tgr$hits_gridss = countBreakpointOverlaps(giab_tier1_pass_gr, gridss_gr, maxgap=100)
tgr$hits_manta =  countBreakpointOverlaps(giab_tier1_pass_gr, cgrs[["manta"]], maxgap=100)
tgr$hits_crest =  countBreakpointOverlaps(giab_tier1_pass_gr, cgrs[["crest"]], maxgap=100)

misseddf = missed_by_gridss_gr %>% as.data.frame()
ggplot(misseddf) +
  aes(x=abs(svLen), fill=behit350) +
  geom_histogram() +
  facet_wrap(~svtype) +
  scale_x_log10() +
  labs(title="TPs missed by GRIDSS")

asm_interval_gr = missed_by_gridss_gr
strand(asm_interval_gr) = "*"
asm_interval_gr = disjoin(flank(missed_by_gridss_gr, 4000))
export(asm_interval_gr, "gridss2_asm_interval_gr.bed")







