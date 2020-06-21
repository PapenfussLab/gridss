source("libbenchmark.R")

truth = list()
truth$hg002 = list()
truth$hg002$giab_tier1_regions = import(paste0(datadir, "hg002/HG002_SVs_Tier1_v0.6.bed"))
truth$hg002$giab_tier1_vcf = readVcf(paste0(datadir, "hg002/HG002_SVs_Tier1_v0.6.vcf.gz"))
truth$hg002$giab_tier1_gr = breakpointRanges(truth$hg002$giab_tier1_vcf)
truth$hg002$giab_tier1_gr$caller = "GIAB Tier 1 Truth"
truth$hg002$giab_tier1_pass_gr = truth$hg002$giab_tier1_gr[truth$hg002$giab_tier1_gr$FILTER == "PASS"]
truth$hg002$gr = truth$hg002$giab_tier1_pass_gr
truth$hg002$region_gr = truth$hg002$giab_tier1_regions
truth$hg002$filter_function = NULL
truth$na12878 = list()
truth$na12878$parikh_vcf = readVcf(paste0(datadir, "na12878/00000000000000000000000000000003.vcf"))
truth$na12878$parikh_gr = breakpointRanges(truth$na12878$parikh_vcf)
truth$na12878$gr = truth$na12878$parikh_gr
truth$na12878$region_gr = NULL
truth$na12878$filter_function = NULL
seqlevelsStyle(truth$na12878$gr) = "UCSC"

apply_gridss293_filters = function(caller_name, vcf) {
  if (caller_name == "gridss2") {
    i = info(vcf)
    VariantAnnotation::fixed(vcf)$FILTER = paste0(VariantAnnotation::fixed(vcf)$FILTER, ifelse(with(i, VF == 0 & BUM == 0), ";NO_RP", ""))
    VariantAnnotation::fixed(vcf)$FILTER = paste0(VariantAnnotation::fixed(vcf)$FILTER, ifelse(with(i, VF == 0 & BSC == 0), ";NO_SR", ""))
    VariantAnnotation::fixed(vcf)$FILTER = paste0(VariantAnnotation::fixed(vcf)$FILTER, ifelse(with(i, VF == 0 & (BASSR + BASRP - BSC - BUM) / (BASSR + BASRP + BSC + BUM) > 0.5), ";BE_BIAS", ""))
  }
  return(vcf)
}
caller_files = list(
  hg002 = list(
  	gridss=paste0(datadir, "hg002/gridss_2.9.3.vcf"),
    #gridss=paste0(datadir, "hg002/gridss_2.9.1_mapq10.vcf"),
    svaba=paste0(datadir, "svaba/HG002_60x/HG002_60x.svaba.unfiltered.germline.sv.vcf"),
    breakdancer=paste0(datadir, "hg002/breakdancer_1.4.5.vcf"),
    crest=paste0(datadir, "hg002/crest.vcf"),
    delly=paste0(datadir, "hg002/delly_0.7.6.vcf"),
    gridss1=paste0(datadir, "hg002/gridss_1.6.1.vcf"),
    hydra=paste0(datadir, "hg002/hydra_master20160129.vcf"),
    manta=paste0(datadir, "hg002/manta_1.1.1.vcf"),
    pindel=paste0(datadir, "hg002/pindel_0.2.5b6.vcf")),
  na12878 = list(
    gridss=paste0(datadir, "gridss_mapq10/NA12878/NA12878.gridss.vcf"),
    svaba=paste0(datadir, "svaba/na12878/NA12878.svaba.unfiltered.germline.sv.vcf"),
    breakdancer=paste0(datadir, "na12878/breakdancer_1.4.5.vcf"),
    crest=paste0(datadir, "na12878/crest.vcf"),
    delly=paste0(datadir, "na12878/delly_0.7.6.vcf"),
    gridss1=paste0(datadir, "na12878/gridss_1.4.1.vcf"),
    hydra=paste0(datadir, "na12878/hydra_master20160129.vcf"),
    manta=paste0(datadir, "na12878/manta_1.1.1.vcf"))
    #pindel=paste0(datadir, "na12878/pindel_0.2.5b6.vcf")) # terminated early. VCF is incomplete
)
if (!exists("cgrs")) {
  cgrs = list(hg002=list(), na12878=list())
}
for (sample_name in names(caller_files)) {
  for (caller_name in names(caller_files[[sample_name]])) {
    if (is.null(cgrs[[caller_name]])) {
      caller_file = caller_files[[sample_name]][[caller_name]]
      write(paste("Processing", sample_name, caller_name), stderr())
      caller_vcf = readVcf(caller_file)
      #caller_vcf = apply_gridss293_filters(caller_name, caller_vcf)
      VariantAnnotation::fixed(caller_vcf)$QUAL = score_for_caller(caller_name, caller_vcf)
      caller_bpgr = breakpointRanges(caller_vcf)
      caller_begr = breakendRanges(caller_vcf)
      caller_gr = c(caller_bpgr, caller_begr)
      caller_gr$caller = caller_name
      caller_gr$sample_name = sample_name
      cgrs[[sample_name]][[caller_name]] = caller_gr
    }
  }
}

# Add breakpoint only variants
#bponly=cgrs[["gridss"]]
#bponly$caller=paste0(bponly$caller, " breakpoints only")
#bponly = bponly[!is.na(bponly$partner)]
#cgrs[[unique(bponly$caller)]] = bponly

# GRIDSS 1.6.1 was run with the then-beta non-default breakend reporting capability - drop it as that's a GRIDSS2 feature
cgrs$hg002$gridss1 = cgrs$hg002$gridss1[!is.na(cgrs$hg002$gridss1$partner)]
cgrs$na12878$gridss1 = cgrs$hg002$gridss1[!is.na(cgrs$na12878$gridss1$partner)]

roclist = lapply(cgrs, function(s) { lapply(s, function(gr) {
  sample_name = unique(gr$sample_name)
  calc_roc_pass_all(
    truth[[sample_name]]$gr,
    gr,
    sample_name = sample_name,
    filter_to_region_gr=truth[[sample_name]]$region_gr,
    additional_filter=truth[[sample_name]]$filter_function)
})})
combined_rocgr = unlist(GRangesList(lapply(roclist, function(x) unlist(GRangesList(lapply(x, function(xx) xx$gr))))))
rocdf = bind_rows(lapply(roclist, function(x) bind_rows(lapply(x, function(xx) xx$roc))))
roc_by_type = bind_rows(lapply(roclist, function(x) bind_rows(lapply(x, function(xx) xx$roc_by_type))))

# save.image(paste0(datadir, "workspace.RData"))

rocdf %>% group_by(sample_name, caller) %>% 
  #mutate(f1score=2*((precision*recall)/(precision+recall))) %>%
  summarise(recall=max(recall))

require("viridis")
ggplot(rocdf %>% 
    mutate(caller=factor(caller, levels=c("gridss", "gridss1", "breakdancer", "crest", "delly", "hydra", "manta", "pindel", "svaba")))) +
  aes(x=recall, y=precision, colour=caller, linetype=subset) +
  #scale_shape_manual(values=c(3,4)) +
  geom_line() +
	facet_wrap(~ sample_name) +
  scale_linetype_manual(values=c("dashed", "solid")) +
  #scale_color_viridis() +
  #scale_color_brewer(palette="Paired") +
  scale_color_manual(values = c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  #scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(limits = c(0,1), expand=c(0, 0), labels = scales::percent) +
  scale_y_continuous(limits = c(0,1.01), expand=c(0, 0), labels = scales::percent) +
  theme(plot.margin = margin(0,0,0,0),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
figsave("hg002_giab_tier1_roc", width=6, height=5)


ggplot(roctypedf %>% mutate(caller=factor(caller, levels=c("gridss2", "gridss2 breakpoints only", "gridss1", "breakdancer", "crest", "delly", "hydra", "manta", "pindel", "svaba")))) +
  aes(x=recall, y=precision, colour=caller, linetype=subset) +
  #scale_shape_manual(values=c(3,4)) +
  geom_line() +
  scale_linetype_manual(values=c("dashed", "solid")) +
  facet_wrap(~ simpletype) +
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







