setwd("../")
source("libgridss.R")
setwd("gridss2_manuscript")
source("libbenchmark.R")
#for f in $(find . -name '*.bgz') ; do gunzip -c $f > $(dirname $f)/$(basename $f .bgz) ; done
#zip colo829_benchmark_results.zip $(find . -name '*.gridss.vcf.somatic.vcf') $(find . -name 'somaticSV.0000.vcf') $(find . -name '*.svaba.somatic.sv.vcf') $(find . -name 'novoBreak.pass.flt.vcf') $(find . -name 'somatic.vcf.bgz')
#colo829truth = readVcf(paste0(datadir, "COLO829.somatic.overlapped.vcf"))
#git clone https://github.com/UMCUGenetics/COLO829_somaticSV
colo829truth = VariantAnnotation::readVcf(paste0(datadir, "COLO829_somaticSV/truthset_somaticSVs_COLO829.vcf"))
info(colo829truth)$SVTYPE[!str_detect(as.character(rowRanges(colo829truth)$ALT), stringr::fixed("<"))] = "BND"
colo829truthgr = breakpointRanges(colo829truth)
colo829truthgr$FILTER = "PASS"
colo829truthgr$QUAL = 20
regression_callers = list.files(paste0(datadir, "regressiontest"), pattern="^v.*")
pipeline_regression_callers = list.files(paste0(datadir, "regressiontest"), pattern="^[0-9]{6}_.*")
callers = c("gridss", "manta", "novobreak", "svaba", "gridss1")
samples = c(
  "colo829_1",
  "colo829_2",
  "colo829_3",
  "colo829_purity_5x_of_60x",
  "colo829_purity_10x_of_60x",
  "colo829_purity_15x_of_60x",
  "colo829_purity_20x_of_60x",
  "colo829_purity_25x_of_60x",
  "colo829_purity_30x_of_60x",
  "colo829_purity_45x_of_60x",
  "colo829_purity_50x_of_60x",
  "colo829_purity_60x_of_60x")
sample_depth=c(
  "colo829_1"=100,
  "colo829_2"=100,
  "colo829_3"=100,
  "colo829_purity_5x_of_60x"=60,
  "colo829_purity_10x_of_60x"=60,
  "colo829_purity_15x_of_60x"=60,
  "colo829_purity_20x_of_60x"=60,
  "colo829_purity_25x_of_60x"=60,
  "colo829_purity_30x_of_60x"=60,
  "colo829_purity_45x_of_60x"=60,
  "colo829_purity_50x_of_60x"=60,
  "colo829_purity_60x_of_60x"=60)
sample_purity=c(
  "colo829_1"=1,
  "colo829_2"=1,
  "colo829_3"=1,
  "colo829_purity_5x_of_60x"=5/60,
  "colo829_purity_10x_of_60x"=10/60,
  "colo829_purity_15x_of_60x"=15/60,
  "colo829_purity_20x_of_60x"=20/60,
  "colo829_purity_25x_of_60x"=25/60,
  "colo829_purity_30x_of_60x"=30/60,
  "colo829_purity_45x_of_60x"=45/60,
  "colo829_purity_50x_of_60x"=50/60,
  "colo829_purity_60x_of_60x"=60/60)
load_somatic = function(caller, sample_name) {
  write(paste("Processing", caller, sample_name), stderr())
  if (caller %in% callers) {
    filename = paste0(datadir, caller, "/", sample_name, "/", switch(caller,
      "gridss"=paste0(sample_name, ".gridss.somatic.vcf.bgz"),
      "gridss1"=paste0(sample_name, ".gridss.vcf.somatic.vcf"),
      "manta"="results/variants/somaticSV.vcf.gz",
      "svaba"=paste0(sample_name,".svaba.somatic.sv.vcf"),
      "novobreak"="novoBreak.pass.flt.vcf",
      "somatic.vcf.bgz" # GRIDSS regression
      ))
  }
  if (caller %in% regression_callers) {
    filename = switch(sample_name,
      "colo829_1"= paste0(datadir, "regressiontest/", caller, "/colo829/somatic.vcf.bgz"),
      "missing")
  }
  if (caller %in% pipeline_regression_callers) {
    hmfsamplename = switch(sample_name,
           "colo829_1"="COLO829v001",
           "colo829_2"="COLO829v002",
           "colo829_3"="COLO829v003",
           sample_name)
    filename = paste0(datadir, "regressiontest/", caller, "/structuralVariants/gridss/", hmfsamplename, "R_", hmfsamplename, "T/", hmfsamplename, "T.gridss.somatic.vcf.gz")
    if (!file.exists(filename)) {
      # HMF pipeline location changed in v5.8
      filename = paste0(datadir, "regressiontest/", caller, "/structural_caller/", hmfsamplename, "T.gridss.somatic.filtered.vcf.gz")
    }
    if (!file.exists(filename)) {
      # HMF pipeline location changed v5
      filename = paste0(datadir, "regressiontest/", caller, "/structural_caller/", hmfsamplename, "T.gridss.somatic.vcf.gz")
    }
  }
  if (file.exists(filename)) {
    vcf = readVcf(filename)
    if (caller == "novobreak") {
      # https://github.com/KChen-lab/novoBreak/issues/1
      names(vcf) = paste0("id", seq_along(names(vcf)), "_", seqnames(rowRanges(vcf)), "_", start(rowRanges(vcf)))
    }
    if (caller == "gridss") {
      # Remove the PON filter field since we explicitly handle that in this analysis
      VariantAnnotation::fixed(vcf)$FILTER = ifelse(VariantAnnotation::fixed(vcf)$FILTER == "PON", "PASS", VariantAnnotation::fixed(vcf)$FILTER)
    }
    gr = c(
      breakpointRanges(vcf),
      breakendRanges(vcf))
    if (caller == "novobreak") {
      # https://github.com/KChen-lab/novoBreak/issues/1
      gr$svLen = ifelse(gr$svLen == 0, NA, gr$svLen)
    }
    gr$caller = caller
    gr$sample_name = sample_name
    roc = calc_roc_pass_all(colo829truthgr, gr, sample_name=sample_name, bpmaxgap=100, bemaxgap=100, minsize=50, maxsizedifference=25, keepInterchromosomal=TRUE)
    pongr = pon_filter(gr)
    rocpon = calc_roc_pass_all(colo829truthgr, pongr, sample_name=sample_name, bpmaxgap=100, bemaxgap=100, minsize=50, maxsizedifference=25, keepInterchromosomal=TRUE)
    roc$roc$ponfiltered = FALSE
    roc$gr$ponfiltered = FALSE
    roc$roc_by_type$ponfiltered = FALSE
    rocpon$roc$ponfiltered = TRUE
    rocpon$gr$ponfiltered = TRUE
    rocpon$roc_by_type$ponfiltered = TRUE
    roc = list(
      roc=bind_rows(roc$roc, rocpon$roc),
      gr=c(roc$gr, rocpon$gr),
      roc_by_type=bind_rows(roc$roc_by_type, rocpon$roc_by_type)
    )
  } else {
    warning(paste("Missing", filename))
    return(NULL)
  }
  return(roc)
}
load_all_somatics = function(caller_list = callers, sample_list = samples) {
  out = list()
  for (caller in caller_list) {
    out[[caller]] = list()
    for (s in sample_list) {
      out[[caller]][[s]] = load_somatic(caller, s)
    }
  }
  return(out)
}
somatics = load_all_somatics(c(callers))
#somatics = load_all_somatics(c(callers, regression_callers, pipeline_regression_callers))
#somatics = load_all_somatics(c(regression_callers, pipeline_regression_callers))
rocdf = bind_rows(lapply(names(somatics), function(c) { bind_rows(lapply(samples, function(s) {somatics[[c]][[s]]$roc })) }))
summarydf = rocdf %>% group_by(caller, sample_name, subset, ponfiltered) %>%
  filter(cumncalls==max(cumncalls)) %>%
  mutate(f1score=2*((precision*recall)/(precision+recall)))
gr = unlist(GRangesList(unlist(lapply(names(somatics), function(c) { unlist(lapply(samples, function(s) {somatics[[c]][[s]]$gr })) }))))

ggplot(summarydf %>% filter(!str_detect(sample_name, "purity") & subset == "PASS only" & caller %in% callers)) + 
  aes(x=recall, y=precision, colour=caller, label=sample_purity[sample_name]*sample_depth[sample_name], shape=ponfiltered) +
  geom_text_repel() +
  #geom_text(vjust="inward",hjust="inward") +
  geom_point() +
  scale_color_brewer(palette="Dark2") +
  scale_x_continuous(limits = c(0,1), expand=c(0, 0), labels = scales::percent) +
  scale_y_continuous(limits = c(0,1.01), expand=c(0, 0), labels = scales::percent) +
  theme(plot.margin = margin(0,0,0,0),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggplot(summarydf %>% filter(str_detect(sample_name, "purity")) %>% filter(subset == "PASS only")) + 
  aes(x=recall, y=precision, colour=caller, label=paste0(sample_purity[sample_name]*sample_depth[sample_name], "x"), shape=ponfiltered) +
  geom_text_repel() +
  #geom_text(vjust="inward",hjust="inward") +
  geom_point() +
  scale_color_manual(values = c("#000000", "#0072B2", "#d95f02", "#1b9e77", "#e7298a")) +
  #scale_color_brewer(palette="Dark2") +
  scale_x_continuous(limits = c(0,1), expand=c(0, 0), labels = scales::percent) +
  scale_y_continuous(limits = c(0,1.01), expand=c(0, 0), labels = scales::percent) +
  theme(plot.margin = margin(0,0,0,0),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
figsave("colo829_purity_downsampling", width=6, height=5)

mergeddf = bind_rows(
  summarydf %>%
    filter(str_detect(sample_name, "purity") & subset == "PASS only") %>%
    mutate(dataset="40x/60x purity downsample"),
  summarydf %>%
    filter(!str_detect(sample_name, "purity") & subset == "PASS only" & caller %in% callers) %>%
    mutate(dataset="40x/100x replicate")) %>%
  mutate(purity=sample_purity[sample_name]) %>%
  mutate(PON=ifelse(ponfiltered, "With Panel of Normals", "Without Panel of Normals")) %>%
  filter(sample_name != "colo829_3")
roc_plot = ggplot(mergeddf ) + 
  aes(x=recall, y=precision, colour=caller, shape=dataset, alpha=ponfiltered) +
  #geom_text_repel() +
  geom_point(aes(size=purity)) +
  geom_line(data=mergeddf %>% filter(dataset=="40x/60x purity downsample")) +
  #scale_color_brewer(palette="Dark2") +
  scale_color_manual(values = c("#000000", "#0072B2", "#d95f02", "#1b9e77", "#e7298a")) +
  scale_x_continuous(limits = c(0,1.02), expand=c(0, 0), labels = scales::percent) +
  scale_y_continuous(limits = c(0,1.02), expand=c(0, 0), labels = scales::percent) +
  scale_shape_manual(values=c(1,16)) +
  scale_alpha_manual(values=c(0.5, 1)) +
  theme(plot.margin = margin(0,0,0,0),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
roc_plot
figsave("colo829_combined_roc_overlay", width=7, height=4)
roc_plot + facet_grid(PON ~ . )
figsave("colo829_combined_roc", width=7, height=8)

ggplot(bind_rows(
  summarydf %>% filter(str_detect(sample_name, "purity")),
  summarydf %>% filter(!str_detect(sample_name, "purity") & caller %in% callers))) + 
  aes(x=recall, y=precision, shape=subset, colour=caller, label=paste0(sample_purity[sample_name]*sample_depth[sample_name], "x")) +
  facet_wrap(~ ponfiltered) +
  geom_text_repel() +
  #geom_text(vjust="inward",hjust="inward") +
  geom_point() +
  scale_color_brewer(palette="Dark2") +
  scale_x_continuous(limits = c(0,1), expand=c(0, 0), labels = scales::percent) +
  scale_y_continuous(limits = c(0,1.01), expand=c(0, 0), labels = scales::percent) +
  theme(plot.margin = margin(0,0,0,0),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title="Effect of GRIDSS PON")

# Paper numbers
mergeddf %>% filter(dataset=="40x/100x replicate") %>%
  group_by(ponfiltered, caller, subset, sample_name) %>%
  filter(cumtp == max(cumtp)) %>%
  group_by(ponfiltered, caller, subset) %>%
  summarise(
    precision=mean(precision),
    recall=mean(recall))

## COLO829 regression tests
ggplot(summarydf %>% filter(sample_depth[sample_name] == 100)) +
  aes(x=recall, y=precision, colour=caller, label=str_replace(caller, "^[0-9]{6}_COLO829v00[1-3]_", "")) +
  facet_grid(sample_name ~ subset) +
  geom_point(alpha=0.5) +
  #geom_jitter(width=0.01) +
  #geom_label() +
  geom_text(alpha=0.5) +
  labs(title="COLO829")

summarydf %>% filter(subset=="PASS only") %>% arrange(desc(f1score)) %>% View()

fngr = gr[gr$sample_name %in% samples[1:3] & gr$caller == "gridss" & (gr$fn %na% FALSE) & overlapsAny(gr, gr[gr$tp & gr$caller != "gridss"])]
names(fngr) = NULL
View(fngr %>% as.data.frame() %>% arrange(seqnames, start))
missedgr = unlist(GRangesList(lapply(samples, function(s) { 
  sgr = gr[!(gr$fn %na% FALSE) & gr$sample_name==s & s %in% gr$sample_name[gr$caller == "gridss"]]
  sgr[sgr$tp & !overlapsAny(sgr, sgr[sgr$subset == "PASS only" & sgr$caller=="gridss"])]})))
names(missedgr) = NULL
View(missedgr %>% as.data.frame() %>% arrange(seqnames, start))

v48gr2 = gr[gr$caller=="190513_COLO829v002_v4.8"]
gridss_gr2 = gr[gr$caller=="gridss" & gr$sample_name=="colo829_2"]
missed_gr = v48gr2[v48gr2$tp & !overlapsAny(v48gr2, gridss_gr2[gridss_gr2$FILTER == "PASS"])]


# COLO829 direct regression
grref = somatics$v2.7.3$colo829_1$gr
grtest = somatics$v2.9.4$colo829_1$gr
colo829_missing = grref[!overlapsAny(grref, grtest)]
colo829_additional = grtest[!overlapsAny(grtest, grref)]
colo928_filtered_truth = grtest[grtest$subset=="All calls" & !overlapsAny(grtest, grtest[grtest$subset=="PASS only"]) & overlapsAny(grtest, grref[grref$subset=="PASS only"])]

