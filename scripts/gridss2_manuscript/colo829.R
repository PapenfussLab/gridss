source("libbenchmark.R")

colo829truth = readVcf(paste0(datadir, "COLO829.somatic.overlapped.vcf"))
colo829truthgr = breakpointRanges(colo829truth, inferMissingBreakends=TRUE)
colo829truthgr$FILTER = "PASS"
colo829truthgr$QUAL = 20
regression_callers = list.files(paste0(datadir, "regressiontest"), pattern="^v.*")
pipeline_regression_callers = list.files(paste0(datadir, "regressiontest"), pattern="^[0-9]{6}_.*")
callers = c("gridss", "manta", "novobreak", "svaba")
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
  if (caller %in% callers) {
    filename = paste0(datadir, caller, "/", sample_name, "/", switch(caller,
      "gridss"=paste0(sample_name, ".gridss.vcf.somatic.vcf"),
      "manta"="workspace/svHyGen/somaticSV.0000.vcf",
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
somatics = load_all_somatics(c(callers, regression_callers, pipeline_regression_callers))
#somatics = load_all_somatics(c(regression_callers, pipeline_regression_callers))
rocdf = bind_rows(lapply(names(somatics), function(c) { bind_rows(lapply(samples, function(s) {somatics[[c]][[s]]$roc })) }))
summarydf = rocdf %>% group_by(caller, sample_name, subset) %>%
  filter(cumncalls==max(cumncalls)) %>%
  mutate(f1score=2*((precision*recall)/(precision+recall)))
gr = unlist(GRangesList(unlist(lapply(names(somatics), function(c) { unlist(lapply(samples, function(s) {somatics[[c]][[s]]$gr })) }))))

ggplot(summarydf %>% filter(!str_detect(sample_name, "purity") & subset == "PASS only" & caller %in% callers)) + 
  aes(x=recall, y=precision, colour=caller, label=sample_purity[sample_name]*sample_depth[sample_name]) +
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
  aes(x=recall, y=precision, colour=caller, label=paste0(sample_purity[sample_name]*sample_depth[sample_name], "x")) +
  geom_text_repel() +
  #geom_text(vjust="inward",hjust="inward") +
  geom_point() +
  scale_color_manual(values = c("#000000", "#0072B2", "#d95f02", "#1b9e77")) +
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
  summarydf %>% filter(str_detect(sample_name, "purity") & subset == "PASS only") %>% mutate(dataset="40x/60x purity downsample"),
  summarydf %>% filter(!str_detect(sample_name, "purity") & subset == "PASS only" & caller %in% callers) %>% mutate(dataset="40x/100x replicate")) %>%
  mutate(purity=sample_purity[sample_name])
ggplot(mergeddf) + 
  aes(x=recall, y=precision, colour=caller, shape=dataset) +
  #geom_text_repel() +
  geom_point(aes(size=purity)) +
  geom_line(data=mergeddf %>% filter(dataset=="40x/60x purity downsample")) +
  #scale_color_brewer(palette="Dark2") +
  scale_color_manual(values = c("#000000", "#0072B2", "#d95f02", "#1b9e77")) +
  scale_x_continuous(limits = c(0,1), expand=c(0, 0), labels = scales::percent) +
  scale_y_continuous(limits = c(0,1.01), expand=c(0, 0), labels = scales::percent) +
  scale_shape_manual(values=c(4,0)) +
  theme(plot.margin = margin(0,0,0,0),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
figsave("colo829_combined_roc", width=7, height=4)

ggplot(bind_rows(
  summarydf %>% filter(str_detect(sample_name, "purity")),
  summarydf %>% filter(!str_detect(sample_name, "purity") & caller %in% callers))) + 
  aes(x=recall, y=precision, shape=subset, colour=caller, label=paste0(sample_purity[sample_name]*sample_depth[sample_name], "x")) +
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
  group_by(caller, subset, sample_name) %>%
  filter(cumtp == max(cumtp)) %>%
  group_by(caller, subset) %>%
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

