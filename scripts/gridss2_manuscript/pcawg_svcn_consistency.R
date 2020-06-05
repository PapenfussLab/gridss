source("libbenchmark.R")

#### PCAWG ####
pcawg_dir=paste0(datadir, "pcawgcnsv/")
pcawg_evaluate_cn_transitions = function(sampleId, ...) {
  write(paste("Processing ", sampleId), stderr())
  cnfile=paste0(pcawg_dir, "/", sampleId, ".consensus.20170119.somatic.cna.txt")
  svfile=paste0(pcawg_dir, "/", sampleId, ".pcawg_consensus_1.6.161116.somatic.sv.bedpe")
  if (!file.exists(cnfile)) {
    write(paste("Missing ", cnfile), stderr())
    return(data.frame())
  }
  if (!file.exists(svfile)) {
    write(paste("Missing ", svfile), stderr())
    return(data.frame())
  }
  sv_bedpe = read_delim(svfile, delim="\t", col_names=TRUE, col_types="cnncnncncccc")
  cndf = read_delim(cnfile, delim="\t", col_names=TRUE, col_types="cnnnnnn", na="NA")
  if (nrow(sv_bedpe) == 0 || nrow(cndf) ==0) {
    return(data.frame())
  }
  svgr = with(sv_bedpe, GRanges(
    seqnames=c(chrom1, chrom2),
    ranges=IRanges(start=c(start1 + 1, start2 + 1), end=c(end1, end2)),
    strand=c(strand1, strand2),
    sampleId=sampleId,
    sourceId=c(paste0(sampleId, sv_id, "_o"), paste0(sampleId, sv_id, "_h")),
    partner=c(paste0(sampleId, sv_id, "_h"), paste0(sampleId, sv_id, "_o"))
  ))
  names(svgr) = svgr$sourceId
  cngr = with(cndf, GRanges(
    seqnames=chromosome,
    ranges=IRanges(start=start, end=end),
    sampleId=sampleId,
    cn=total_cn,
    cn_major=major_cn,
    cn_minor=minor_cn,
    star=star))
  cn_consistency = evaluate_cn_transitions(cngr, svgr, ...)
  cn_consistency$sampleId = sampleId
  as.data.frame(cn_consistency)
}

sampleIds = unique(str_extract(list.files(path=pcawg_dir), "^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}"))
pcawg_cn_transitions_list = lapply(sampleIds, pcawg_evaluate_cn_transitions, distance="cn_transition")
pcawg_sv_transitions_list = lapply(sampleIds, pcawg_evaluate_cn_transitions, distance="sv")
# batched version so we don't lose everything when a file is bad
#pcawg_cn_transitions_list = list()
#batch_size = 20
#while (length(sampleIds) - length(pcawg_cn_transitions_list) > 0) {
#  batch_sampleIds = sampleIds[!(sampleIds %in% names(pcawg_cn_transitions_list))]
#  batch_sampleIds = batch_sampleIds[1:min(length(batch_sampleIds), batch_size)]
#  batch_pcawg_evaluate_cn_transitions = lapply(batch_sampleIds, pcawg_evaluate_cn_transitions)
#  names(batch_pcawg_evaluate_cn_transitions) = batch_sampleIds
#  pcawg_cn_transitions_list = c(pcawg_cn_transitions_list, batch_pcawg_evaluate_cn_transitions)
#}
pcawg_cn_transitions = bind_rows(pcawg_cn_transitions_list[unlist(sapply(pcawg_cn_transitions_list, nrow)) > 0])
pcawg_sv_transitions = bind_rows(pcawg_sv_transitions_list[unlist(sapply(pcawg_sv_transitions_list, nrow)) > 0])

#### Hartwig ####
lnx_svs = read_csv(paste0(privatedatadir, "/hartwig/LNX_SVS.csv"))
lnx_cns = read_tsv(paste0(privatedatadir, "/hartwig/LNX_VIS_COPY_NUMBER.tsv"))
hartwig_svgr = lnx_svs %>%
  filter(!(ChrEnd == "0" & is.na(InsertSeq))) %>% # strip purple placeholder calls
  filter(!Recovered) %>% # strip rescued calls since they're not direct GRIDSS output
  lnx_to_gr()
hartwig_cngr = with(lnx_cns, GRanges(
  seqnames=Chromosome,
  ranges=IRanges(start=Start, end=End),
  cn=CopyNumber,
  cn_major=CopyNumber * BAF,
  cn_minor=CopyNumber * (1 - BAF),
  sampleId=SampleId))
hartwig_evaluate_cn_transitions = function(sampleId, ...) {
  write(paste("Processing ", sampleId), stderr())
  svgr = hartwig_svgr[hartwig_svgr$SampleId == sampleId,]
  cngr = hartwig_cngr[hartwig_cngr$sampleId == sampleId,]
  cn_consistency = evaluate_cn_transitions(cngr, svgr, ...)
  cn_consistency$sampleId = sampleId
  as.data.frame(cn_consistency)
}
hartwig_cn_transitions_list = lapply(unique(lnx_svs$SampleId), hartwig_evaluate_cn_transitions)
hartwig_cn_transitions = bind_rows(hartwig_cn_transitions_list[unlist(sapply(hartwig_cn_transitions_list, nrow)) > 0])
hartwig_sv_transitions_list = lapply(unique(lnx_svs$SampleId), hartwig_evaluate_cn_transitions, distance="sv")
hartwig_sv_transitions = bind_rows(hartwig_sv_transitions_list[unlist(sapply(hartwig_sv_transitions_list, nrow)) > 0])

#### Snapshot ####

save.image(paste0(privatedatadir, "pcawg_svcn_consistency.RData"))
#load(paste0(privatedatadir, "pcawg_svcn_consistency.RData"))

#### Analysis ####

cn_transitions = bind_rows(
    pcawg_cn_transitions %>% mutate(cohort="PCAWG"),
    hartwig_cn_transitions %>% mutate(cohort="Hartwig")) %>%
  mutate(
    inGap = overlapsAny(GRanges(seqnames=seqnames, ranges=IRanges(start=start, end=end)), hg19_gaps, maxgap=100000),
    inCentromere = overlapsAny(GRanges(seqnames=seqnames, ranges=IRanges(start=start, end=end)), hg19_centromeres, maxgap=100000),
    nearGapOrCentromere=inGap | inCentromere)

sv_transitions = bind_rows(
  pcawg_sv_transitions %>% mutate(cohort="PCAWG"),
  hartwig_sv_transitions %>% mutate(cohort="Hartwig")) %>%
  mutate(
    inGap = overlapsAny(GRanges(seqnames=seqnames, ranges=IRanges(start=start, end=end)), hg19_gaps, maxgap=100000),
    inCentromere = overlapsAny(GRanges(seqnames=seqnames, ranges=IRanges(start=start, end=end)), hg19_centromeres, maxgap=100000),
    nearGapOrCentromere=inGap | inCentromere,
    distance=ifelse(is.na(distance), 100000, distance))

ggplot(sv_transitions) +
  aes(x=distance, fill=cohort) +
  geom_histogram(position=position_dodge()) +
  scale_x_log10() +
  scale_y_log10()

sv_transitions %>%
  filter(!nearGapOrCentromere) %>%
  mutate(match_cn=distance <= 2 & !is.na(distance)) %>%
  group_by(cohort) %>%
  summarise(
    percent=1-sum(match_cn) / n(),
    match_cn=sum(match_cn),
    n=n())


missed_summary =
  cn_transitions %>%
  filter(!is.na(cn_left) & !is.na(cn_right)) %>%
  mutate(missed_sv = sv_matches == 0) %>%
  group_by(cohort, nearGapOrCentromere) %>%
  summarise(portion_missed=sum(missed_sv) / n(), missed=sum(missed_sv), transitions=n())

missed_clonal_summary =
  cn_transitions %>%
  filter(!is.na(cn_left) & !is.na(cn_right) & abs(cn_left - cn_right) > 0.5) %>%
  mutate(missed_sv = sv_matches == 0) %>%
  group_by(cohort, nearGapOrCentromere) %>%
  summarise(portion_missed=sum(missed_sv) / n(), transitions=n())

cn_transitions %>%
  group_by(cohort) %>%
  filter(!is.na(distance)) %>%
  summarise(percent_matching_position=sum(distance <= 1) / n())

ggplot(cn_transitions %>% filter(can_evaluate_cn_error)) +
  aes(x=pmin(cn_error, 8)) +
  facet_wrap(~ cohort, scales="free") +
  geom_histogram() +
  labs("Magnitude of copy number transition inconsistency")

cn_transitions %>%
  filter(can_evaluate_cn_error) %>%
  filter(!is.na(cn_error)) %>%
  mutate(rounded_cn_error=pmin(round(cn_error), 8)) %>%
  group_by(cohort, rounded_cn_error) %>%
  summarise(count=n()) %>%
  group_by(cohort) %>%
  mutate(percent=1-count/sum(count))

ggplot(cn_transitions %>%
    filter(can_evaluate_cn_error) %>%
    mutate(rounded_cn_error=pmin(round(cn_error), 8)) %>%
    group_by(cohort, rounded_cn_error) %>%
    summarise(count=n()) %>%
    group_by(cohort) %>%
    mutate(percent=1-count/sum(count))) +
  aes(x=rounded_cn_error, y=percent, fill=cohort) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(
    title="Copy number consistency",
    x="Magnitude of copy number transition inconsistency",
    y="Portion of copy number transitions\nwith a single supporting structural variant")

ggplot(cn_transitions %>% filter(!is.na(distance))) +
  aes(x=distance) +
  facet_wrap(~ cohort, scales="free") +
  geom_histogram() +
  scale_x_log10()

cn_transitions %>%
  filter(!is.na(cn_left) & !is.na(cn_right)) %>%
  group_by(cohort, sv_matches, nearGapOrCentromere) %>%
  summarise(n=n()) %>%
  mutate(portion=n/sum(n))

cn_transitions %>%
  filter(!is.na(cn_left) & !is.na(cn_right) & abs(cn_delta) > 0) %>%
  group_by(sv_matches, cohort) %>%
  summarise(n=n()) %>%
  group_by(cohort) %>%
  mutate(portion=n/sum(n))

cn_transitions %>%
  filter(!is.na(cn_left) & !is.na(cn_right) & !nearGapOrCentromere) %>%
  group_by(cohort, sv_matches) %>%
  summarise(n=n()) %>%
  mutate(portion=n/sum(n))

missed_by_sample = cn_transitions %>%
  filter(!is.na(cn_left) & !is.na(cn_right)) %>%
  mutate(missed_sv = sv_matches == 0) %>%
  group_by(sampleId) %>%
  summarise(portion_missed=sum(missed_sv) / n(), transitions=n()) %>%
  arrange(desc(portion_missed))

ggplot(missed_by_sample) +
  aes(x=portion_missed) +
  geom_histogram() +
  labs(title="Portion of CN transitions with no corresponding SV")
ggplot(missed_by_sample) +
  aes(x=portion_missed, y=transitions) +
  geom_point() +
  scale_y_log10() +
  labs(title="Portion of CN transitions with no corresponding SV")

missed_by_sample_excluding_cndelta0 = pcawg_cn_transitions %>%
  filter(!is.na(cn_left) & !is.na(cn_right) & abs(cn_delta) > 0) %>%
  mutate(missed_sv = sv_matches == 0) %>%
  group_by(sampleId) %>%
  summarise(portion_missed=sum(missed_sv) / n(), transitions=n()) %>%
  arrange(desc(portion_missed))
ggplot(missed_by_sample_excluding_cndelta0) +
  aes(x=portion_missed, y=transitions) +
  geom_point() +
  scale_y_log10() +
  labs(title="Portion of CN transitions with no corresponding SV")

ggplot(pcawg_cn_transitions %>% filter(sv_matches)) +
  aes(x=distance, fill="") +
  geom_histogram() +
  scale_y_log10() +
  scale_x_log10()


#### Per-sample breakdown ####

cn_summary_by_sample = cn_transitions %>%
  filter(!is.na(cn_left) & !is.na(cn_right)) %>%
  filter(!inCentromere  & !inGap) %>%
  group_by(cohort, sampleId) %>%
  summarise(
    transitions=n(),
    missing_svs=sum(sv_matches == 0),
    cn_error_evaluable_tranisions=sum(can_evaluate_cn_error),
    mean_cn_error=mean(cn_error, na.rm=TRUE))
# mean CN requires loading the PCAWG CNA tracks jointly instead of individually

ggplot(cn_summary_by_sample %>% arrange(cohort)) +
  aes(x=transitions, y=missing_svs, colour=cohort) +
  geom_point(alpha=0.5) +
  scale_color_brewer(palette="Set2") +
  coord_cartesian(xlim=c(0, 2000), ylim=c(0, 300)) +
  # TODO: add in 2.1% and 21.1% trend lines
  labs(x="Copy number transitions", y="Transitions missing SVs")


#### Split by CN ####
# TODO:
# - bin by integer CN delta
# - bin by floating CN delta, with fixed 100,000 record width (should see more badness at low CN-delta)
ggplot(cn_transitions %>%
         filter(!nearGapOrCentromere) %>%
         mutate(
           classification=factor(ifelse(sv_match_rescued, "Rescued", sv_match_classification), levels=c("Missing SV", "Rescued", "Single Breakend", "Breakpoint")),
           cnbin=round(pmin(abs(cn_delta), 8))) %>%
         group_by(classification, cnbin, cohort) %>%
         summarise(n=n())) +
  aes(x=cnbin, fill=classification, y=n) +
  geom_bar(position="fill", stat="identity") +
  facet_wrap(~cohort)

hartwig_bin_width = 0.25
cnsvmatchdf = cn_transitions %>%
  filter(!nearGapOrCentromere) %>%
  mutate(
    classification=factor(ifelse(sv_match_rescued, "Rescued", sv_match_classification), levels=rev(c("Missing SV", "Rescued", "Single Breakend", "Breakpoint"))))
cnsvmatchdf = bind_rows(
  cnsvmatchdf %>%
    filter(cohort=="PCAWG") %>%
    mutate(
      cnbin=round(pmin(abs(cn_delta), 8)),
      xmin=cnbin-0.5,
      xmax=cnbin+0.5,
      xwidth=xmax-xmin),
  cnsvmatchdf %>%
    filter(cohort=="Hartwig") %>%
    mutate(
      cnbin=round(pmin(abs(cn_delta), 8)/hartwig_bin_width)*hartwig_bin_width,
      xmin=cnbin-hartwig_bin_width/2,
      xmax=cnbin+hartwig_bin_width/2,
      xwidth=xmax-xmin))

cn_transistions_missing_breakpoints_by_cndelta_plot_df = cnsvmatchdf %>%
  group_by(cohort, cnbin, xmin, xmax, xwidth, classification) %>%
  summarise(n=n()) %>%
  group_by(cohort, cnbin, xmin, xmax, xwidth) %>%
  mutate(
    total_n = sum(n),
    percentage=n/total_n) %>%
  ungroup() %>%
  mutate(cohort=ifelse(cohort=="Hartwig", "Hartwig cohort\nGRIDSS/PURPLE", "PCAWG cohort\nConsensus CN/SV"))

ggplot(cn_transistions_missing_breakpoints_by_cndelta_plot_df) +
  aes(x=cnbin, fill=classification, y=percentage, width=xwidth) +
  scale_x_continuous(expand=c(0,0), breaks = 0:8, labels = 0:8) +
  scale_y_continuous(expand=c(0,0, 0, 0.01), breaks=seq(0, 0.3, 0.05), labels=paste0(seq(0, 0.3, 0.05)*100, "%")) +
  coord_cartesian(ylim=c(0, 0.30), xlim=c(0, 8)) +
  geom_col(position="stack") +
  scale_fill_manual(values=c("#DDEEDD", "#1b9e77", "#fb6a4a")) +
  facet_wrap(~cohort) +
  labs(x="Magnitude of copy number change", y="Portion of copy number transitions") +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
figsave("cn_transistions_missing_breakpoints_by_cndelta_percentage", width=5, height=4)

ggplot(cn_transistions_missing_breakpoints_by_cndelta_plot_df) +
  aes(x=cnbin, fill=classification, weight=n, width=xwidth) +
  geom_bar(position="stack") + 
  scale_fill_manual(values=c("#DDEEDD", "#1b9e77", "#fb6a4a")) +
  scale_x_continuous(expand=c(0,0), breaks = 0:8, labels = 0:8) +
  scale_y_continuous(expand=c(0,0, 0, 0.01), breaks=seq(0, 800000, 100000), labels=format(seq(0, 800000, 100000), big.mark=",", scientific=FALSE)) +
  coord_cartesian(xlim=c(0, 8)) +
  facet_wrap(~cohort) +
  labs(x="Magnitude of copy number change", y="Copy number transitions") +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
figsave("cn_transistions_missing_breakpoints_by_cndelta_counts", width=5, height=4)


# Numbers in paper
with(cn_transitions %>% filter(!nearGapOrCentromere), table(cohort, sv_match_classification))
153376/(1954871+67703+153376)
cn_transitions %>%
  filter(!nearGapOrCentromere) %>%
  filter((cohort == "Hartwig" & cn_delta < 0.5) | (cohort == "PCAWG" & cn_delta > 8)) %>%
  group_by(cohort, sv_match_classification) %>%
  summarise(n=n()) %>%
  group_by(cohort) %>%
  mutate(percent=n/sum(n))
