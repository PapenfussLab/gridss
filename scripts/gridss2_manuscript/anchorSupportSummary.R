source("libbenchmark.R")

query_somatic_structuralVariants = function(dbConnect, table="structuralVariant") {
  query = paste(
    "SELECT id, sampleId, startChromosome, endChromosome, startPosition, endPosition, startOrientation, endOrientation, startHomologySequence, endHomologySequence, insertSequence, startIntervalOffsetStart,startIntervalOffsetEnd, endIntervalOffsetStart,endIntervalOffsetEnd, startLinkedBy,endLinkedBy, startAnchoringSupportDistance, endAnchoringSupportDistance",
    " FROM ", table,
    " WHERE filter = 'PASS'",
    sep = "")
  return (DBI::dbGetQuery(dbConnect, query) %>%
            mutate(id=as.character(id)))
}
sv_gr <- function(dbdf, include.homology=TRUE) {
  grs = GRanges(
    seqnames=dbdf$startChromosome,
    ranges=IRanges(start=dbdf$startPosition + ifelse(include.homology, ifelse(is.na(dbdf$startIntervalOffsetStart), 0, dbdf$startIntervalOffsetStart), 0),
                   end=dbdf$startPosition + ifelse(include.homology, ifelse(is.na(dbdf$startIntervalOffsetEnd), 0, dbdf$startIntervalOffsetEnd), 0)),
    strand=ifelse(dbdf$startOrientation == 1, "+", "-"),
    insertSequence=dbdf$insertSequence,
    partner=ifelse(is.na(dbdf$endChromosome), NA_character_, paste0(dbdf$id, "h")),
    id=dbdf$id,
    sampleid=dbdf$sampleId,
    anchorSupportDistance=dbdf$startAnchoringSupportDistance,
    beid=paste0(dbdf$id, ifelse(is.na(dbdf$endChromosome), "b",  "o")),
    linkedBy=dbdf$startLinkedBy)
  names(grs)=grs$beid
  dbdf = dbdf %>% filter(!is.na(endChromosome))
  rc_insert_sequence = dbdf$insertSequence
  rc_insert_sequence[!str_detect(rc_insert_sequence, "[^ACGTN]")] = as.character(reverseComplement(DNAStringSet(rc_insert_sequence[!str_detect(rc_insert_sequence, "[^ACGTN]")])))
  grh = GRanges(
    seqnames=dbdf$endChromosome,
    ranges=IRanges(start=dbdf$endPosition + ifelse(include.homology, ifelse(is.na(dbdf$endIntervalOffsetStart), 0, dbdf$endIntervalOffsetStart), 0),
                   end=dbdf$endPosition + ifelse(include.homology, ifelse(is.na(dbdf$endIntervalOffsetEnd), 0, dbdf$endIntervalOffsetEnd), 0)),
    strand=ifelse(dbdf$endOrientation == 1, "+", "-"),
    insertSequence=ifelse(dbdf$startOrientation != dbdf$endOrientation, dbdf$insertSequence, rc_insert_sequence),
    partner=paste0(dbdf$id, "o"),
    id=dbdf$id,
    sampleid=dbdf$sampleId,
    anchorSupportDistance=dbdf$endAnchoringSupportDistance,
    beid=paste0(dbdf$id, "h"),
    linkedBy=dbdf$endLinkedBy)
  names(grh)=grh$beid
  return(c(grs, grh))
}
#db = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
#dbdf = query_somatic_structuralVariants(db)
#dbdf = dbdf %>% filter(sampleId %in% read_csv(paste0(privatedatadir, "/hartwig/LNX_SVS.csv"))$SampleId)
#save(dbdf, file = paste0(privatedatadir, "anchorsupport.RData"))
load(file = paste0(privatedatadir, "anchorsupport.RData"))
source("libbenchmark.R")
gr = sv_gr(dbdf)
#pgr = sgr[ifelse(is.na(sgr$partner), names(sgr), sgr$partner)]

asm_linked_breakends = as.data.frame(gr) %>%
  dplyr::select(sampleid, beid, linkedBy) %>%
  dplyr::mutate(linkedBy = str_split(as.character(linkedBy), stringr::fixed(","))) %>%
  tidyr::unnest(linkedBy) %>%
  group_by(sampleid) %>%
  dplyr::filter(str_detect(linkedBy, "^asm"))
asm_links = asm_linked_breakends %>%
  inner_join(asm_linked_breakends, by=c("sampleid", "linkedBy"), suffix = c("1", "2")) %>%
  filter(beid1 != beid2) %>%
  group_by(sampleid, beid1, beid2) %>%
  summarise(linkedBy=paste0(linkedBy, collapse=","))

maxdistance=1000
anchor_df = findOverlaps(gr, gr, maxgap=maxdistance, ignore.strand=TRUE) %>%
  as.data.frame() %>%
  filter(
    gr$id[queryHits] != gr$id[subjectHits],
    gr$sampleid[queryHits] == gr$sampleid[subjectHits],
    as.logical(strand(gr)[queryHits] == "-"),
    as.logical(strand(gr)[subjectHits] == "+"),
    start(gr)[queryHits] <= start(gr)[subjectHits]) %>%
  mutate(
    beid1=gr$beid[queryHits],
    beid2=gr$beid[subjectHits],
    distance=abs(start(gr[queryHits])-start(gr[subjectHits])),
    anchorSupportDistance1=gr$anchorSupportDistance[queryHits],
    anchorSupportDistance2=gr$anchorSupportDistance[subjectHits]) %>%
  dplyr::select(-queryHits, -subjectHits)


adj_df = anchor_df %>%
  group_by(beid1) %>%
  mutate(isClosest=distance==min(distance)) %>%
  ungroup() %>%
  left_join(asm_links, by=c("beid1"="beid1", "beid2"="beid2")) %>%
  mutate(
    is_asm_linked=!is.na(linkedBy),
    isbp1 = !is.na(gr[beid1]$partner),
    isbp2 = !is.na(gr[beid2]$partner)) %>%
  mutate(phasing=ifelse(is_asm_linked, "cis", ifelse(pmax(anchorSupportDistance1, anchorSupportDistance2) > distance + 10, "trans", "unphased")))

ggplot(adj_df %>% filter(is_asm_linked)) +
  aes(x=distance, fill=pmax(anchorSupportDistance1, anchorSupportDistance2) > distance + 10) +
  geom_histogram(bins=50) +
  labs(title="Assembly linked variants", fill="supportDistance > distance")

ggplot() +
  aes(x=distance, y=pmax(anchorSupportDistance1, anchorSupportDistance2), colour=is_asm_linked) +
  geom_point(data=adj_df %>% filter(!is_asm_linked) %>% sample_frac(size=0.1), colour="red", size=0.1) +
  geom_point(data=adj_df %>% filter(is_asm_linked), colour="blue", size=0.1) +
  scale_x_continuous(limits=c(0, 600))


ggplot() +
  aes(x=distance, y=pmax(anchorSupportDistance1, anchorSupportDistance2), colour=is_asm_linked) +
  geom_point(data=adj_df %>% filter(!is_asm_linked) %>% sample_frac(size=0.1), colour="red", size=0.1) +
  geom_point(data=adj_df %>% filter(is_asm_linked), colour="blue", size=0.1) +
  scale_x_continuous(limits=c(0, 600))

ggplot(
    adj_df %>% filter(isClosest & phasing == "unphased") %>%
    mutate(probably_cis=pmax(anchorSupportDistance1, anchorSupportDistance2) > distance - 5)) +
  aes(x=distance, y=pmax(anchorSupportDistance1, anchorSupportDistance2)) +
  geom_point(size=0.1) +
  scale_x_continuous(limits=c(0, 1000)) +
  scale_y_continuous(limits=c(0, 1000)) +
  labs(title="Unphased variants", x="Distance to nearest SV", y="(max) anchor support distance")
figsave("unphased_anchoring_distance", width=8, height=8)
adj_df %>% filter(isClosest & phasing == "unphased") %>%
  mutate(probably_cis=pmax(anchorSupportDistance1, anchorSupportDistance2) > distance - 5) %>%
  pull(probably_cis) %>%
  table()

ggplot(adj_df %>% filter(isClosest | is_asm_linked) %>% mutate(type=ifelse(isbp1 & isbp1, "Both breakpoints", ifelse(isbp1 | isbp2, "Breakpoint/breakend", "Both breakends")))) +
  aes(x=distance, fill=phasing) +
  geom_histogram(binwidth=10, boundary=0) +
  scale_x_continuous(
    expand = c(0,0),
    limits = c(30, 800),
    breaks=c(30, seq(100, 800, 100)),
    labels=c("30", "", "200", "", "400", "", "600", "", "800")) +
  coord_cartesian(xlim=c(30, 800)) +
  facet_wrap( ~ type) +
  scale_y_continuous(expand=expand_scale(mult=c(0, 0.02))) +
  scale_fill_manual(values=c("#b2df8a", "#984ea3", "#ff7f00")) +
  theme(
    axis.line = element_line(colour = "black"),
    panel.border = element_blank()) +
  labs(y="structural variants", x="distance to adjacent SV", title="Phasing by breakpoint/breakend status")
figsave("assembly_phasing_by_type", width=5, height=4)

plot_assembly_phasing = ggplot(adj_df %>% filter(isClosest | is_asm_linked) %>% filter(isbp1 & isbp2)) +
  aes(x=distance, fill=phasing) +
  geom_histogram(binwidth=10, boundary=0) +
  scale_x_continuous(
    expand = c(0,0),
    limits = c(30, 800),
    breaks=c(30, seq(100, 800, 100)),
    labels=c("30", "", "200", "", "400", "", "600", "", "800")) +
  coord_cartesian(xlim=c(30, 800)) +
  scale_y_continuous(expand=expand_scale(mult=c(0, 0.02))) +
  scale_fill_manual(values=c("#b2df8a", "#984ea3", "#ff7f00")) +
  theme(
    axis.line = element_line(colour = "black"),
    panel.border = element_blank()) +
  labs(y="structural variants", x="distance to adjacent SV") +
  theme(legend.position=c(.9,.9))
plot_assembly_phasing
figsave("assembly_phasing_breakpoint", width=5, height=4)

ggplot(adj_df %>% filter(isClosest | is_asm_linked)) +
  aes(x=distance, fill=phasing) +
  geom_histogram(binwidth=10, boundary=0) +
  scale_x_continuous(
    expand = c(0,0),
    limits = c(0, 800),
    breaks=seq(0, 800, 100),
    labels=c("0", "", "200", "", "400", "", "600", "", "800")) +
  coord_cartesian(xlim=c(0, 800)) +
  scale_y_continuous(expand=expand_scale(mult=c(0, 0.02))) +
  scale_fill_manual(values=c("#b2df8a", "#984ea3", "#ff7f00")) +
  theme(
    axis.line = element_line(colour = "black"),
    panel.border = element_blank()) +
  labs(y="structural variants", x="distance to adjacent SV")
figsave("assembly_phasing", width=5, height=4)
# Counts:
adj_df_30_1000 = adj_df %>% filter(distance >= 30)
length(unique(c(adj_df_30_1000$beid1,adj_df_30_1000$beid2)))
gr$hasNearby = gr$beid %in% unique(c(adj_df_30_1000$beid1, adj_df_30_1000$beid2)) | gr$partner %in% unique(c(adj_df_30_1000$beid1, adj_df_30_1000$beid2))
adj_df %>% filter(isClosest) %>% filter(distance > 30) %>% pull(phasing) %>% table()
table(gr$hasNearby)

nearbyBySample = gr %>% as.data.frame() %>%
  group_by(sampleid) %>%
  summarise(n=n(), hasNearby=sum(hasNearby)) %>%
  mutate(portionWithNearby=100 * hasNearby / n)
median(nearbyBySample$portionWithNearby)
ggplot(nearbyBySample) +
  aes(x=portionWithNearby) +
  geom_histogram(bins=50) +
  labs(y="Samples", x="Percentage of SVs with nearby SV", title="Portion of SVs with an adjacent SV within 35-1000bp")
figsave("portion_of_nearby_svs", width=5, height=4)

adj_df %>% group_by(beid1) %>%
  summarise(state=max(ifelse(phasing=="cis", 2, ifelse(phasing=="trans", 1, 0)))) %>%
  group_by(state) %>%
  summarise(n=n()) %>%
  mutate(percentage=n/length(gr))
# TODO: fix percentage as we're:
# - only counting the left side
# - not filtering to > 50bp

nearby_summary_df = gr %>% as.data.frame() %>% dplyr::select(beid, sampleid) %>%
  left_join(rbind(
        adj_df %>% filter(distance >= 50 & distance <= 1000 & str_sub(adj_df$beid1, end=8) != str_sub(adj_df$beid2, end=8)) %>% mutate(beid=beid1) %>% dplyr::select(beid, distance),
        adj_df %>% filter(distance >= 50 & distance <= 1000 & str_sub(adj_df$beid1, end=8) != str_sub(adj_df$beid2, end=8)) %>% mutate(beid=beid2) %>% dplyr::select(beid, distance)) %>%
      group_by(beid) %>% summarise(distance=min(distance)),
  by="beid") %>%
group_by(sampleid) %>%
  summarise(n=n(), hasNearby=sum(!is.na(distance)))

nearby_summary_df = nearby_summary_df %>%
  filter(hasNearby > 1) %>%
  group_by(sampleid, n, hasNearby) %>%
  do({data.frame(
      pvalue=prod(.$n - seq(1, .$hasNearby - 1)) * ((1000/3000000000) ** (.$hasNearby - 1)),
      log10pvalue=sum(log10(.$n - seq(1, .$hasNearby - 1))) + (.$hasNearby - 1) * log10(1000/3000000000)
  )})
sum(nearby_summary_df$log10pvalue)

# process each sample separately to avoid hotspot pairing combinatoric explosion
long_anchor_df =
  data.frame(sampleid=unique(gr$sampleid)) %>%
  group_by(sampleid) %>%
  do({
    samplegr = gr[gr$sampleid==.$sampleid]
    findOverlaps(samplegr, samplegr, maxgap=20000, ignore.strand=TRUE) %>%
    as.data.frame() %>%
    filter(
      samplegr$id[queryHits] != samplegr$id[subjectHits],
      samplegr$sampleid[queryHits] == samplegr$sampleid[subjectHits],
      as.logical(strand(samplegr)[queryHits] == "-"),
      as.logical(strand(samplegr)[subjectHits] == "+"),
      start(samplegr)[queryHits] <= start(samplegr)[subjectHits]) %>%
    mutate(
      beid1=samplegr$beid[queryHits],
      beid2=samplegr$beid[subjectHits],
      distance=abs(start(samplegr[queryHits])-start(samplegr[subjectHits])),
      anchorSupportDistance1=samplegr$anchorSupportDistance[queryHits],
      anchorSupportDistance2=samplegr$anchorSupportDistance[subjectHits]) %>%
    dplyr::select(-queryHits, -subjectHits)
  })

long_anchor_df %>%
  group_by(beid1) %>%
  summarise(
    within500=sum(distance <= 500),
    within20000=n()) %>%
  ungroup() %>%
  summarise(
    has500 = sum(within500 > 0),
    has20000 = sum(within20000 > 0),
    has500percent=has500/length(gr),
    has20000percent=has20000/length(gr),
    longreadimprovementpercent=has20000percent - has500percent)





########
# Theoretical phasability
distance_to_phasable = function(gr, maxdistance=10000000, maxbreakpoints=10) {
  leftgr = gr[as.character(strand(gr)) == "-"]
  rightgr = gr[as.character(strand(gr)) == "+"]
  left_overlap_gr = resize(leftgr, fix="start", width=maxdistance, ignore.strand=TRUE)
  hits = findOverlaps(left_overlap_gr, rightgr, ignore.strand=TRUE)
  resultdf = data.frame(
    sampleid = leftgr$sampleid[queryHits(hits)],
    startid=names(leftgr)[queryHits(hits)],
    endid=names(rightgr)[subjectHits(hits)],
    distance = start(rightgr)[subjectHits(hits)] - start(leftgr)[queryHits(hits)] + 1) %>%
    filter(distance > 0) %>%
    group_by(startid) %>%
    arrange(distance) %>%
    mutate(breakpoint=row_number()) %>%
    filter(breakpoint <= maxbreakpoints)
  return(resultdf)
}
phasabledf = bind_rows(lapply(unique(gr$sampleid), function(s) distance_to_phasable(gr[gr$sampleid == s])))
# Calculate the theoretical distribution
# we have enough data points that we can just simulate by creating breakpoints with total # matching our cohort
hg19_chr_len = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:24]
hg19_chr_offset = c(0, cumsum(as.numeric(hg19_chr_len)))
names(hg19_chr_offset) = names(hg19_chr_len)
generate_random_breaks = function(breakendCount, sampleid) {
  rdf = data.frame(linearpos = runif(breakendCount, 1, sum(hg19_chr_len))) %>%
    mutate(
      strand=ifelse(row_number() %% 2 == 0, "+", "-"),
      chr=as.character(cut(linearpos, hg19_chr_offset, names(hg19_chr_len))),
      pos=linearpos - hg19_chr_offset[chr])
  randgr = with(rdf, GRanges(seqnames=chr, strand=strand, ranges=IRanges(start=pos, width=1)))
  randgr$sampleid=paste0("random_", sampleid)
  names(randgr) = paste0("rand_", sampleid, "_", seq_along(randgr$sampleid))
  return(randgr)
}
phasable_random_baseline_df = bind_rows(lapply(unique(gr$sampleid), function(s) distance_to_phasable(generate_random_breaks(sum(gr$sampleid == s), s))))

phasablesummarydf = bind_rows(
    phasabledf %>% mutate(dataset="actual"),
    phasable_random_baseline_df %>% mutate(dataset="expected")) %>%
  mutate(distance = signif(distance, 3)) %>% # round to 3 significant figures so we don't overplot
  group_by(dataset, breakpoint, distance) %>%
  summarise(count=n()) %>%
  group_by(dataset, breakpoint) %>%
  arrange(distance) %>%
  mutate(cumcount=cumsum(count)) %>%
  mutate(cdf=2 * cumcount / length(gr)) # x2 because phasing two breakpoints makes the other side of the breakpoints also phased

ggplot(phasablesummarydf %>% filter(breakpoint == 1)) +
  aes(x = distance, y=cdf, linetype=dataset) +
  scale_x_log10(limits=c(100, 10000000), expand=c(0,0), breaks=c(100, 1000, 10000, 100000, 1000000, 10000000), labels=c("100", "1kb", "10kb", "100kb", "1Mb", "10Mb")) +
  scale_y_continuous(expand=c(0,0, 0, 0.01), labels = scales::percent_format()) +
  geom_line() +
  labs(x="Distance to next breakpoint", y="Percentage of SVs in cohort")
figsave("hartwig_distance_to_breakpoint", width=3, height=3)

ggplot(phasablesummarydf) +
  aes(x=distance, y=cdf, colour=as.factor(breakpoint), linetype=dataset) +
  scale_x_log10(limits=c(100, 10000000), expand=c(0,0), breaks=c(100, 1000, 10000, 100000, 1000000, 10000000), labels=c("100", "1kb", "10kb", "100kb", "1Mb", "10Mb")) +
  scale_y_continuous(expand=c(0,0, 0, 0.01), labels = scales::percent_format()) +
  geom_line() +
  labs(x="Distance to next breakpoint", y="Percentage of SVs in cohort")

phasablesummarydf %>% filter(breakpoint == 1 & distance == 10000)
    
    
    
    
    
    
    
    
    
    
    
    


