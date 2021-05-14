source("libbenchmark.R")
library(chromPlot)
dataDir = 'protecteddata/centromeres'

sglDestChrs = read.csv(paste(dataDir,'sgl_destination_chromosomes.csv',sep='/'))
sglLocalVsTrans = read.csv(paste(dataDir,'sgl_local_vs_trans.csv',sep='/'))
sglOrigPositions = read.csv(paste(dataDir,'sgl_originating_positions.csv',sep='/'))
sglFbCentroSummary = read.csv(paste(dataDir,'sgl_arm_cn_gain_summary.csv',sep='/'))
sglRateOfGain = read.csv(paste(dataDir,'sgl_rate_of_gain.csv',sep='/'))
sgl_centro_sv_data = read.csv(paste(dataDir,'sgl_centro_sv_data.csv',sep='/'))


# additional annotations
sgl_centro_sv_data = sgl_centro_sv_data %>%
   group_by(SampleId, CentroChr) %>%
   mutate(
      breaksToCentromere=n()) %>%
   ungroup() %>%
   mutate(
      hasCNChange=CentroCNChg > 0 & !is.na(CentroCNChg),
      isIntraChromosomal=ChrStart==CentroChr,
      hasMultipleBreaksToCentromere=breaksToCentromere > 1,
      centroArm=paste0(CentroChr, ifelse(PCentroGain, "p", ifelse(QCentroGain, "q", "?"))))




# Plot 1: Destination of centromeric SGLs
#write.csv(sglDestChrs,'protecteddata/centromeres/sgl_destination_chromosomes.csv',row.names = F,quote = F)


print(ggplot(sglDestChrs, aes(x=reorder(Chromosome,ChrIndex), y=n))
      + geom_bar(stat = "identity", colour = "black")
      + labs(x='Chromosome', y='# Centromeric SGLs', title='Desintation for Centromeric SGLs'))


# Plot 2: Frequency of SVs ending on each chromosome and arm by whether a local or translocation SV
#write.csv(sglLocalVsTrans,'protecteddata/centromeres/sgl_local_vs_trans.csv',row.names = F,quote = F)


print(ggplot(sglLocalVsTrans, aes(x=reorder(CentroChrText,-Total), y=Percent, fill=SvType))
      + geom_bar(stat = "identity", colour = "black")
      + labs(x='Centromere Chromosome (total SGLs)', y='% of Centromeric SGLs per Chromosome',title='Proportion of Local vs Translocation Centromeric SGLs by Chromosome')
      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10)))



# Plot 3: Location of local SGLs by chromosome
#write.csv(sglOrigPositions,'protecteddata/centromeres/sgl_originating_positions.csv',row.names = F,quote = F)


centromerebounddf = range(hg19_centromeres) %>%
   as.data.frame() %>%
   mutate(
      StartBucket=start / 1000000,
      EndBucket=end / 1000000,
      Chromosome=seqnames)
centromeredf = hg19_centromeres %>% as.data.frame() %>%
   mutate(
      StartPosBucket=start / 1000000,
      EndPosBucket=end / 1000000,
      Chromosome=seqnames)
gapsdf = hg19_gaps %>% as.data.frame() %>%
   filter(seqnames %in% c(1:22, "X", "Y")) %>%
   mutate(
      StartPosBucket=start / 1000000,
      EndPosBucket=end / 1000000,
      Chromosome=seqnames)
print(ggplot(sglOrigPositions, aes(x=reorder(SourcePosLabel,PosBucket), y=n))
      + geom_bar(stat = "identity", colour = "black")
      + geom_rect(data=gapsdf, aes(x=StartPosBucket, xmin=StartPosBucket, xmax=EndPosBucket, y=-1, ymin=-5, ymax=-3))
      + facet_wrap(~Chromosome)
      + labs(x='Chromosomal Position', y='# SGLs', title = "Originating position of Local Centromeric SGLs by Chromosome")
      + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()))
#+ geom_rect(data=centromeredf, aes(x=StartPosBucket, xmin=StartPosBucket, xmax=EndPosBucket, y=-1, ymin=-3, ymax=-1))
#+ geom_vline(data=centromeredf, aes(xintercept=StartBucket), colour="grey")
#+ geom_vline(data=centromeredf, aes(xintercept=EndBucket), colour="grey")


# Plot 4: Centromeric SGLs vs Arms with CN Gain
#write.csv(sglFbCentroSummary,'protecteddata/centromeres/sgl_arm_cn_gain_summary.csv',row.names = F,quote = F)


print(ggplot(sglFbCentroSummary %>% filter(CentroGain>0|SglCentro>0|FBArms>0), aes(x=CentroGain,y=SglCentro))
      + geom_point(position="jitter")
      + geom_smooth(,method=lm,se=FALSE, fullrange=F)
      + labs(x='Arms with Centromeric Copy Number Change', y='Arms with Centromeric SGLs',title = "Chromosomal Arms per Sample with Centromeric Copy Number Change vs Centromeric SGLs"))

# Plot 5: Rate of CN Gain vs Centromeric SGLs


print(ggplot(sglRateOfGain, aes(x=reorder(Chromosome,ChrIndex),y=Rate,fill=RateOfCentromericGain))
      + geom_bar(position="dodge",stat='identity')
      + labs(x='Centromere Chromosome', y='% of Samples with Centromeric Gain',title='Rate of Centromeric CN Gain with and without Centromeric SGLs'))

ggplot(sglRateOfGain) +
   aes(x=reorder(Chromosome,ChrIndex), y=Rate, fill=RateOfCentromericGain) +
   geom_bar(position="stack", stat='identity') +
   scale_fill_brewer(palette="Dark2") +
   labs(
      x='Chromosome',
      y='% of Samples with Centromeric Gain',
      fill="",
      title='NOT SAME DENOMINATOR - CANNOT STACK\nRate of Centromeric CN Gain with and without Centromeric SGLs')




grBuckets = with(sglOrigPositions, GRanges(seqnames=Chromosome, ranges=IRanges(start=1000000*PosBucket-999999, width=1000000), n=n))
grNBuckets = with(sglOrigPositions %>% group_by(Chromosome, PosBucket) %>% do({data.frame(hitno=seq_len(.$n))}),
                  GRanges(seqnames=Chromosome, ranges=IRanges(start=1000000*PosBucket-999999, width=1000000)))
hg19_gaps$Colors="Black"
pdf(paste0(figdir, "/raw_centromeric_single_location.pdf"))
chromPlot(bands=hg19_gaps, annot1=grNBuckets, figCols=24)
dev.off()

#######################
## Number inline in paper
inner_join(sglFbCentroSummary, sgl_centro_sv_data %>% filter(PCentroGain | QCentroGain) %>% group_by(SampleId) %>% summarise(singleWithCNChange=length(unique(CentroChr)))) %>%
   summarise(cnChange=sum(CentroGain), singleWithCNChange=sum(singleWithCNChange))
sgl_centro_sv_data %>%
   mutate(label=paste(ifelse(hasCNChange, "with CN change", "without CN change"), ifelse(hasMultipleBreaksToCentromere, "multipleBreakends", "singleBreakend"))) %>%
   group_by(label) %>%
   summarise(n=n()) %>%
   ungroup() %>%
   mutate(pct=n/sum(n)*100)
sgl_centro_sv_data %>%
   ungroup() %>%
   filter(!hasCNChange) %>%
   summarise(
      is13_14_15_21_22=sum(ifelse(CentroChr %in% c(13, 14, 15, 21, 22), 1, 0)),
      total=n())

sgl_centro_sv_data %>% 
   filter(CentroChr ==1, ChrStart != 1) %>%
   group_by(CentroChr, centroArm) %>%
   summarise(n=n())
# samples with 1q gain
sgl_centro_sv_data %>% filter(centroArm == "1q") %>% summarise(n=length(unique(SampleId)))


ggplot(sgl_centro_sv_data %>%
      mutate(
         label=paste(ifelse(isIntraChromosomal, "Intra-chromosomal", "Inter-chromosomal"), ifelse(hasCNChange, " with CN change", " without CN change"))
         ) %>%
      group_by(ChrStart, hasCNChange, isIntraChromosomal, label) %>%
      summarise(events=n()) %>%
      group_by(ChrStart) %>%
      mutate(
         total_events=sum(events),
         portion=events/total_events)) +
   aes(x=factor(ChrStart, levels=c(1:22, "X", "Y")), y=portion, fill=label) +
   geom_bar(position="stack", stat="identity") +
   scale_fill_brewer(palette="Paired") +
   scale_y_continuous(limits = c(0,1), expand=c(0,0), labels = scales::percent) +
   labs(x="Originating chromosome", y="Portion of single breakends to centromeric sequence", title="Single breakends to centromeric sequence") +
   theme(plot.margin = margin(0,0,0,0),
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank())

ggplot(sgl_centro_sv_data %>% filter(hasCNChange)) +
   aes(x=factor(CentroChr, levels=c(1:22, "X", "Y")), fill=ifelse(isIntraChromosomal, "Intra-chromosomal", "Inter-chromosomal")) +
   geom_bar(position="stack", stat="count") +
   labs(x="Most likely centromere", y="single breakends", fill="") +
   scale_fill_manual(values=c("#8da0cb", "#66c2a5")) # brewer Set2

ggplot(sgl_centro_sv_data) +
   aes(x=factor(CentroChr, levels=c(1:22, "X", "Y")), fill=paste(ifelse(isIntraChromosomal, "Intrachromosomal", "Interchromosomal"), ifelse(hasCNChange, " with CN change", paste(" without CN change", ifelse(hasMultipleBreaksToCentromere, "multiple breaks", "single break"))))) +
   geom_bar(position="stack", stat="count") +
   labs(x="Most likely centromere", y="single breakends", fill="") +
   scale_fill_manual(values=c("#a6cee3", "#1f78b4", "#124466", "#b2df8a", "#33a02c", "#21661c"))  +# Paired, extended to 3 per colour
   #scale_fill_brewer(palette="Paired") +
   theme(plot.margin = margin(0,0,0,0)) +
   scale_y_continuous(expand=c(0,0, 0.05, 0))
figsave("fig3b_single_breakend_chr_summary", width=8, height=4)
#
sgl_centro_sv_data


ggplot(sgl_centro_sv_data %>%
    mutate(
       hasCNChange=CentroCNChg > 0 & !is.na(CentroCNChg),
       isIntraChromosomal=ChrStart==CentroChr,
       label=paste(ifelse(isIntraChromosomal, "Intrachromosomal", "Interchromosomal"), ifelse(hasCNChange, " with CN change", " without CN change"))
    ) %>%
    group_by(CentroChr, hasCNChange, isIntraChromosomal, label) %>%
    summarise(events=n()) %>%
    group_by(CentroChr) %>%
    mutate(
       total_events=sum(events),
       portion=events/total_events)) +
   aes(x=factor(CentroChr, levels=c(1:22, "X", "Y")), y=portion, fill=label) +
   geom_bar(position="stack", stat="identity") +
   scale_fill_brewer(palette="Paired") +
   scale_y_continuous(limits = c(0,1), expand=c(0,0), labels = scales::percent) +
   labs(x="Centromeric chromosome", y="Portion of single breakends to centromeric sequence", title="Single breakends to centromeric sequence") +
   theme(plot.margin = margin(0,0,0,0),
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank())

chrconfheatmapdf = sgl_centro_sv_data %>%
   group_by(ChrStart, CentroChr, HighConf) %>%
   mutate(hasCNChange=CentroCNChg > 0 & !is.na(CentroCNChg)) %>%
   summarise(nCNChange=sum(hasCNChange), n=n()) %>%
   mutate(cnChangeRate=nCNChange/n)
chrheatmapdf = sgl_centro_sv_data %>%
   group_by(ChrStart, CentroChr) %>%
   mutate(hasCNChange=CentroCNChg > 0 & !is.na(CentroCNChg)) %>%
   summarise(nCNChange=sum(hasCNChange), n=n()) %>%
   mutate(cnChangeRate=nCNChange/n)


ggplot(chrheatmapdf) +
   aes(y=factor(CentroChr, levels=c(1:22, "X", "Y")), x=factor(ChrStart, levels=c(1:22, "X", "Y")), fill=n, label=n) +
   geom_tile() +
   scale_fill_gradient(low="#FFFFFF", high="#000000") +
   labs(y="Most likely centromere", x="Single breakend chromosome", fill="events")
figsave("fig3c_centromeric_single_breakend_heatmap", width=5, height=4)

ggplot(chrheatmapdf) +
   aes(y=factor(CentroChr, levels=c(1:22, "X", "Y")), x=factor(ChrStart, levels=c(1:22, "X", "Y")), fill=cnChangeRate, label=round(cnChangeRate, digits=2)) +
   geom_tile() +
   geom_text() + 
#   scale_fill_gradient(low="#FFFFFF", high="#000000") +
   labs(y="Most likely centromere", x="Single breakend chromosome", fill="CN change rate")


armlevels= paste0(rep(c(1:22, "X", "Y"), each=3), c("p", "?", "q"))
armconfheatmapdf = sgl_centro_sv_data %>%
   mutate(
      breakendArm=paste0(ChrStart, tolower(ArmStart)),
      centroArm=paste0(CentroChr, ifelse(PCentroGain, "p", ifelse(QCentroGain, "q", "?")))) %>%
   group_by(breakendArm, centroArm, HighConf) %>%
   mutate(hasCNChange=CentroCNChg > 0 & !is.na(CentroCNChg)) %>%
   summarise(nCNChange=sum(hasCNChange), n=n()) %>%
   mutate(cnChangeRate=nCNChange/n)
armheatmapdf = sgl_centro_sv_data %>%
   mutate(
      breakendArm=paste0(ChrStart, tolower(ArmStart)),
      centroArm=paste0(CentroChr, ifelse(PCentroGain, "p", ifelse(QCentroGain, "q", "?")))) %>%
   group_by(breakendArm, centroArm) %>%
   mutate(hasCNChange=CentroCNChg > 0 & !is.na(CentroCNChg)) %>%
   summarise(nCNChange=sum(hasCNChange), n=n()) %>%
   mutate(cnChangeRate=nCNChange/n)

ggplot(armheatmapdf %>% filter(!str_detect(centroArm, "Y"), !str_detect(breakendArm, "Y"))) +
   aes(y=factor(centroArm, levels=armlevels), x=factor(breakendArm, levels=armlevels), fill=n, label=n, color=
          ifelse(str_detect(centroArm, stringr::fixed("?")), "No CN change",
            ifelse(str_detect(centroArm, stringr::fixed("p")), "P gain", "Q gain"))) +
   geom_tile() +
   geom_text() + 
   scale_fill_gradient(low="#FFFFFF", high="#000000") +
   labs(y="Most likely centromere", x="Single breakend location", fill="events", colour="Arm")

ggplot(armheatmapdf %>% filter(!str_detect(centroArm, stringr::fixed("?")))) +
   aes(y=factor(centroArm, levels=armlevels), x=factor(breakendArm, levels=armlevels), fill=n, label=n) +
   geom_tile() +
   geom_text() + 
   scale_fill_gradient(low="#FFFFFF", high="#000000") +
   labs(y="Most likely centromere", x="Single breakend chromosome arm", fill="events")
figsave("supp_figure_centromeric_single_heatmap", width=12, height=10)

grBuckets = with(sglOrigPositions, GRanges(seqnames=Chromosome, ranges=IRanges(start=1000000*PosBucket-999999, width=1000000), n=n))
grNBuckets = with(sglOrigPositions %>% group_by(Chromosome, PosBucket) %>% do({data.frame(hitno=seq_len(.$n))}),
                  GRanges(seqnames=Chromosome, ranges=IRanges(start=1000000*PosBucket-999999, width=1000000)))

pdf(paste0(figdir, "/raw_fig_centromeric_single_location_intra_inter.pdf"))
chromPlot(bands=hg19_gaps, 
   annot1=with(sgl_centro_sv_data %>% filter(isIntraChromosomal), GRanges(seqnames=ChrStart, ranges=IRanges(start=PosStart, width=1), strand=OrientStart)),
   annot2=with(sgl_centro_sv_data %>% filter(!isIntraChromosomal), GRanges(seqnames=ChrStart, ranges=IRanges(start=PosStart, width=1), strand=OrientStart)),
   chrSide=c(-1,1,1,1,1,1,1,1),
   figCols=24)
dev.off()

# Single breakends to chr1
chromPlot(bands=hg19_gaps, 
   annot1=with(sgl_centro_sv_data %>%
      mutate(
       breakendArm=paste0(ChrStart, tolower(ArmStart)),
       centroArm=paste0(CentroChr, ifelse(PCentroGain, "p", ifelse(QCentroGain, "q", "?")))) %>%
      filter(CentroChr==1, QCentroGain),
      GRanges(seqnames=ChrStart, ranges=IRanges(start=PosStart, width=1), strand=OrientStart)),
   annot2=with(sgl_centro_sv_data %>%
      mutate(
         breakendArm=paste0(ChrStart, tolower(ArmStart)),
         centroArm=paste0(CentroChr, ifelse(PCentroGain, "p", ifelse(QCentroGain, "q", "?")))) %>%
      filter(CentroChr==1, PCentroGain),
      GRanges(seqnames=ChrStart, ranges=IRanges(start=PosStart, width=1), strand=OrientStart)),
   chrSide=c(-1,1,1,1,1,1,1,1),
   figCols=24)
# Single breakends to chr8
chromPlot(bands=hg19_gaps, 
   annot1=with(sgl_centro_sv_data %>%
                mutate(
                   breakendArm=paste0(ChrStart, tolower(ArmStart)),
                   centroArm=paste0(CentroChr, ifelse(PCentroGain, "p", ifelse(QCentroGain, "q", "?")))) %>%
                filter(CentroChr==8, QCentroGain),
             GRanges(seqnames=ChrStart, ranges=IRanges(start=PosStart, width=1), strand=OrientStart)),
   annot2=with(sgl_centro_sv_data %>%
                mutate(
                   breakendArm=paste0(ChrStart, tolower(ArmStart)),
                   centroArm=paste0(CentroChr, ifelse(PCentroGain, "p", ifelse(QCentroGain, "q", "?")))) %>%
                filter(CentroChr==8, PCentroGain),
             GRanges(seqnames=ChrStart, ranges=IRanges(start=PosStart, width=1), strand=OrientStart)),
   chrSide=c(-1,1,1,1,1,1,1,1),
   figCols=24)
# chr1, chr8
for (current_chr in c(1:22, "X", "Y")) {
   pdf(paste0(figdir, "/details_centromeric_single_location_chr", current_chr, ".pdf"))
   chromPlot(bands=hg19_gaps, 
          annot2=with(sgl_centro_sv_data %>% filter(CentroChr==current_chr),
                      GRanges(seqnames=ChrStart, ranges=IRanges(start=PosStart, width=1), strand=OrientStart)),
          annot1=with(sgl_centro_sv_data,
                      GRanges(seqnames=ChrStart, ranges=IRanges(start=PosStart, width=1), strand=OrientStart)),
          chrSide=c(-1,1,1,1,1,1,1,1),
          figCols=24)
   dev.off()
}

sgl_centro_sv_data %>%
   group_by(SampleId, CentroChr) %>%
   summarise(singles=n()) %>%
   group_by(singles, CentroChr) %>%
   summarise(count=n()) %>%
ggplot() +
   aes(x=factor(CentroChr, levels=c(1:22, "X", "Y")), y=count, fill=factor(pmin(singles, 4), levels=4:0)) +
   geom_bar(position="fill",stat='identity') +
   labs(title="Single breakends to arm", x="Most likely centromere", fill="Number of single breakends")
figsave("supp_figure_number_of_centromeric_singles_per_chr")

sgl_centro_sv_data %>%
   group_by(SampleId, ChrStart) %>%
   summarise(singles=n()) %>%
   group_by(singles, ChrStart) %>%
   summarise(count=n()) %>%
   ggplot() +
   aes(x=factor(ChrStart, levels=c(1:22, "X", "Y")), y=count, fill=factor(pmin(singles, 4), levels=4:0)) +
   geom_bar(position="fill",stat='identity') +
   labs(title="Single breakends from arm")













