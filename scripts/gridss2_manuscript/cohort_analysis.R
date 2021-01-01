source("libbenchmark.R")
library(openxlsx)
library(ggExtra)

#### Hartwig ####
lnx_svs = read_csv(paste0(privatedatadir, "/hartwig/LNX_SVS.csv"))
lnx_cns = read_tsv(paste0(privatedatadir, "/hartwig/LNX_VIS_COPY_NUMBER.tsv"))
hartwig_svgr = lnx_svs %>%
	filter(!(ChrEnd == "0" & is.na(InsertSeq))) %>% # strip purple placeholder calls
#	filter(!Recovered) %>%
	lnx_to_gr()
hartwig_svgr$sampleId = hartwig_svgr$SampleId
hartwig_cngr = with(lnx_cns, GRanges(
	seqnames=Chromosome,
	ranges=IRanges(start=Start, end=End),
	cn=CopyNumber,
	cn_major=CopyNumber * BAF,
	cn_minor=CopyNumber * (1 - BAF),
	sampleId=SampleId))
hartwig_evaluate_cn_transitions = function(sampleId, ...) {
	write(paste("Processing ", sampleId), stderr())
	svgr = hartwig_svgr[hartwig_svgr$sampleId == sampleId,]
	cngr = hartwig_cngr[hartwig_cngr$sampleId == sampleId,]
	cn_consistency = evaluate_cn_transitions(cngr, svgr, ...)
	return(cn_consistency)
}
hartwig_sv_cn_transitions_list = lapply(unique(lnx_svs$SampleId), hartwig_evaluate_cn_transitions)
gpl_hartwig_cn_transitions_gr = unlist(GRangesList(lapply(hartwig_sv_cn_transitions_list, function(x) { x$cn_transition } )))
gpl_hartwig_sv_transitions_gr = unlist(GRangesList(lapply(hartwig_sv_cn_transitions_list, function(x) { x$sv } )))

#### PCAWG Consensus ####
purple_qc_status = read_tsv(paste0(privatedatadir, "/pcawg/qcstatus.txt"), col_types="cc", col_names=c("donor_ID", "QCStatus"))

fullmddf = read.xlsx(paste0(privatedatadir, "/pcawg/PCAWG_meta.xlsx"))
#names(fullmddf) = make.names(names(fullmddf), unique=TRUE)
#raw_metadata_df = read.xlsx(paste0(privatedatadir, "/pcawg/PCAWG_meta.xlsx"), sheet="GPL")
metadata_df = fullmddf %>%
	inner_join(purple_qc_status, by="donor_ID") %>%
	filter(
		!(SAMPLE %in% SAMPLE[duplicated(SAMPLE)]),
		!(donor_ID %in% donor_ID[duplicated(donor_ID)]),
		file.exists(paste0(pcawg_dir, "/", SAMPLE, ".consensus.20170119.somatic.cna.txt")),
		file.exists(paste0(pcawg_dir, "/", SAMPLE, ".pcawg_consensus_1.6.161116.somatic.sv.bedpe")),
		file.exists(paste0(privatedatadir, "/pcawg/", donor_ID, "T.purple.sv.vcf.gz")),
		file.exists(paste0(privatedatadir, "/pcawg/", donor_ID, "T.purple.cnv.somatic.tsv")),
	)
write(paste("Found files and unique donor_ID<->SAMPLE mapping for", nrow(metadata_df), "of", "PCAWG samples."), stderr())
metadata_df = metadata_df %>% filter(QCStatus == "PASS")
write(paste(nrow(metadata_df), "passing PURPLE QC"), stderr())

consensus_pcawg_sv_cn_transitions_list = lapply(metadata_df$SAMPLE, pcawg_evaluate_cn_transitions)
consensus_pcawg_cn_transitions_gr = unlist(GRangesList(lapply(consensus_pcawg_sv_cn_transitions_list, function(x) { x$cn_transition } )))
consensus_pcawg_sv_transitions_gr = unlist(GRangesList(lapply(consensus_pcawg_sv_cn_transitions_list, function(x) { x$sv } )))

#### Hartwig Pipeline on PCAWG data ####
pcawg_gpl_evaluate_cn_transitions = function(sampleId) {
	donorId = metadata_df$donor_ID[metadata_df$SAMPLE==sampleId]
	write(paste("Processing ", sampleId, donorId), stderr())
	gridss_vcf = readVcf(paste0(privatedatadir, "/pcawg/", donorId, "T.purple.sv.vcf.gz"))
	purple_df = read_table2(paste0(privatedatadir, "/pcawg/", donorId, "T.purple.cnv.somatic.tsv"), col_types="ciiddddcccidiidd")
	gridss_vcf = gridss_vcf[rowRanges(gridss_vcf)$FILTER == "PASS",]
	purple_gr = with(purple_df, GRanges(
		seqnames=chromosome,
		ranges=IRanges(start=start, end=end),
		cn=minorAllelePloidy + majorAllelePloidy,
		cn_major=minorAllelePloidy,
		cn_minor=majorAllelePloidy,
		sampleId=sampleId))
	gridss_gr = c(
		breakpointRanges(gridss_vcf, nominalPosition=TRUE),
		breakendRanges(gridss_vcf, nominalPosition=TRUE))
	if (length(gridss_gr) > 0) {
		names(gridss_gr) = paste0(sampleId, "_", names(gridss_gr))
		gridss_gr$partner = ifelse(is.na(gridss_gr$partner), NA, paste0(sampleId, "_", gridss_gr$partner))
	}
	result = evaluate_cn_transitions(purple_gr, gridss_gr)
	result$sv$sampleId=rep(sampleId, length(result$sv))
	result$cn_transition$sampleId=rep(sampleId, length(result$cn_transition))
	return(result)
}

gpl_pcawg_sv_cn_transitions_list = lapply(metadata_df$SAMPLE, pcawg_gpl_evaluate_cn_transitions)
gpl_pcawg_cn_transitions_gr = unlist(GRangesList(lapply(gpl_pcawg_sv_cn_transitions_list, function(x) { x$cn_transition } )))
gpl_pcawg_sv_transitions_gr = unlist(GRangesList(lapply(gpl_pcawg_sv_cn_transitions_list, function(x) { x$sv } )))

#### Snapshot ####

#save.image(paste0(privatedatadir, "cohort_svcn_consistency.RData"))
#load(paste0(privatedatadir, "cohort_svcn_consistency.RData"))

#### Analysis ####

cn_transitions = bind_rows(
	gpl_hartwig_cn_transitions_gr %>% as.data.frame(row.names=NULL) %>% mutate(cohort="Hartwig GRIDSS2"),
	gpl_pcawg_cn_transitions_gr %>% as.data.frame(row.names=NULL) %>% mutate(cohort="PCAWG GRIDSS2"),
	consensus_pcawg_cn_transitions_gr %>% as.data.frame(row.names=NULL) %>% mutate(cohort="PCAWG Consensus")) %>%
	mutate(
		inGap = overlapsAny(GRanges(seqnames=seqnames, ranges=IRanges(start=start, end=end)), hg19_gaps, maxgap=100000),
		inCentromere = overlapsAny(GRanges(seqnames=seqnames, ranges=IRanges(start=start, end=end)), hg19_centromeres, maxgap=100000),
		nearGapOrCentromere=inGap | inCentromere)

# evaluate_cn_transitions doesn't annotate breakend cn so we can instead load from the files
lnx_svs %>%
	as.data.frame(row.names=NULL) %>%
	mutate(
		sv_match_classification=ifelse(ChrEnd==0, "Single Breakend", "Breakpoint"),
		cohort="Hartwig GRIDSS2") %>%
	dplyr::select(seqnames=ChrStart, start=PosStart, end=PosStart, strand=OrientStart, beid=Id, sampleId=SampleId, cn_left=CNStart, cn_delta=CNChgStart, ploidy=Ploidy, sv_match_classification, cohort)

sv_transitions = bind_rows(
	gpl_hartwig_sv_transitions_gr %>% as.data.frame(row.names=NULL) %>% mutate(cohort="Hartwig GRIDSS2"),
	gpl_pcawg_sv_transitions_gr %>% as.data.frame(row.names=NULL) %>% mutate(cohort="PCAWG GRIDSS2"),
	consensus_pcawg_sv_transitions_gr %>%
		as.data.frame() %>%
		mutate(
			cohort="PCAWG Consensus",
			cn_delta=abs(cn_left - cn_right))) %>%
	mutate(
		inGap = overlapsAny(GRanges(seqnames=seqnames, ranges=IRanges(start=start, end=end)), hg19_gaps, maxgap=100000),
		inCentromere = overlapsAny(GRanges(seqnames=seqnames, ranges=IRanges(start=start, end=end)), hg19_centromeres, maxgap=100000),
		nearGapOrCentromere=inGap | inCentromere,
		distance=ifelse(is.na(distance), 100000, distance))

## Exploratory plots

## Figure 2
hartwig_bin_width = 0.25
cnsvmatchdf = cn_transitions %>%
	filter(!nearGapOrCentromere) %>%
	mutate(
		classification=factor(ifelse(sv_match_classification=="Rescued", "Missing SV", sv_match_classification), levels=rev(c("Missing SV", "Rescued", "Single Breakend", "Breakpoint"))))
		#classification=factor(ifelse(sv_match_rescued, "Rescued", sv_match_classification), levels=rev(c("Missing SV", "Rescued", "Single Breakend", "Breakpoint"))))
cnsvmatchdf = bind_rows(
	cnsvmatchdf %>%
		filter(cohort=="PCAWG Consensus") %>%
		mutate(
			cnbin=round(pmin(abs(cn_delta), 8)),
			xmin=cnbin-0.5,
			xmax=cnbin+0.5,
			xwidth=xmax-xmin),
	cnsvmatchdf %>%
		filter(cohort!="PCAWG Consensus") %>%
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
	mutate(cohort=factor(cohort, levels=c("PCAWG Consensus", "PCAWG GRIDSS2", "Hartwig GRIDSS2")))

ggplot(cn_transistions_missing_breakpoints_by_cndelta_plot_df) +
	aes(x=cnbin, fill=classification, y=percentage, width=xwidth) +
	scale_x_continuous(expand=c(0,0), breaks = 0:8, labels = 0:8) +
	scale_y_continuous(expand=c(0,0, 0, 0.01), breaks=seq(0, 0.3, 0.05), labels=paste0(seq(0, 0.3, 0.05)*100, "%")) +
	coord_cartesian(ylim=c(0, 0.31), xlim=c(0, 8)) +
	geom_col(position="stack") +
	scale_fill_manual(values=c("#DDEEDD", "#1b9e77", "#fb6a4a")) +
	facet_wrap(~cohort) +
	labs(x="Magnitude of copy number change", y="Portion of copy number transitions") +
	theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
figsave("cn_transistions_missing_breakpoints_by_cndelta_percentage", width=8, height=4)

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

svcnmatchdf = sv_transitions %>%
	group_by(cohort, sampleId, seqnames) %>%
	arrange(start) %>%
	mutate(
		isolated_l = is.na(lead(start)) | abs(start - lead(start)) > isolation_distance,
		isolated_r = is.na( lag(start)) | abs(start -  lag(start)) > isolation_distance,
		isolated=isolated_l & isolated_r) %>%
	filter(!nearGapOrCentromere, !is.na(cn_left) & !is.na(cn_right), isolated) %>%
	mutate(
		classification=factor(ifelse(is.na(partner), "Single Breakend", "Breakpoint"), levels=rev(c("Missing SV", "Rescued", "Single Breakend", "Breakpoint"))),
		cn_delta=abs(cn_left - cn_right))
svcnmatchdf = bind_rows(
	svcnmatchdf %>%
		filter(cohort=="PCAWG Consensus") %>%
		mutate(
			cnbin=round(pmin(abs(cn_delta), 8)),
			xmin=cnbin-0.5,
			xmax=cnbin+0.5,
			xwidth=xmax-xmin),
	svcnmatchdf %>%
		filter(cohort!="PCAWG Consensus") %>%
		mutate(
			cnbin=round(pmin(abs(cn_delta), 8)/hartwig_bin_width)*hartwig_bin_width,
			xmin=cnbin-hartwig_bin_width/2,
			xmax=cnbin+hartwig_bin_width/2,
			xwidth=xmax-xmin))
sv_transistions_by_cndelta_plot_df = svcnmatchdf %>%
	group_by(cohort, cnbin, xmin, xmax, xwidth, classification) %>%
	summarise(
		count = n()) %>%
	group_by(cohort) %>%
	mutate(pct = count/sum(count))

ggplot(sv_transistions_by_cndelta_plot_df) +
	aes(x=cnbin, fill=classification, y=pct, width=xwidth) +
	geom_col(position="stack") +
	facet_wrap(~cohort) +
	scale_x_continuous(expand=c(0,0), breaks = 0:8, labels = 0:8) +
	scale_y_continuous(expand=c(0,0, 0, 0.01)) +
	facet_wrap(~cohort) +
	labs(x="Magnitude of copy number change", y="Portion of structural variants transitions") +
	theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
#figsave("sv_transistions_by_cndelta_percentage", width=5, height=4)


### Numbers in paper
cnsvmatchdf %>%
	mutate(classification=factor(ifelse(sv_match_rescued, "Missing SV", sv_match_classification), levels=rev(c("Missing SV", "Rescued", "Single Breakend", "Breakpoint")))) %>%
	group_by(cohort, classification) %>%
	summarise(count=n()) %>%
	group_by(cohort) %>%
	mutate(pct=100*count/sum(count), total=sum(count))

cnsvmatchdf %>%
	filter(cn_delta > 0.75) %>% # edge of CN=1 cluster
	mutate(classification=factor(ifelse(sv_match_rescued, "Missing SV", sv_match_classification), levels=rev(c("Missing SV", "Rescued", "Single Breakend", "Breakpoint")))) %>%
	group_by(cohort, classification) %>%
	summarise(count=n()) %>%
	group_by(cohort) %>%
	mutate(pct=100*count/sum(count))

cnsvmatchdf %>%
	filter(cn_delta <= 0.75) %>% # edge of CN=1 cluster
	mutate(classification=factor(ifelse(sv_match_rescued, "Missing SV", sv_match_classification), levels=rev(c("Missing SV", "Rescued", "Single Breakend", "Breakpoint")))) %>%
	group_by(cohort, classification) %>%
	summarise(count=n()) %>%
	group_by(cohort) %>%
	mutate(pct=100*count/sum(count))


