source("libbenchmark.R")
library(openxlsx)
library(ggExtra)

#gsutil -m cp gs://patient-report-bucket-umc/**/*.purple.qc .
#gsutil -m cp gs://patient-report-bucket-umc/**/*.purple.sv.vcf.gz .
#gsutil -m cp gs://patient-report-bucket-umc/**/*.purple.cnv.somatic.tsv .
#grep QCStatus *.qc | sed -e 's/T.purple.qc:QCStatus//' > qcstatus.txt
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
consensus_pcawg_sv_cn_transitions_list = lapply(metadata_df$SAMPLE, pcawg_evaluate_cn_transitions)
consensus_pcawg_cn_transitions_gr = unlist(GRangesList(lapply(consensus_pcawg_sv_cn_transitions_list, function(x) { x$cn_transition } )))
consensus_pcawg_sv_transitions_gr = unlist(GRangesList(lapply(consensus_pcawg_sv_cn_transitions_list, function(x) { x$sv } )))

save.image(paste0(privatedatadir, "pcawg1528_svcn_consistency.RData"))
#load(paste0(privatedatadir, "pcawg1528_svcn_consistency.RData"))


cn_transitions = bind_rows(
	consensus_pcawg_cn_transitions_gr %>% as.data.frame() %>% mutate(cohort="Consensus"),
	gpl_pcawg_cn_transitions_gr %>% as.data.frame() %>% mutate(cohort="GRIDSS2")) %>%
	mutate(
		inGap = overlapsAny(GRanges(seqnames=seqnames, ranges=IRanges(start=start, end=end)), hg19_gaps, maxgap=100000),
		inCentromere = overlapsAny(GRanges(seqnames=seqnames, ranges=IRanges(start=start, end=end)), hg19_centromeres, maxgap=100000),
		nearGapOrCentromere=inGap | inCentromere) %>%
	inner_join(metadata_df, by=c("sampleId"="SAMPLE"))

sv_transitions = bind_rows(
	consensus_pcawg_sv_transitions_gr %>% as.data.frame() %>% mutate(cohort="Consensus"),
	gpl_pcawg_sv_transitions_gr %>% as.data.frame() %>% mutate(cohort="GRIDSS2")) %>%
	mutate(
		inGap = overlapsAny(GRanges(seqnames=seqnames, ranges=IRanges(start=start, end=end)), hg19_gaps, maxgap=100000),
		inCentromere = overlapsAny(GRanges(seqnames=seqnames, ranges=IRanges(start=start, end=end)), hg19_centromeres, maxgap=100000),
		nearGapOrCentromere=inGap | inCentromere,
		distance=ifelse(is.na(distance), 100000, distance)) %>%
	inner_join(metadata_df, by=c("sampleId"="SAMPLE"))

hartwig_bin_width = 1
cnsvmatchdf = cn_transitions %>%
	filter(
		!nearGapOrCentromere,
		!is.na(cn_left),
		!is.na(cn_right)) %>%
	mutate(
		classification=factor(ifelse(sv_match_rescued, "Rescued", sv_match_classification), levels=rev(c("Missing SV", "Rescued", "Single Breakend", "Breakpoint"))))
cnsvmatchdf = bind_rows(
	cnsvmatchdf %>%
		filter(cohort=="Consensus") %>%
		mutate(
			cnbin=round(pmin(abs(cn_delta), 8)),
			xmin=cnbin-0.5,
			xmax=cnbin+0.5,
			xwidth=xmax-xmin),
	cnsvmatchdf %>%
		filter(cohort=="GRIDSS2") %>%
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
	ungroup()

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
figsave("pcawg_cn_transistions_missing_breakpoints_by_cndelta_percentage", width=5, height=4)

ggplot(cn_transistions_missing_breakpoints_by_cndelta_plot_df) +
	aes(x=cnbin, fill=classification, weight=n, width=xwidth) +
	geom_bar(position="stack") + 
	scale_fill_manual(values=c("#DDEEDD", "#1b9e77", "#fb6a4a")) +
	scale_x_continuous(expand=c(0,0), breaks = 0:8, labels = 0:8) +
	coord_cartesian(xlim=c(0, 8)) +
	facet_wrap(~cohort) +
	labs(x="Magnitude of copy number change", y="Copy number transitions") +
	theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
figsave("pcawg_cn_transistions_missing_breakpoints_by_cndelta_counts", width=5, height=4)

# Figures in paper
cnsvmatchdf %>% group_by(cohort, classification) %>%
	summarise(count=n()) %>%
	group_by(cohort) %>%
	mutate(
		total=sum(count),
		pct=100*count/total)




test_bin_width = 0.25
good_cn_transitions = cn_transitions %>%
	filter(!nearGapOrCentromere & !is.na(cn_delta)) %>%
	mutate(
		classification=factor(ifelse(sv_match_rescued, "Rescued", sv_match_classification), levels=rev(c("Missing SV", "Rescued", "Single Breakend", "Breakpoint"))),
		isfn=sv_matches == 0,
		cnbin=round(pmin(abs(cn_delta), 8)/test_bin_width)*test_bin_width) %>%
	inner_join(metadata_df, by=c("sampleId"="SAMPLE"))

ggplot(good_cn_transitions) + 
	aes(x=cnbin, fill=classification) +
	geom_histogram() +
	facet_grid(QCStatus ~ cohort, scales="free")

ggplot(good_cn_transitions) + 
	aes(x=cnbin, fill=classification) +
	geom_histogram() +
	facet_grid(Project_Code_Count ~ cohort, scales="free")

ggplot(good_cn_transitions) + 
	aes(x=cnbin, fill=classification) +
	geom_histogram() +
	facet_grid(sampleId ~ cohort, scales="free")

