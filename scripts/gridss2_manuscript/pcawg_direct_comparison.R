setwd("../")
source("libgridss.R")
setwd("gridss2_manuscript")
source("libbenchmark.R")
library(openxlsx)
library(ggExtra)

#gsutil -m cp gs://pcawg/**/*.purple.qc .
#gsutil -m cp gs://pcawg/**/*.purple.sv.vcf.gz .
#gsutil -m cp gs://pcawg/**/*.purple.cnv.somatic.tsv .
#grep QCStatus *.qc | sed -e 's/T.purple.qc:QCStatus//' > qcstatus.txt
purple_qc_status = read_tsv(paste0(privatedatadir, "/pcawg/qcstatus.txt"), col_types="cc", col_names=c("donor_ID", "QCStatus"))

fullmddf = read.xlsx(paste0(privatedatadir, "/pcawg/PCAWG_meta.xlsx"))
#names(fullmddf) = make.names(names(fullmddf), unique=TRUE)
#raw_metadata_df = read.xlsx(paste0(privatedatadir, "/pcawg/PCAWG_meta.xlsx"), sheet="GPL")
metadata_df = fullmddf %>%
	inner_join(purple_qc_status, by="donor_ID") %>%
	mutate(
		purple_qc_pass=QCStatus == "PASS",
		purple_qc_ok=
			!str_detect(QCStatus, "FAIL_NO_TUMOR") &
			!str_detect(QCStatus, "WARN_HIGH_COPY_NUMBER_NOISE")& 
			!str_detect(QCStatus, "WARN_GENDER_MISMATCH")&
			!str_detect(QCStatus, "WARN_DELETED_GENES"),
		has_sv_consensus=file.exists(paste0(pcawg_dir, "/", SAMPLE, ".pcawg_consensus_1.6.161116.somatic.sv.bedpe")),
		has_cnv_consensus=file.exists(paste0(pcawg_dir, "/", SAMPLE, ".consensus.20170119.somatic.cna.txt")),
		has_purple_sv=file.exists(paste0(privatedatadir, "/pcawg/", donor_ID, "T.purple.sv.vcf.gz")),
		has_purple_cnv=file.exists(paste0(privatedatadir, "/pcawg/", donor_ID, "T.purple.cnv.somatic.tsv")),
		sample_multi=SAMPLE %in% SAMPLE[duplicated(SAMPLE)],
		donor_multi=donor_ID %in% donor_ID[duplicated(donor_ID)])
#metadata_df %>%
#	group_by(purple_qc_pass, purple_qc_ok, has_sv_consensus, has_cnv_consensus, has_purple_sv, has_purple_cnv, sample_multi, donor_multi) %>%
#	summarise(n=n()) %>%
#	View()
write(paste(nrow(purple_qc_status), "PURPLE results"), stderr())
write(paste(nrow(metadata_df), "with PURPLE results"), stderr())
metadata_df = metadata_df %>% filter(!donor_multi)
write(paste(nrow(metadata_df), "after removing donors with multiple specimens sequenced"), stderr())
metadata_df = metadata_df %>% filter(has_sv_consensus & has_cnv_consensus)
write(paste(nrow(metadata_df), "after removing missing PCAWG consensus"), stderr())
metadata_df = metadata_df %>% filter(has_purple_sv & has_purple_cnv)
write(paste(nrow(metadata_df), "after removing missing PURPLE results"), stderr())
metadata_df = metadata_df %>% filter(purple_qc_ok)
write(paste(nrow(metadata_df), "after removing PURPLE QC not ok"), stderr())
metadata_df = metadata_df %>% filter(purple_qc_pass)
write(paste(nrow(metadata_df), "after removing PURPLE QC not PASS"), stderr())

pcawg_compare_sample_run = function(sampleId, svgr, svcaller, cngr, cncaller, ...) {
	caller = paste(svcaller, cncaller)
	result = evaluate_cn_transitions(cngr, svgr, ...)
	result$sv$sampleId=rep(sampleId, length(result$sv))
	result$sv$caller=rep(caller, length(result$sv))
	result$cn_transitions$sampleId=rep(sampleId, length(result$cn_transitions))
	result$cn_transitions$caller=rep(caller, length(result$cn_transitions))
	return(result)
}
pcawg_compare_sample = function(sampleId) {
	donorId = metadata_df$donor_ID[metadata_df$SAMPLE==sampleId]
	write(paste("Processing ", sampleId, donorId), stderr())
	gridss_gr = load_gr_sv_gridss(sampleId, donorId)
	purple_gr = load_gr_cn_purple(sampleId, donorId)
	pcawgsv_gr = load_gr_sv_pcawg(sampleId, donorId)
	pcawgcn_gr = load_gr_cn_pcawg(sampleId, donorId)
	results = list(
		pcawg_compare_sample_run(sampleId, gridss_gr, "gridss", purple_gr, "purple"),
		pcawg_compare_sample_run(sampleId, gridss_gr, "gridss", pcawgcn_gr, "pcawg"),
		pcawg_compare_sample_run(sampleId, pcawgsv_gr, "pcawg", purple_gr, "purple"),
		pcawg_compare_sample_run(sampleId, pcawgsv_gr, "pcawg", pcawgcn_gr, "pcawg"))
	result = list(
		sv = unlist(GRangesList(lapply(results, function(x) x$sv))),
		cn_transitions = unlist(GRangesList(lapply(results, function(x) x$cn_transitions))))
	return(result)
}
icgc_sv_cn_transitions_list = lapply(metadata_df$SAMPLE, pcawg_compare_sample)
save.image(paste0(privatedatadir, "pcawg3_100purple3_", nrow(metadata_df), "_svcn_consistency.RData"))
#load(paste0(privatedatadir, "pcawg3_100purple3_", nrow(metadata_df), "_svcn_consistency.RData"))


cn_transitions = unlist(GRangesList(lapply(icgc_sv_cn_transitions_list, function(x) x$cn_transitions))) %>%
	as.data.frame() %>%
	mutate(
		inGap = overlapsAny(GRanges(seqnames=seqnames, ranges=IRanges(start=start, end=end)), hg19_gaps, maxgap=100000),
		inCentromere = overlapsAny(GRanges(seqnames=seqnames, ranges=IRanges(start=start, end=end)), hg19_centromeres, maxgap=100000),
		nearGapOrCentromere=inGap | inCentromere) %>%
	inner_join(metadata_df, by=c("sampleId"="SAMPLE"))
sv_transitions = unlist(GRangesList(lapply(icgc_sv_cn_transitions_list, function(x) x$sv))) %>%
	as.data.frame(row.names=NULL) %>%
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
		classification=factor(ifelse(sv_match_rescued, "Rescued", sv_match_classification), levels=rev(c("Missing SV", "Rescued", "Single Breakend", "Breakpoint"))),
		classification_100k=factor(ifelse(sv_match_rescued, "Rescued", sv_match_classification_100k), levels=rev(c("Missing SV", "Rescued", "Single Breakend", "Breakpoint"))))
cnsvmatchdf = bind_rows(
	cnsvmatchdf %>%
		filter(!str_detect(caller, "purple")) %>%
		mutate(
			cnbin=round(pmin(abs(cn_delta), 8)),
			xmin=cnbin-0.5,
			xmax=cnbin+0.5,
			xwidth=xmax-xmin),
	cnsvmatchdf %>%
	filter(str_detect(caller, "purple")) %>%
		mutate(
			cnbin=round(pmin(abs(cn_delta), 8)/hartwig_bin_width)*hartwig_bin_width,
			xmin=cnbin-hartwig_bin_width/2,
			xmax=cnbin+hartwig_bin_width/2,
			xwidth=xmax-xmin)) %>%
	mutate(cohort=caller)

cn_transistions_missing_breakpoints_by_cndelta_plot_df = cnsvmatchdf %>%
	group_by(cohort, cnbin, xmin, xmax, xwidth, classification) %>%
	summarise(n=n()) %>%
	group_by(cohort, cnbin, xmin, xmax, xwidth) %>%
	mutate(
		total_n = sum(n),
		percentage=n/total_n) %>%
	ungroup() %>%
	mutate(cohort=factor(cohort, levels=c("pcawg pcawg", "gridss purple", "pcawg purple", "gridss pcawg")))

ggplot(cn_transistions_missing_breakpoints_by_cndelta_plot_df %>% filter(cohort %in% c("pcawg pcawg", "gridss purple"))) +
	aes(x=cnbin, fill=classification, y=percentage, width=xwidth) +
	scale_x_continuous(expand=c(0,0), breaks = 0:8, labels = 0:8) +
	scale_y_continuous(expand=c(0,0, 0, 0.01), breaks=seq(0, 0.3, 0.05), labels=paste0(seq(0, 0.3, 0.05)*100, "%")) +
	coord_cartesian(ylim=c(0, 0.30), xlim=c(0, 8)) +
	geom_col(position="stack") +
	scale_fill_manual(values=c("#DDEEDD", "#1b9e77", "#fb6a4a")) +
	facet_wrap(~cohort, scale="free_y") +
	labs(x="Magnitude of copy number change", y="Portion of copy number transitions") +
	theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
figsave("pcawg_cn_transistions_missing_breakpoints_by_cndelta_percentage", width=5, height=4)

ggplot(cn_transistions_missing_breakpoints_by_cndelta_plot_df %>% filter(cohort %in% c("pcawg purple", "gridss pcawg"))) +
	aes(x=cnbin, fill=classification, y=percentage, width=xwidth) +
	scale_x_continuous(expand=c(0,0), breaks = 0:8, labels = 0:8) +
	scale_y_continuous(expand=c(0,0, 0, 0.01), breaks=seq(0, 0.6, 0.05), labels=paste0(seq(0, 0.6, 0.05)*100, "%")) +
	coord_cartesian(ylim=c(0, 0.60), xlim=c(0, 8)) +
	geom_col(position="stack") +
	scale_fill_manual(values=c("#DDEEDD", "#1b9e77", "#fb6a4a")) +
	facet_wrap(~cohort, scale="free_y") +
	labs(x="Magnitude of copy number change", y="Portion of copy number transitions") +
	theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
figsave("cross_pcawg_cn_transistions_missing_breakpoints_by_cndelta_percentage", width=5, height=4)


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
# Clonal transitions
cnsvmatchdf %>% group_by(cohort, classification) %>%
	filter(abs(cn_delta) > 0.5) %>%
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
	facet_grid(QCStatus.x ~ caller, scales="free")

ggplot(good_cn_transitions) + 
	aes(x=cnbin, fill=classification) +
	geom_histogram() +
	facet_grid(Project_Code_Count ~ cohort, scales="free")

ggplot(good_cn_transitions) + 
	aes(x=cnbin, fill=classification) +
	geom_histogram() +
	facet_grid(sampleId ~ cohort, scales="free")


# Analysis: why is the cross-correlation between the calls sets so bad?
# Is it the clustered calling?
# let's inspect the first sample
z = cnsvmatchdf %>% filter(sampleId == "f393ba16-9361-5df4-e040-11ac0d4844e8") %>% arrange(seqnames, start) %>% filter(cohort %in% c("gridss pcawg", "pcawg pcawg"))
z %>% dplyr::select(sampleId, donor_ID, seqnames, start, end, cn_left, cn_right, cn_delta, nearGapOrCentromere) %>%
	distinct() %>%
	left_join(z %>% dplyr::select(sampleId, seqnames, start, cohort, classification, distance) %>% filter(cohort=="gridss pcawg"), by=c("sampleId", "seqnames", "start")) %>%
	left_join(z %>% dplyr::select(sampleId, seqnames, start, cohort, classification, distance) %>% filter(cohort=="pcawg pcawg"), by=c("sampleId", "seqnames", "start"), suffix=c("",".pcawg")) %>%
	filter(classification != classification.pcawg) %>%
	View()
# Ok, so evaluate_cn_transitions doesn't like inexact matches. Let's replot with distance:
cnsvmatchdf %>%
	group_by(cohort, cnbin, xmin, xmax, xwidth, classification_100k) %>%
	summarise(n=n()) %>%
	group_by(cohort, cnbin, xmin, xmax, xwidth) %>%
	mutate(
		total_n = sum(n),
		percentage=n/total_n) %>%
	ungroup() %>%
ggplot() +
	aes(x=cnbin, fill=classification_100k, y=percentage, width=xwidth) +
	scale_x_continuous(expand=c(0,0), breaks = 0:8, labels = 0:8) +
	scale_y_continuous(expand=c(0,0, 0, 0.01), breaks=seq(0, 0.3, 0.05), labels=paste0(seq(0, 0.3, 0.05)*100, "%")) +
	coord_cartesian(ylim=c(0, 1), xlim=c(0, 8)) +
	geom_col(position="stack") +
	scale_fill_manual(values=c("#DDEEDD", "#1b9e77", "#fb6a4a")) +
	facet_wrap(~cohort) +
	labs(x="Magnitude of copy number change", y="Portion of copy number transitions") +
	theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
cnsvmatchdf %>%
	mutate(classification=ifelse(is.na(distance), "3Missing", ifelse(distance==1, "1Exact", "2Inexact"))) %>%
	group_by(cohort, cnbin, xmin, xmax, xwidth, classification) %>%
	summarise(n=n()) %>%
	group_by(cohort, cnbin, xmin, xmax, xwidth) %>%
	mutate(
		total_n = sum(n),
		percentage=n/total_n) %>%
ggplot() +
	aes(x=cnbin, fill=classification, weight=n, width=xwidth) +
	geom_bar(position="stack") + 
	scale_fill_manual(values=c("#DDEEDD", "#1b9e77", "#fb6a4a")) +
	scale_x_continuous(expand=c(0,0), breaks = 0:8, labels = 0:8) +
	coord_cartesian(xlim=c(0, 8)) +
	facet_wrap(~cohort) +
	labs(x="Magnitude of copy number change", y="Copy number transitions") +
	theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
cnsvmatchdf %>% group_by(caller) %>%
	summarise(
		total=n(),
		missing=sum(is.na(distance)),
		fnr=missing/total)
# TODO
# - fix imprecise matches
#  - greedy assign sv -> cn based on min distance
# - output flanking CN segment lengths
# - output distance to nearest SV?

crossmatchdf = cnsvmatchdf %>%
	mutate(cn_caller = str_extract(caller, " .*$"))
crossmatchdf = inner_join(crossmatchdf, crossmatchdf, by=c("sampleId", "seqnames", "start", "cn_caller"), suffix=c(".x", ".y")) %>%
	filter(
		((caller.x == "gridss purple" & caller.y == "pcawg purple") |
		 	(caller.x == "pcawg pcawg" & caller.y == "gridss pcawg")) &
			classification_100k.x == "Missing SV") 
crossmatchdf %>%
	group_by(caller.x, caller.y, classification_100k.x, classification_100k.y) %>%
	summarise(
		events=n(),
		sv_partner_matches.y=sum(sv_partner_matches.y)) %>%
	group_by(caller.x, caller.y) %>%
	mutate(pct=100*events/sum(events))
	
#crossmatchdf %>% filter(caller.x == "gridss purple") %>% View()
crossmatchdf %>%
	filter(classification_100k.x=="Missing SV" & classification_100k.y != classification_100k.x) %>%
	group_by(sampleId, donor_ID.x, caller.x, caller.y, classification_100k.x, classification_100k.y) %>%
	summarise(n=n()) %>%
ggplot() +
	aes(x=n) +
	geom_histogram(binwidth=1) +
	facet_wrap(caller.x ~ caller.y)

crossmatchdf %>%
	filter(classification_100k.x=="Missing SV" & classification_100k.y != classification_100k.x) %>%
	group_by(sampleId, donor_ID.x, caller.x, caller.y, classification_100k.x, classification_100k.y) %>%
	summarise(total=n()) %>%
	group_by(caller.x) %>%
	arrange(caller.x, desc(total)) %>%
	mutate(
		rn=row_number(),
		csum=cumsum(total)) %>%
	View()
cnsvmatchdf %>% filter(classification_100k=="Missing SV" & caller == "gridss purple" & donor_ID=="DO52686") %>% View()

crossmatchdf = cnsvmatchdf %>% filter(cohort %in% c("gridss pcawg", "pcawg pcawg"))
crossmatchdf %>% dplyr::select(sampleId, donor_ID, seqnames, start, end, cn_left, cn_right, cn_delta, nearGapOrCentromere, caller) %>%
	distinct() %>%
	left_join(crossmatchdf %>% dplyr::select(sampleId, seqnames, start, cohort, classification, distance) %>% filter(cohort=="gridss pcawg"), by=c("sampleId", "seqnames", "start")) %>%
	left_join(crossmatchdf %>% dplyr::select(sampleId, seqnames, start, cohort, classification, distance) %>% filter(cohort=="pcawg pcawg"), by=c("sampleId", "seqnames", "start"), suffix=c("",".pcawg")) %>%
	filter(classification != classification.pcawg) %>%
	group_by(cohort, classification, classification.pcawg) %>%
	summarise(total=n()) %>%
	View()

# What's wrong with DO1006? So many missing GRIDSS SVs. Is it tumour contamination in the normal?
z = unlist(GRangesList(lapply(icgc_sv_cn_transitions_list, function(x) x$cn_transitions)))
z = z[z$caller == "gridss purple" & is.na(z$distance)]
z = resize(z, width=50, fix="center", ignore.strand=TRUE)
start(z) = pmax(1, start(z))
export(con="missing_gridss_calls.bed", z)












