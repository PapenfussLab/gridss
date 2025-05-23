source("libbenchmark.R")
#
# Probe validation results
#
rawprobeResult = read_csv(paste0(figdir, "supptable_probe_validation_results.csv"))
probeResult = rawprobeResult %>% filter(is.na(exclusion))

rbind(
	probeResult %>% filter(scope!="SharedBoth"),
	probeResult %>% filter(scope=="SharedBoth") %>% mutate(scope="SharedManta"),
	probeResult %>% filter(scope=="SharedBoth") %>% mutate(scope="SharedStrelka")) %>%
	group_by(sampleId, callset, scope, supported) %>%
	summarise(n=n()) %>%
	mutate(category=factor(paste(callset, scope),
												 levels=c("Manta Private", "Gridss SharedManta", "Gridss Private", "Gridss SharedStrelka", "Strelka Private"),
												 labels=c("manta", "gridss+manta", "gridss", "gridss+strelka", "strelka (32bp+)"))) %>%
ggplot() +
	aes(x=sampleId, fill=supported, y=n) +
	geom_bar(stat="identity") +
	facet_grid(category ~ ., scales="free") +
	scale_fill_manual(values=c(gridss_fig_fp_colours[2], gridss_fig_tp_colours[1])) +
	scale_y_continuous(expand=expansion(mult=c(0, 0.02))) +
	labs(fill="Validated", x="", y="breakpoint calls") +
	theme(axis.text.x = element_text(angle = 90))
figsave("probe_results_by_sample", width=5, height=6)


# GRIDSS vs manta
gridss_manta_over50_df = rbind(
	probeResult %>% filter(scope!="SharedBoth"),
	probeResult %>% filter(scope=="SharedBoth") %>% mutate(scope="SharedManta"),
	probeResult %>% filter(scope=="SharedStrelka") %>% mutate(callset="Gridss", scope="Private")) %>%
	mutate(category=factor(paste(callset, scope),
												 levels=c("Gridss Private", "Gridss SharedManta", "Manta Private"),
												 labels=c("gridss", "shared", "manta"))) %>%
	mutate(under50bp=!is.na(length) & abs(length) <= 50) %>%
	filter(!under50bp & !is.na(category))
ggplot(gridss_manta_over50_df %>% replace_na(list(length=1))) +
	aes(x=abs(length), fill=supported) +
	geom_histogram() +
	scale_x_log10() +
	facet_grid(category ~ type) +
	labs(title="Breakdown of SVs over 50bp") +
	theme(axis.text.x = element_text(angle = 90))
probe_results_vs_manta_over_50bp = gridss_manta_over50_df %>%
	group_by(category, supported) %>%
	summarise(n=n()) %>%
	ggplot() +
	aes(x=category, fill=supported, y=n) +
	geom_bar(stat="identity") +
	scale_y_continuous(expand=expand_scale(mult=c(0, 0.02))) +
	scale_fill_manual(values=c(gridss_fig_fp_colours[2], gridss_fig_tp_colours[1])) +
	labs(fill="Probe\nValidated", x="", y="Variant calls (50bp+)") +
	theme(legend.position="none") +
	theme(axis.text.x = element_text(angle = 90))
probe_results_vs_manta_over_50bp
figsave("probe_results_vs_manta_over_50bp", width=3, height=4)

plot_probe_results_vs_stelka_under_50bp = rbind(
	probeResult %>% filter(scope!="SharedBoth"),
	probeResult %>% filter(scope=="SharedBoth") %>% mutate(scope="SharedStrelka"),
	probeResult %>% filter(scope=="SharedManta") %>% mutate(callset="Gridss", scope="Private")) %>%
	mutate(category=factor(paste(callset, scope),
												 levels=c("Gridss Private", "Gridss SharedStrelka", "Strelka Private"),
												 labels=c("gridss", "shared", "strelka"))) %>%
	mutate(under50bp=!is.na(length) & abs(length) <= 50) %>%
	filter(under50bp & !is.na(category)) %>%
	group_by(category, supported) %>%
	summarise(n=n()) %>%
	ggplot() +
	aes(x=category, fill=supported, y=n) +
	geom_bar(stat="identity") +
	scale_y_continuous(expand=expand_scale(mult=c(0, 0.02))) +
	scale_fill_manual(values=c(gridss_fig_fp_colours[2], gridss_fig_tp_colours[1])) +
	labs(fill="Probe\nValidated", x="", y="Variant calls (32-50bp)") +
	theme(legend.position=c(0.6, 0.6)) +
	theme(axis.text.x = element_text(angle = 90))
plot_probe_results_vs_stelka_under_50bp
figsave("probe_results_vs_stelka_under_50bp", width=3, height=4)

sumdf = probeResult %>%
	group_by(callset, scope) %>%
	summarise(calls=n(), validated=sum(supported)) %>%
	ungroup()
rbind(
	sumdf %>%
		filter(str_detect(paste0(callset,scope), "Manta")) %>%
		summarise(calls=sum(calls), validated=sum(validated)) %>%
		mutate(caller="Manta"),
	sumdf %>%
		filter(str_detect(paste0(callset,scope), "Gridss")) %>%
		summarise(calls=sum(calls), validated=sum(validated)) %>%
		mutate(caller="gridss"),
	sumdf %>%
		filter(str_detect(paste0(callset,scope), "Strelka")) %>%
		summarise(calls=sum(calls), validated=sum(validated)) %>%
		mutate(caller="Strelka")) %>%
	mutate(
		prec=validated/calls,
		fdr=1-validated/calls)

shortsumdf = probeResult %>%
	mutate(under50bp=!is.na(length) & abs(length) <= 50) %>%
	filter(under50bp) %>%
	group_by(callset, scope) %>%
	summarise(calls=n(), validated=sum(supported)) %>%
	ungroup()
rbind(
	shortsumdf %>%
		filter(str_detect(paste0(callset,scope), "Manta")) %>%
		summarise(calls=sum(calls), validated=sum(validated)) %>%
		mutate(caller="Manta"),
	shortsumdf %>%
		filter(str_detect(paste0(callset,scope), "Gridss")) %>%
		summarise(calls=sum(calls), validated=sum(validated)) %>%
		mutate(caller="gridss"),
	shortsumdf %>%
		filter(str_detect(paste0(callset,scope), "Strelka")) %>%
		summarise(calls=sum(calls), validated=sum(validated)) %>%
		mutate(caller="Strelka")) %>%
	mutate(
		prec=validated/calls,
		fdr=1-validated/calls)

# 32-100bp DUP
rbind(
	probeResult %>% filter(scope!="SharedBoth"),
	probeResult %>% filter(scope=="SharedBoth") %>% mutate(scope="SharedManta"),
	probeResult %>% filter(scope=="SharedBoth") %>% mutate(scope="SharedStrelka")) %>%
	filter(type == "DUP" & abs(endPosition - startPosition) <= 100 & abs(endPosition - startPosition) >= 32) %>%
	group_by(callset, scope, supported) %>%
	summarise(n=n()) %>%
	# hack to ensure even 0 event counts are plotted
	ungroup() %>%
	rbind(probeResult %>% dplyr::select(callset, scope, supported) %>% distinct() %>% mutate(n=0)) %>%
	group_by(callset, scope, supported) %>%
	summarise(n=sum(n)) %>%
	filter(!(callset== "Gridss" & scope=="SharedBoth" & n == 0)) %>%
	mutate(category=factor(paste(callset, scope),
												 levels=c("Manta Private", "Gridss SharedManta", "Gridss Private", "Gridss SharedStrelka", "Strelka Private"),
												 labels=c("manta", "gridss+manta", "gridss", "gridss+strelka", "strelka"))) %>%
	ggplot() +
	scale_y_continuous(expand=expand_scale(mult=c(0, 0.02))) +
	aes(x=category, fill=supported, y=n) +
	geom_bar(stat="identity") +
	scale_fill_manual(values=c(gridss_fig_fp_colours[2], gridss_fig_tp_colours[1])) +
	labs(fill="Validated", x="", y="Duplications (32-100bp)") +
	theme(axis.text.x = element_text(angle = 90))
figsave("probe_results_32-100bp_DUP", width=5, height=4)

# Single breakend probe support
probeResult %>%
	filter(is.na(endChromosome)) %>%
	ggplot() +
	aes(x=scope, fill=supported) +
	geom_bar()

probeResult %>%
	filter(is.na(endChromosome)) %>%
	ungroup() %>%
	summarise(
		total=n(),
		supported=sum(supported),
		pct=supported/total * 100)

# export to fasta so we can run RM to identify single breakend repeat type
sglprobes = probeResult %>%
	filter(is.na(endChromosome))
sglfasta = DNAStringSet(x=sglprobes$insertSequence)
names(sglfasta) = paste0(sglprobes$uniqueId, "_", sglprobes$supported)
export(sglfasta, "breakend_probe_sequences.fasta")	
#RepeatMasker -species human -dir probe_repeat_masker breakend_probe_sequences.fasta

###
if (!exists("db")) db = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
samplelookup = DBI::dbGetQuery(db, "Select hmfSampleId, SampleId from clinical WHERE hmfSampleId IN ('HMF000214A','HMF001434A','HMF000566A','HMF000823A','HMF000952A','HMF003333A','HMF001559A','HMF002624A','HMF002224A','HMF000379A','HMF002760A','HMF003228A','HMF001184A')")
lnx_drivers = read_csv(paste0(privatedatadir, "/hartwig/LNX_DRIVERS.csv"))
lnx_disruptions = read_csv(paste0(privatedatadir, "/hartwig/LNX_DISRUPTIONS.csv"))
lnx_svs = read_csv(paste0(privatedatadir, "/hartwig/LNX_SVS.csv"))
emblgenes=read_csv = read_csv(paste0(datadir, "ensembl_gene_data.csv"), col_types="cccnnncc")

probe_sample_disruptions = lnx_disruptions %>%
	inner_join(samplelookup, by=c("SampleId"="sampleId")) %>%
	inner_join(emblgenes, by=c("GeneId", "GeneName")) %>%
	dplyr::select(SampleId, Reportable, GeneName, GeneStart, GeneEnd) %>%
	distinct()
probe_sample_drivers = lnx_drivers %>%
	inner_join(samplelookup, by=c("SampleId"="sampleId")) %>%
	inner_join(emblgenes, by=c("Gene"="GeneName"))

svs_in_driver_events = lnx_svs %>%
	inner_join(samplelookup, by=c("SampleId"="sampleId")) %>%
	inner_join(lnx_drivers, by=c("SampleId", "ClusterId"))

driver_sv_links = inner_join(svs_in_driver_events, probeResult, by=c("hmfSampleId"="sampleId", "ChrStart"="startChromosome", "PosStart"="startPosition"))
View(driver_sv_links %>%
	group_by(hmfSampleId, Gene, callset, scope) %>%
	summarise(count=n()) %>%
	mutate(callset_scope=paste(callset, scope)) %>%
	dplyr::select(-callset, -scope) %>%
	spread(callset_scope, count))

View(driver_sv_links %>%
		 	group_by(hmfSampleId, Gene, DriverType, EventType, callset, scope) %>%
		 	summarise(count=n()) %>%
		 	mutate(callset_scope=paste(callset, scope)) %>%
		 	dplyr::select(-callset, -scope) %>%
		 	spread(callset_scope, count))

driver_sv_links %>%
	group_by(hmfSampleId, Gene, ClusterId, ClusterCount.x, callset, scope) %>%
	summarise(count=n()) %>%
	mutate(callset_scope=paste0(callset, "_", scope)) %>%
	dplyr::select(-callset, -scope) %>%
	spread(callset_scope, count) %>% 
	group_by(hmfSampleId, ClusterId, ClusterCount.x, Gridss_Private, Gridss_SharedManta, Manta_Private) %>%
	summarise(genes=n(), gene_names=paste(Gene, collapse=", ")) %>%
	View()

driver_sv_links %>% filter(callset=="Manta" & scope =="Private") %>% View()



