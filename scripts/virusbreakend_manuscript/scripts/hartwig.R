library(tidyverse)
library(ggExtra)
library(VariantAnnotation)
theme_set(theme_bw())
setwd("../virusbreakend_manuscript/scripts")
if (!file.exists("../protecteddata/hmf_cancertype.csv")) {
	library(RMySQL)
	if (!exists("db")) db = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
	cancertypedf = DBI::dbGetQuery(db, "select sampleId, primaryTumorLocation, primaryTumorSubLocation, primaryTumorType, primaryTumorSubType, primaryTumorExtraDetails FROM clinical")
	write_tsv(cancertypedf, "../protecteddata/hmf_cancertype.csv")
} else {
	cancertypedf = read_tsv("../protecteddata/hmf_cancertype.csv")
}
# krakendf = read_tsv(
# 		"../protecteddata/hmf/viral.out",
# 		col_name=c("file", "pct", "treereads", "directreads", "level", "taxid", "name"),
# 		col_types="cniicic") %>%
# 	mutate(sample=str_extract(str_replace(file, stringr::fixed("./"), ""), "^[^.]+"))
# extracteddf = read_tsv(
# 		"../protecteddata/hmf/extracted.out",
# 		col_name=c("file", "pct", "treereads", "directreads", "level", "taxid", "name"),
# 		col_types="cniicic") %>%
# 	mutate(sample=str_extract(str_replace(file, stringr::fixed("./"), ""), "^[^.]+"))
sampledf = read_tsv("../protecteddata/hmfv3/samples.csv", col_names="sample", col_types="c")
cancertypedf = cancertypedf %>% filter(sampleId %in% sampledf$sample)
cancertype_summary = cancertypedf %>% group_by(primaryTumorLocation) %>% summarise(n=n())
summarydf = read_tsv("../protecteddata/hmfv3/merged.vcf.summary.csv") %>%
	left_join(cancertypedf, by=c("sample"="sampleId")) %>%
	mutate(hasIntegration=integrations > 0)
names(summarydf)[6] = "name_species"
names(summarydf)[7] = "reads_species"
mergedvcf = readVcf("../protecteddata/hmfv3/merged.vcf")
rownames(mergedvcf) = NULL
alt = str_match(as.character(rowRanges(mergedvcf)$ALT), "([^\\[\\]]*)([\\[\\]])adjusted_kraken_taxid_([0-9]+)_(.*):([0-9]+)[\\[\\]]([^\\[\\]]*)")
colnames(alt) = c("ALT", "leftins", "ori", "taxid", "virus_chr", "virus_pos", "rightins")
sitedf = as.data.frame(alt) %>%
	mutate(
		host_chr=as.character(seqnames(mergedvcf)),
		host_pos=start(mergedvcf),
		QUAL=rowRanges(mergedvcf)$QUAL,
		FILTER=rowRanges(mergedvcf)$FILTER,
		taxid=as.integer(taxid)) %>%
	bind_cols(as.data.frame(info(mergedvcf))) %>%
	mutate(
		sample=SAMPLE,
		bealn=unlist(lapply(BEALN, function(x) paste0(x, collapse=";"))),
		mapq=as.integer(str_extract(unlist(lapply(BEALN, function(x) x[[1]])), "[0-9]+$")),
		bealn_hits=elementNROWS(BEALN))

friendly_species_name = c(
	"Adeno-associated dependoparvovirus A"="AAV-2",
	"Alphapapillomavirus 7"="HPV-18",
	"Alphapapillomavirus 9"="HPV-16",
	"Gammapapillomavirus 16"="HPV", #-137",
	"Gammapapillomavirus 19"="HPV", #"HPV-166",
	"Gammapapillomavirus 8"="HPV", #"HPV-112"
	"Gammapapillomavirus 9"="HPV", #"HPV-112"
	"Betapapillomavirus 4"="HPV",
	"Torque teno virus"="Other",
	"Torque teno virus 1"="Other",
	"Torque teno virus 6"="Other",
	"Torque teno virus 9"="Other",
	"Human alphaherpesvirus 1"="HSV-1",
	"Human alphaherpesvirus 3"="HHV-3",
	"Human betaherpesvirus 5"="HHV-5", # aka HCMV
	"Human betaherpesvirus 6A"="HHV-6",
	"Human betaherpesvirus 6B"="HHV-6",
	"Human betaherpesvirus 7" = "HHV-7",
	"Human gammaherpesvirus 4"="EBV",
	"Human gammaherpesvirus 8"="HHV-8",
	"Human polyomavirus 1"="BKPyV",
	"Human polyomavirus 5"="MCPyV",
	"Human polyomavirus 6"="HPyV6",
	"Human polyomavirus 7"="HPyV7",
	"Hepatitis B virus"="HBV",
	"Baboon endogenous virus"="Other",
	"Primate erythroparvovirus 1"="Other",
	"Porcine type-C oncovirus"="Other",
	"Vaccinia virus"="Other",
	"Human erythrovirus V9"="Other")

friendly_summary_species_name = c(
	"Adeno-associated dependoparvovirus A"="Other",
	"Alphapapillomavirus 7"="HPV",
	"Alphapapillomavirus 9"="HPV",
	"Gammapapillomavirus 16"="HPV", #-137",
	"Gammapapillomavirus 19"="HPV", #"HPV-166",
	"Gammapapillomavirus 8"="HPV", #"HPV-112"
	"Gammapapillomavirus 9"="HPV", #"HPV-112"
	"Betapapillomavirus 4"="HPV",
	"Human alphaherpesvirus 1"="HSV",
	"Human alphaherpesvirus 3"="HHV",
	"Human betaherpesvirus 5"="HHV", # aka HCMV
	"Human betaherpesvirus 6A"="HHV",
	"Human betaherpesvirus 6B"="HHV",
	"Human betaherpesvirus 7" = "HHV",
	"Human gammaherpesvirus 4"="EBV",
	"Human gammaherpesvirus 8"="HHV",
	"Human polyomavirus 1"="BKPyV",
	"Human polyomavirus 5"="MCPyV",
	"Human polyomavirus 6"="HPyV",
	"Human polyomavirus 7"="HPyV",
	"Hepatitis B virus"="HBV",
	"Torque teno virus"="Other",
	"Torque teno virus 1"="Other",
	"Torque teno virus 6"="Other",
	"Torque teno virus 9"="Other",
	"Baboon endogenous virus"="Other",
	"Primate erythroparvovirus 1"="Other",
	"Porcine type-C oncovirus"="Other",
	"Vaccinia virus"="Other",
	"Human erythrovirus V9"="Other")

# Filter:
exclude_from_analysis = c(
	44561, # Murine type C retrovirus
	11764, # Baboon endogenous virus strain M7
	168238, #	Porcine endogenous retrovirus E
	11757, # Mouse mammary tumor virus
	summarydf %>% filter(str_detect(name_species, "Torque teno")) %>% pull(taxid)
)
summarydf = summarydf %>% filter(!(taxid %in% exclude_from_analysis))
sitedf = sitedf %>% filter(!(taxid %in% exclude_from_analysis))
	
sampledf %>%
	mutate(
		hasVirus=sample %in% summarydf$sample,
		hasIntegration=sample %in% sitedf$sample) %>%
	summarise(
		n=n(),
		hasVirus=sum(hasVirus),
		hasIntegration=sum(hasIntegration))

summarydf %>%
	group_by(name_genus, name_species) %>%
	summarise(
		n=n(),
		withIntegration=sum(hasIntegration),
		maxcoverage=max(coverage),
		meancoverage=mean(coverage),
		mediandepth=median(meandepth),
		avgdepth=mean(meandepth)) %>%
	View()

summarydf %>% filter(meandepth > 10) %>% summarise(n=n(), withIntegration=sum(hasIntegration))
summarydf %>% filter(!hasIntegration) %>% View()

# What's up with Baboon endogenous virus?
# 5/6 insertions are into the exact same location in BCL7C, all sub-clonal
sitedf %>%
	filter(taxid %in% (summarydf %>% filter(name_species=="Baboon endogenous virus") %>% pull(taxid))) %>% 
	inner_join(cancertypedf, by=c("sample"="sampleId"))

p = ggplot(summarydf) +
	aes(x=coverage, y=meandepth, colour=integrations > 0, shape=integrations > 0) +
	geom_point() +
	labs(title="Viral coverage") +
	scale_y_log10()
ggMarginal(p, margins="y", type="histogram", groupFill=TRUE)

ggplot(summarydf) +
	aes(x=coverage, y=meandepth, colour=integrations > 0, shape=integrations > 0) +
	geom_point() +
	labs(title="Viral coverage") +
	facet_wrap(~friendly_species_name[name_species]) +
	scale_y_log10()

# what about the non-clonal HPV-16 integrations?
summarydf %>% filter( friendly_species_name[name_species] == "HPV-16") %>% group_by(primaryTumorLocation) %>% filter(integrations==0, meandepth<10) %>% pull(sample)

plot_hmf_depth = ggplot(summarydf) +
	aes(x=meandepth, fill=integrations > 0) +
	geom_histogram(bins=50) +
	scale_x_log10(expand=c(0,0), breaks=10^(-1:5), labels=as.character(10^(-1:5))) +
	scale_y_continuous(expand=c(0,0), limits=c(0, 50)) +
	coord_cartesian(xlim=c(0.1, 50000)) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
	labs(x="Viral genome coverage", y="samples", fill="Viral\nIntegration\nDetected")

# Viral prevalence by cancer type
fulldf = expand.grid(sample=cancertypedf$sampleId, friendly_summary_species_name=unique(friendly_summary_species_name)) %>%
	inner_join(cancertypedf, by=c("sample"="sampleId")) %>%
	left_join(summarydf %>% mutate(friendly_summary_species_name=friendly_summary_species_name[name_species])) %>%
	mutate(
		status=ifelse(is.na(integrations), "Not found", ifelse(integrations > 0, "Integrated", "Found")),
		primary=paste(primaryTumorLocation, primaryTumorSubLocation))
ggplot(fulldf) +
	aes(x=friendly_summary_species_name, fill=status) +
	geom_bar() +
	facet_wrap(~ primary, scales="free") +
	theme(axis.text.x = element_text(angle = 90)) +
	labs(fill="Viral\nIntegration\nFound", x="Virus", y="")

ggplot(summarydf %>% filter(name_genus=="Alphapapillomavirus")) +
	aes(x=meandepth, fill=integrations > 0) +
	facet_wrap(~ name_species) +
	geom_histogram(bins=50) +
	scale_x_log10() +
	labs(title="Alphapapillomavirus viral coverage")


genusfulldf = expand.grid(sample=cancertypedf$sampleId, name_genus=unique(summarydf$name_genus)) %>%
	inner_join(cancertypedf, by=c("sample"="sampleId")) %>%
	left_join(summarydf) %>%
	mutate(
		status=ifelse(is.na(integrations), "Not found", ifelse(integrations > 0, "Integrated", "Found")),
		primary=paste(primaryTumorLocation, primaryTumorSubLocation))

genusfulldf %>%
	group_by(primary, name_genus) %>%
	summarise(
		n=n(),
		found=sum(status != "Not found"),
		found_integrated=sum(status=="Integrated")) %>%
	filter(
		(name_genus =="Orthohepadnavirus" & str_detect(primary, "Liver ")) | 
		(name_genus =="Alphapapillomavirus" & str_detect(primary, "Cervix"))
	)
# Cervical samples without HPV
cervical = cancertypedf %>% filter(primaryTumorSubLocation=="Cervix") %>%
	left_join(summarydf %>% filter(name_genus=="Alphapapillomavirus"), by=c("sampleId"="sample"))
cervical_with_hpv = cervical %>% filter(!is.na(integrations))
cervical_without_hpv = cervical %>% filter(is.na(integrations)) %>% pull(sampleId)
summarydf %>% filter(sample %in% cervical_without_hpv)

ggplot(summarydf) +
	aes(x=reads_species, fill=integrations > 0) +
	geom_bar() +
	theme(axis.text.x = element_text(angle = 90)) +
	labs(fill="Viral\nIntegration\nFound", x="Virus", y="")

#####
# Gene annotation helpers
collapse_sample_nearby = function(gr, maxgap) {
	hits = findOverlaps(gr, gr, maxgap=100000) %>%
		as.data.frame() %>%
		filter(queryHits != subjectHits) %>%
		filter(gr$sample[queryHits] == gr$sample[subjectHits])
	end(gr)[hits$subjectHits] = pmax(end(gr)[hits$subjectHits], end(gr)[hits$queryHits])
	start(gr)[hits$subjectHits] = pmin(start(gr)[hits$subjectHits], start(gr)[hits$queryHits])
	gr = gr[-(hits %>% filter(subjectHits >= queryHits) %>% pull(subjectHits))]
	return(gr)
}
nearby_genes = function(gr, distance, txdb=TxDb.Hsapiens.UCSC.hg19.knownGene) {
	library(TxDb.Hsapiens.UCSC.hg19.knownGene)
	library(org.Hs.eg.db)
	gr$nearby=""
	hitdf = findOverlaps(gr, genes(txdb), maxgap=distance, ignore.strand=TRUE) %>%
		as.data.frame() %>%
		mutate(
			gene_id=genes(txdb)$gene_id[subjectHits],
			symbol=AnnotationDbi::select(org.Hs.eg.db, keys=gene_id, columns="SYMBOL", keytype="ENTREZID")$SYMBOL) %>%
		group_by(queryHits) %>%
		summarise(symbols=paste(symbol, collapse=","))
	gr$nearby[hitdf$queryHits] = hitdf$symbols
	return(gr$nearby)
}
ann_nearby = function(gr) {
	seqlevelsStyle(gr) = "UCSC"
	gr$gene = nearby_genes(gr, 0)
	gr$gene10k = nearby_genes(gr, 10000)
	gr$gene100k = nearby_genes(gr, 100000)
	gr$gene250k = nearby_genes(gr, 250000)
	return(gr)
}

######
# HBV integrations
hbv_gr = with(
		sitedf %>%
			inner_join(summarydf, by=c("sample", "taxid")) %>%
			filter(name_genus=="Orthohepadnavirus"),
		GRanges(seqnames=host_chr, ranges=IRanges(start=as.numeric(host_pos), width=1), sample=sample)) %>%
	collapse_sample_nearby(100000) %>%
	ann_nearby()
tert_hbv_gr = hbv_gr[overlapsAny(hbv_gr, GRanges(seqnames=5, ranges=IRanges(start=1253287-300000, end=1295162+300000)))]
non_tert_hbv_gr=hbv_gr[!(hbv_gr$sample %in% tert_hbv_gr$sample)]
View(non_tert_hbv_gr %>% as.data.frame() %>% inner_join(cancertypedf, by=c("sample"="sampleId")))

######
# Merkel cell
summarydf %>% filter(name_species == "Human polyomavirus 5") %>% pull(sample)
merkeldf = summarydf %>%
	filter(name_species == "Human polyomavirus 5" |
		str_detect(paste(primaryTumorSubType, primaryTumorType, primaryTumorLocation, primaryTumorSubLocation), "Merkel"))
merkel_gr = with(
	sitedf %>% filter(sample %in% merkeldf$sample),
	GRanges(seqnames=host_chr, ranges=IRanges(start=as.numeric(host_pos), width=1), sample=sample)) %>%
	collapse_sample_nearby(100000) %>%
	ann_nearby()
View(merkel_gr %>% as.data.frame())

cancertypedf %>%
	filter(str_detect(primaryTumorSubType, "Merkel")) %>%
	left_join(summarydf, by=c("sampleId"="sample")) %>%
	View()

######
# Human gammaherpesvirus 8

kaposidf = summarydf %>%
	filter(name_species == "Human gammaherpesvirus 8" |
				 	str_detect(paste(primaryTumorSubType, primaryTumorType, primaryTumorLocation, primaryTumorSubLocation), "aposi"))

######
# BKPyV 
BKPyVdf = summarydf %>%
	filter(name_species == "Human polyomavirus 1" |
				 	str_detect(paste(primaryTumorSubType, primaryTumorType, primaryTumorLocation, primaryTumorSubLocation), "zzzz"))
BKPyVdf %>% pull(sample)

######
# Repeat distributions
not_extracted_viruses = krakendf %>%
	filter(!(taxid %in% c(1, 10239, 28883, 687329, 40272, 1507401))) %>%
	filter(directreads >= 50) %>%
	group_by(sample) %>%
	filter(any(!extracted) & any(extracted)) %>%
	filter(length(unique(str_trunc(name, 20))) > 1) %>%
	arrange(sample, desc(extracted), desc(directreads), name)

shortRepeatClass = function(x) ifelse(x %in% c("Satellite/centr"), x, str_extract(x, "^[^/]+"))
friendlyRepeatClass = function(x) sapply(x, function(xx) factor(ifelse(is.na(xx), "No repeat", switch(str_extract(xx, "^[^/]+"),
	"Low_complexity"="Simple",
	"LINE"="LINE",
	"Simple_repeat"="Simple",
	"Satellite"="Satellite",
	"SINE"="SINE",
	"DNA"="DNA",
	"LTR"="LTR",
	"rRNA"="Other",
	"Retroposon"="Other",
	"srpRNA"="Other",
	"Unknown"="Other",
	"No repeat"="No repeat",
	"No repeat")),
	levels=c("DNA", "LINE", "LTR", "SINE", "Satellite", "Simple", "Other", "No repeat")))

ggplot(sitedf %>% inner_join(summarydf, by=c("sample", "taxid"))) +
	aes(x=mapq, fill=shortRepeatClass(INSRMRC)) +
	facet_wrap(~ name_genus) +
	geom_histogram(bins=6) +
	labs(title="Insertion site assembly mapping quality", y="Sites", x="mapq", fill="Repeat Class")

ggplot(sitedf %>% inner_join(summarydf, by=c("sample", "taxid"))) +
	aes(x=mapq, fill=name_genus) +
	geom_histogram(bins=6) +
	labs(title="Insertion site assembly mapping quality", y="Sites", x="mapq", fill="Repeat Class")

plot_rmmap = sitedf %>%
	inner_join(summarydf, by=c("sample", "taxid")) %>%
	mutate(mappable=ifelse(mapq>10, "Mappable", "Unmappable")) %>%
	ggplot() +
	aes(x=friendly_summary_species_name[name_species], fill=friendlyRepeatClass(INSRMRC)) +
	geom_bar() +
	facet_wrap(~ mappable) +
	scale_y_continuous(expand=c(0,0)) +
	scale_fill_brewer(palette="Dark2") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
	theme(axis.text.x = element_text(angle = 90)) +
	labs(x="", y="Viral integrations", fill="Repeat")

plot_grid(plot_hmf_depth, plot_rmmap, labels=c("a", "b"))
ggsave("figure3.pdf", height=4, width=10)
ggsave("figure3.png", height=4, width=10)

# Telomeric insertions
sitedf %>% filter(INSRMRT=="(TAACCC)n") %>%
	group_by(sample) %>%
	summarise(sites=n()) %>%
	inner_join(cancertypedf, by=c("sample"="sampleId")) %>%
	group_by(primaryTumorLocation) %>%
	summarise(
		infected_samples=n(),
		sites=sum(sites)) %>%
	left_join(cancertypedf %>%
		filter(sampleId %in% krakendf$sample) %>%
		group_by(primaryTumorLocation) %>%
		summarise(cohort_size=n()))
# Centromeric insertions
sitedf %>%
	group_by(taxid, sample) %>%
	mutate(isSatellite=INSRMRC %in% c("Satellite", "Satellite/centr")) %>%
	filter(any(isSatellite)) %>%
	summarise(satellite_sites=sum(isSatellite), total_sites=n()) %>%
	inner_join(cancertypedf, by=c("sample"="sampleId"))

summarydf %>% inner_join(
	sitedf %>% filter(sample %in%
		(sitedf %>% group_by(sample) %>%
			summarise(
				hasMappable=any(mapq>=10),
				hasUnmappable=any(mapq<10)) %>%
			filter(hasUnmappable & !hasMappable) %>%
			inner_join(summarydf) %>%
			filter(hasIntegration) %>%
			filter(name_genus == "Alphapapillomavirus") %>%
			pull(sample)))) %>%
	View()

######
# Cancer type associations
by_label_summary = function(summarydf, typedf, category) {
	#typedf = typedf %>% filter(sampleId %in% summarydf$sample)
	summarydf %>% 
		mutate(virus=friendly_summary_species_name[name_species]) %>%
		replace_na(list(virus="Other")) %>%
		group_by(label, virus) %>%
		summarise(
			found=n(),
			integrated=sum(integrations > 0)) %>%
		inner_join(typedf %>% group_by(label) %>% summarise(samples=n())) %>%
		mutate(category=category)
}
pancancerratedf = by_label_summary(
	summarydf %>% mutate(label="pancancer"),
	cancertypedf %>% mutate(label="pancancer"),
	"pancancer")
typeratedf = bind_rows(
	by_label_summary(
		summarydf %>% mutate(label=primaryTumorLocation),
		cancertypedf %>% mutate(label=primaryTumorLocation),
		"primaryTumorLocation"),
	by_label_summary(
		summarydf %>% mutate(label=paste(primaryTumorLocation, primaryTumorSubLocation)),
		cancertypedf %>% mutate(label=paste(primaryTumorLocation, primaryTumorSubLocation)),
		"primaryTumorSubLocation"),
	by_label_summary(
		summarydf %>% mutate(label=primaryTumorType),
		cancertypedf %>% mutate(label=primaryTumorType),
		"primaryTumorType"),
	by_label_summary(
		summarydf %>% mutate(label=paste(primaryTumorType, primaryTumorSubType)),
		cancertypedf %>% mutate(label=paste(primaryTumorType, primaryTumorSubType)),
		"primaryTumorSubType"))
require(R.cache)
memoized.fisher.test <- addMemoization(fisher.test)
rowwise_fisher.test <- function(df,
																hits_a_column, total_a_column, a_name,
																hits_b_column, total_b_column, b_name) {
	df %>% bind_cols(
		data.frame(
			a_hits=df[[hits_a_column]],
			b_hits=df[[hits_b_column]],
			total_a=df[[total_a_column]],
			total_b=df[[total_b_column]],
			per_row_hack=seq_len(nrow(df))) %>%
			group_by(per_row_hack) %>%
			do({
				contingency_table <- data.frame(
					hits=c(.$a_hits, .$b_hits),
					misses=c(.$total_a - .$a_hits, .$total_b - .$b_hits))
				test_obj <- memoized.fisher.test(contingency_table)
				df <- data.frame(
					estimate_odds_ratio=test_obj$estimate,
					estimate_odds_ratio_95_lower=test_obj$conf.int[1],
					estimate_odds_ratio_95_upper=test_obj$conf.int[2],
					p_value=test_obj$p.value)
				df[[paste0(a_name, "_likelihood")]] <- .$a_hits / .$total_a
				df[[paste0(b_name, "_likelihood")]] <- .$b_hits / .$total_b
				df$likelihood_ratio = (.$a_hits / .$total_a) / (.$b_hits / .$total_b)
				df
			}) %>%
			ungroup() %>%
			dplyr::select(-per_row_hack)) %>%
		mutate(adjusted_p_value=p.adjust(p_value))
}
rate_dela = typeratedf %>%
	inner_join(pancancerratedf, by="virus", suffix=c("", ".pancancer")) %>%
	rowwise_fisher.test("found", "samples", "label", "found.pancancer", "samples.pancancer", "pancancer")
write_csv(path="cancer_type_associations.csv", rate_dela %>% arrange(p_value))
#rate_dela %>% arrange(p_value) %>% View()

ggplot(pancancerratedf) +
	aes(x=virus) +
	geom_bar(stat="identity", aes(y=found), colour="blue") +
	geom_bar(stat="identity", aes(y=integrated), colour="green") +
	theme(axis.text.x = element_text(angle = 90))

######
#
genus_cohort_summary = krakendf %>%
	mutate(hasIntegration=sample %in% sitedf$sample) %>%
	filter(treereads >= 50) %>%
	filter(str_detect(level, "G")) %>%
	inner_join(cancertypedf, by=c("sample"="sampleId")) %>%
	full_join(cancertypedf %>%
							filter(sampleId %in% krakendf$sample) %>%
							group_by(primaryTumorLocation) %>%
							summarise(cohort_size=n())) %>%
	group_by(name, primaryTumorLocation, cohort_size) %>%
	summarise(integration=sum(hasIntegration), virus=n())

genus_cohort_summary %>%
	group_by(name) %>%
	summarise(integration=sum(integration), virus=sum(virus)) %>%
	mutate(pct=100*integration/virus) %>%
	arrange(desc(pct))

ggplot(genus_cohort_summary %>% filter(name=="Alphapapillomavirus")) +
	aes(x=primaryTumorLocation) +
	geom_point(aes(y=virus)) +
	geom_point(aes(y=integration, colour="green")) + 
labs(title="HPV", x="primary tumour location")


extracteddf %>%
	mutate(integration=sample %in% sitedf$sample) %>%
	group_by(name) %>%
	summarise(integration=sum(integration), virus=n())

extracteddf %>%
	mutate(integration=sample %in% sitedf$sample) %>%
	inner_join(wgsdf, by="sample") %>%
	filter(name %in% c("Merkel cell polyomavirus", "Adeno-associated virus - 2", "Alphapapillomavirus 7", "Human gammaherpesvirus 4", "", "", "", "", "", "", "", "", "") | str_detect(name, "Human papillomavirus") | str_detect(name, "HBV")) %>%
ggplot() +
	aes(x=MEAN_COVERAGE, fill=integration) +
	scale_x_log10() +
	geom_histogram() +
	labs(title="Viral coverage", x="Mean coverage")



# Virus-Host-DB filter
vhdb_human = c(
	10804,46350,57579,82300,68558,68742,335103,338079,1313215,172148,333767,337043,337042,2169991,568715,858367,683172,1247114,767521,645687,
	683174,1247113,90961,622416,649604,38837,11020,1987017,1803394,1980459,341053,1263720,2169992,1391667,565995,35304,35310,2169993,94433,
	1678225,2170194,10317,50292,11272,629725,499556,99565,37124,466217,1346815,1346816,46839,1330491,2003650,2003651,2003652,1582094,10243,
	42779,42769,42770,42771,42772,42773,31704,42775,42776,42782,42783,12089,42785,42786,86107,42788,12067,12071,82639,12072,12073,12074,74561,
	1980519,1867125,742922,742923,742915,742916,742917,742918,742919,1520935,742925,742924,1348500,12637,11053,11060,11069,11070,114654,11318,
	1980521,38767,2170195,11021,128952,103914,103915,12078,35293,47501,47502,47503,47504,47505,47506,47507,35294,47508,47509,47510,45101,47511,
	45999,47512,47516,41846,47513,47514,47515,35295,40280,12062,46018,12060,1452750,1987020,1987021,12104,156647,150846,138948,1435148,39054,
	306588,138949,325447,408681,1248136,1442164,171720,200154,222889,282345,222887,318563,318564,318565,318566,318567,318568,318569,318570,318571,
	325446,138950,861519,325445,1295563,152632,138951,42789,310907,123738,57482,57483,273345,64307,2049444,54290,1673681,1862824,1862825,59300,204269,
	45219,80941,2170197,1214955,1163715,1415627,1415628,1980471,93830,489450,489451,489460,489466,489483,1219353,1219352,489484,1219350,1219351,
	489493,1052603,1052604,1216928,63330,2008762,11103,10407,1282342,106820,10408,356391,63746,41856,356114,33745,33746,42182,1544901,31646,31647,
	484894,31649,31651,356466,1208063,1094901,356426,42792,1094895,31653,1094899,1094898,1094900,693426,745709,693427,1208062,31655,438880,413256,
	413257,438881,467339,12475,12461,509628,12092,682382,1960046,943272,1233675,180587,435487,434023,435489,397342,2021738,10533,1461395,343463,
	28282,10521,31544,46922,10528,28278,10515,32608,1521028,46924,46925,46930,459785,260047,10529,10548,10522,46936,10524,46938,46941,218120,28285,
	107462,332179,651580,714978,880565,980087,10534,1069441,1124795,1145295,1094363,1337398,1101371,1096763,10519,1643649,326124,1972755,45659,28275,
	52275,31545,28280,10298,10310,10335,2038728,1985415,2169934,2169935,2169936,1868658,12456,626160,626164,626180,626181,12701,35300,35741,37130,43358,
	1518575,1235996,1298362,1306931,10359,32603,32604,10372,689403,1511882,638313,1511883,1108925,1108926,1108927,394504,430048,444473,158465,258453,
	1525173,11137,290028,277944,31631,1233383,586420,12070,1345637,1904876,166122,166124,627439,1193974,1006008,318562,460162,650039,1230646,1247477,
	310756,310755,72197,2017081,2017082,2017083,1820160,1820158,1820159,11641,10376,37296,1792832,1488574,270642,988776,1704090,10338,12509,12721,11676,
	11709,129875,108098,129951,130310,130308,130309,162145,11250,10566,587349,518628,587350,427343,915426,518629,518630,587351,915428,720708,765052,
	915429,765054,765055,673323,1055684,746832,909329,909331,909332,909333,1070408,1070409,1070411,1070412,1070413,1070414,1070415,1070416,1070418,
	1070419,942038,909328,10606,743812,990302,1195796,1165934,1209820,1420545,10607,1434986,1434987,1434988,1434782,1478160,1472342,1472343,1851130,
	10608,1542134,1682340,1650736,37954,37955,37111,37112,10614,37958,37959,10617,333923,28312,10621,260717,587347,587348,563190,1647924,333759,915425,
	338327,338323,338326,565537,10580,915427,735496,696746,927771,765051,10604,765056,931209,931210,10573,909330,1070410,1070417,1070420,31546,743811,
	1110710,1248396,333760,1315264,1315262,1315259,333761,10583,333751,31547,2093789,31548,37956,10609,333762,333752,325568,10584,10611,10585,333763,
	10586,333764,10587,37957,334514,10588,10615,10589,10590,10591,10592,10593,10594,40538,10616,40539,10595,10618,333765,1671798,37114,10596,333753,
	186160,471358,10598,37115,10599,31552,40540,37116,334210,28311,37119,37120,45240,338322,940838,37121,37122,10600,10620,39457,120686,1484958,1297549,
	10579,333771,129724,333772,150546,652810,171370,120381,337054,333769,211787,247269,338324,1288116,188538,11224,11226,12063,39085,195055,362784,376148,
	411152,602074,289365,1511919,10798,1729141,145856,12080,12081,12083,12086,1891762,1208309,1303334,1891764,746830,746831,943908,988774,253182,12730,
	11216,12129,573824,185892,185946,39767,147684,185893,185894,31708,147674,185895,147675,44128,185896,12135,185897,185898,147672,44129,185901,185902,
	185903,185904,185905,147673,185906,185907,185908,185909,167321,185910,147676,167324,185911,44131,44132,185912,185913,185914,185915,185916,185917,
	44133,185918,185919,185920,44134,185921,185922,44135,185923,185924,185925,44136,185927,185928,185929,185930,185931,185932,147685,185890,185934,
	167322,185937,185939,12132,185891,185940,185941,185942,185943,150904,185899,185900,44130,167329,147680,185889,147683,147714,147677,147681,185926,
	44137,185933,185936,167328,185938,167325,167326,167327,167330,185945,992230,655428,10941,10942,1979160,1595998,11963,743300,11908,36368,11909,
	28332,318279,511755,1930509,1450749,79899,273341,59563,1755290,641501,641809,93838,130760,488241,335341,370519,370522,370128,211044,183764,
	1332244,518987,98832,11553,290008,11072,10632,1419710,64310,423447,423448,423446,1224522,44024,33743,11577,38766,378830,378831,1126254,11085,
	11620,318848,40057,1965344,172317,172316,82658,11086,649188,11623,642022,10325,338478,11628,348013,1239565,1046251,11269,33727,59301,35279,
	11234,170525,170528,1056492,763792,70149,884096,884100,884099,884094,884098,884093,884101,884097,884095,1089795,419697,12107,493803,943084,
	1335626,64300,12538,10280,10244,619591,1516081,1979165,1384672,11168,301186,1048854,334203,11079,1203539,1231709,1497391,121791,122928,1529909,
	1529910,122929,552592,490039,1514999,1529911,1529912,1529913,1529914,1529915,1529916,1529917,1529918,1529919,1529920,1529921,1529922,1529923,
	1529924,1529925,1529926,1529927,1529928,1529929,1529930,1529931,1529932,262897,336619,1383577,1184498,1383578,1383579,1383580,1383581,1383582,
	1383583,1406147,658397,1173217,1175308,290280,707384,707385,707386,707387,546970,707388,546971,546960,546961,546962,546963,546956,546957,
	546958,546959,1184496,1195613,1195614,1195615,1195616,1195617,1195618,1195619,1195620,1195621,1195622,1195623,1195624,1195625,1195626,
	1195627,1195628,1195629,1195630,1195631,1195632,1195633,1195634,1195635,1195636,1195637,1195638,1195639,1195640,587432,546977,546978,
	546979,714567,546975,546976,546951,546952,546953,546954,546955,546983,546984,546985,546986,546987,546964,546965,546966,546980,546981,
	546982,546972,546973,546974,546967,546968,546969,1380573,1380574,1380575,1380576,1380577,1046049,934293,748172,1458701,863240,1260941,
	1406148,1325533,1195791,1241966,1331832,1331833,1331834,1331835,1331836,1331837,1331838,1331839,1331840,1331841,1331842,1331843,1331844,
	1331845,1331846,1331847,1331848,1241967,660656,660657,660658,660659,1184497,1241968,1566828,1379486,1379487,1282423,1282425,1262401,
	1337425,1268361,1458702,1458703,1458704,1241969,880562,1432852,1432853,1406149,392172,1241970,1133122,1383602,1406809,1406799,1383591,
	1383599,1383601,1406805,1406801,1406802,1383590,1406797,1406798,1383596,1406800,1383609,1383606,1406796,1316726,1221622,1221623,1241971,
	1241972,1183312,1312076,1498296,1241973,1241974,1260944,1282430,1458705,1183319,1183320,1175309,1175310,1175311,1175312,1291872,1291873,
	1486684,1486685,1486687,1486689,1486691,1486692,1486694,1486697,1336082,1336084,1336095,1336100,1336131,756495,1338688,1183318,1261195,
	1261196,1261197,1261198,1261199,1261200,1261201,1261202,1261203,1261204,1261205,1261206,1261207,1261208,1261209,1261210,1261211,1261212,
	1261230,1261231,1261232,1261233,1261234,1261235,1261236,1261237,1261238,1261239,1261240,1261241,1261242,1261243,1261244,1261245,1261246,
	1261247,1261248,1261249,1261250,1261251,1261252,1261253,1261254,1261255,1261256,1261257,1261258,1261259,1261260,1261261,1261262,1261263,
	756497,756496,64292,2025360,31699,12542,2169701,10258,118655,40058,340907,1803956,1511785,1352534,1395611,1395615,1395620,1341019,1608453,
	61673,64314,11083,1511871,11587,1980486,11292,378809,12814,1439707,186539,129003,147711,185935,185944,147712,12131,147679,147678,463676,
	11588,64285,64315,1511807,11029,28875,36427,11041,434309,59303,11080,1330524,1547495,651733,206160,688699,290314,334584,291175,1346530,
	985691,758863,1280920,234601,228407,321147,511429,511430,511431,511432,511433,321149,228415,388737,228404,627442,253433,253434,253435,
	228330,247149,353145,40012,11033,944645,1980490,12557,2697049,694009,1003835,992212,585082,585060,585083,585084,585085,585086,585063,
	585064,585087,585088,585066,585067,585068,585089,585090,585070,585091,585092,585072,585073,585093,585094,585095,585096,585097,585098,
	585099,585077,585078,585079,585080,585081,11642,1980491,11034,44419,84272,289366,289367,52276,1452514,1277649,44027,186540,1216927,
	186541,99000,11084,994672,687379,2065051,2065052,2065053,2065054,2065055,687380,2065044,2065045,2065046,2065047,2065048,2065049,2065050,
	687369,2065036,2065037,2065038,1859149,687370,687371,687373,687374,687375,687376,687377,1535290,1535291,68887,687340,687349,687351,687354,
	687355,687358,687364,687366,687367,687342,687345,687346,687347,687348,11590,862909,93678,1980494,1862826,64286,11591,10245,126794,10249,
	10255,11036,11277,164416,11082,449280,449279,449277,449278,11039,440266,749340,1049224,356663,356664,373193,11837,38804,132475,11089,617102,186538,64320) 
krakendf %>% filter(!(taxid %in% vhdb_human) & directreads > 50) %>% left_join(cancertypedf, by=c("sample"="sampleId")) %>% View()

######
# mapping rate
mapdf = read_tsv("flagstat.tsv", col_names=c("sample", "set", "reads"), col_types="cci") %>%
	spread(set, reads) %>%
	mutate(pct=adjusted/unadjusted) %>%
	inner_join(summarydf)
	# or if we don't double-count samples:
	#inner_join(summarydf %>% group_by(sample) %>% summarise(
	#	genus=paste(name_genus, collapse=" "),
	#	name=paste(name, collapse=" "),
	#	species=paste(name_species, collapse=" "),
	#	))
View(mapdf %>% filter(pct != 1))
ggplot(mapdf) +
	aes(x=pct) +
	geom_histogram() +
	facet_wrap(~genus)

mapdf %>%
	group_by(name_genus) %>%
	summarise(
		no_change=sum(pct==1),
		meanpct=mean(pct),
		maxpct=max(pct))
mapdf %>%
	group_by(name_genus) %>%
	filter(pct > 1) %>%
	summarise(
		meanpct=mean(pct),
		maxpct=max(pct)) %>%
	arrange(meanpct)

mapdf %>%
	group_by(name_genus) %>%
	filter(pct > 1) %>%
	summarise(
		meanpct=mean(pct),
		maxpct=max(pct)) %>%
	arrange(meanpct)
