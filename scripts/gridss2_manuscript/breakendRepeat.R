source("libbenchmark.R")
library(VariantAnnotation)
# gunzip -c DO7358T.purple.sv.vcf.gz | grep "^#" > allbreakends.vcf
# gunzip -c *.gz | grep -E "gridss[0-9]+[fb]_[0-9]+b" >> allbreakends.vcf

bevcf = readVcf(paste0(privatedatadir, "/pcawg/allbreakends.vcf"))
bevcf = bevcf[rowRanges(bevcf)$FILTER == "PASS"]
lnx_svs = read_csv(paste0(privatedatadir, "/hartwig/LNX_SVS.csv")) %>%
	filter(ChrEnd == 0 & QualScore > 0)

repeatdf = bind_rows(
	data.frame(
		seq=as.character(rowRanges(bevcf)$ALT),
		length=nchar(as.character(rowRanges(bevcf)$ALT)) - 1,
		rc=info(bevcf)$INSRMRC,
		rt=info(bevcf)$INSRMRT,
		cohort="PCAWG"),
	data.frame(
		seq=lnx_svs$InsertSeq,
		length=nchar(lnx_svs$InsertSeq),
		rc=lnx_svs$RepeatClass,
		rt=lnx_svs$RepeatType,
		cohort="Hartwig"))

rt_lookup = c(
	"L1HS"="L1HS",
	"(T)n"="poly A",
	"(A)n"="poly A",
	"ALR/Alpha"="ALR/Alpha",
	"HSATII"="Satellite",
	"(CATTC)n"="Satellite", # pericentromeric
	"(GAATG)n"="Satellite",
	"(CCCTAA)n"="Satellite", # telomeric
	"(TTAGGG)n"="Satellite")
rc_lookup = c(
	"SINE/Alu"="SINE/Alu")
rcc_lookup = c(
	"Low_complexity"="Simple/Low complexity",
	"Simple_repeat"="Simple/Low complexity",
	"Satellite"="Satellite",
	"LINE"="LINE",
	"SINE"="SINE",
	"LTR"="Other",
	"Other"="Other",
	"DNA"="Other",
	"rRNA"="Other",
	"scRNA"="Other",
	"tRNA"="Other",
	"snRNA"="Other",
	"LTR?"="Other",
	"RC?"="Other",
	"SINE?"="SINE",
	"srpRNA"="Other",
	"RNA"="Other")
simple_repeat_label = function(rc, rt) {
	x = ifelse(rt %in% names(rt_lookup), rt_lookup[rt],
				 ifelse(rc %in% names(rc_lookup), rc_lookup[rc],
				 			 rcc_lookup[str_replace(rc, "[/].*$", "")]))
	x[is.na(x)] = "No Repeat"
	return(factor(x, levels=c(
		"No Repeat", "Other",
		"LINE", "L1HS",
		"SINE", "SINE/Alu",
		"Satellite", "ALR/Alpha",
		"Simple/Low complexity", "poly A")))
}

ggplot(repeatdf) +
	aes(x=length, fill=simple_repeat_label(rc, rt)) +
	geom_histogram(bins=30) +
	scale_x_continuous(limits=c(0, 800), expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0)) +
	scale_fill_manual(values=c(
		"#000000", "#666666", #1b9e77
		"#d95f02", "#853900", #d95f02
		"#66a61e", "#29420c", #66a61e
		"#8b85d6", "#534f80", ##7570b3
		"#911756", "#3d0a25" #e7298a 
		)) +
	facet_wrap( ~ cohort, scales="free") +
	labs(y="Single breakend variant calls", x="Breakend assembly length", fill="RepeatMasker Annotation")
figsave("breakend_repeat_annotation_by_length", height=4, width=10)


# grep clusterReason DO7196T.linx.svs.tsv > pcawg_lnx_svs.tsv
# grep -v clusterReason *.linx.svs.tsv >> pcawg_lnx_svs.tsv
lnx_svs %>% filter(simple_repeat_label(RepeatClass, RepeatType)== "poly A") %>% group_by(ResolvedType) %>% summarise(n=n()) %>% ungroup() %>% mutate(pct=n/sum(n))
hartwig_line_sv = read_csv(paste0(privatedatadir, "/hartwig/LNX_SVS.csv")) %>% filter(ResolvedType=="LINE") %>% summarise(n=n())

pcawg_linx_svs = read_tsv(paste0(privatedatadir, "/pcawg/pcawg_lnx_svs.tsv"))
pcawg_linx_svs %>% filter(simple_repeat_label(RepeatClass, RepeatType)== "poly A") %>% group_by(ResolvedType) %>% summarise(n=n()) %>% ungroup() %>% mutate(pct=n/sum(n))








