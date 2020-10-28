library(tidyverse)
library(VariantAnnotation)

krakendf = read_tsv(
	"../protecteddata/hmf/viral.out",
	col_name=c("file", "pct", "treereads", "directreads", "level", "taxid", "name"),
	col_types="cniicic")

extracteddf = read_tsv(
	"../protecteddata/hmf/extracted.out",
	col_name=c("file", "pct", "treereads", "directreads", "level", "taxid", "name"),
	col_types="cniicic")

wgsdf = read_tsv(
	"../protecteddata/hmf/wgs_metrics.out",
	col_names=c("file","GENOME_TERRITORY","MEAN_COVERAGE","SD_COVERAGE","MEDIAN_COVERAGE","MAD_COVERAGE","PCT_EXC_ADAPTER","PCT_EXC_MAPQ","PCT_EXC_DUPE","PCT_EXC_UNPAIRED","PCT_EXC_BASEQ","PCT_EXC_OVERLAP","PCT_EXC_CAPPED","PCT_EXC_TOTAL","PCT_1X","PCT_5X","PCT_10X","PCT_15X","PCT_20X","PCT_25X","PCT_30X","PCT_40X","PCT_50X","PCT_60X","PCT_70X","PCT_80X","PCT_90X","PCT_100X","HET_SNP_SENSITIVITY","HET_SNP_Q"),
	col_types="cddddddddddddddddddddddddddddd")

ggplot(wgsdf) +
	aes(x=PCT_1X) +
	geom_histogram() +
	labs(title="Coverage of extracted virus")

mergedvcf = readVcf("../protecteddata/hmf/merged.vcf")
rownames(mergedvcf) = NULL
alt = str_match(as.character(rowRanges(mergedvcf)$ALT), "([^\\[\\]]*)([\\[\\]])adjusted_kraken_taxid_([0-9]+)_(.*):([0-9]+)[\\[\\]]([^\\[\\]]*)")
colnames(alt) = c("ALT", "leftins", "ori", "taxid", "virus_chr", "virus_pos", "rightins")
sitedf = as.data.frame(alt) %>%
	mutate(
		host_chr=as.character(seqnames(mergedvcf)),
		host_pos=start(mergedvcf),
		QUAL=rowRanges(mergedvcf)$QUAL,
		FILTER=rowRanges(mergedvcf)$FILTER) %>%
	bind_cols(as.data.frame(info(mergedvcf)))
