library(StructuralVariantAnnotation)
for (file in list.files(paste0(datadir, "/gridss1"), full.names=TRUE, recursive=TRUE, pattern=".*.gridss.vcf")) {
	g1vcf = VariantAnnotation::readVcf(file)
	g1normalSupport = geno(g1vcf)$SR[,"normal"] + geno(g1vcf)$RP[,"normal"]
	g1tumour_support = geno(g1vcf)$SR[,"tumour"] + geno(g1vcf)$RP[,"tumour"]
	g1somaticvcf = g1vcf[g1normalSupport == 0]
	writeVcf(filename = paste0(file, ".somatic.vcf"), g1somaticvcf)
}


g1bpgr = breakpointRanges(g1somaticvcf)
g1hits = findBreakpointOverlaps(g1bpgr, colo829truthgr)
g1bpgr$tp = FALSE
g1bpgr$tp[queryHits(g1hits)] = TRUE

g1somaticpassvcf = g1somaticvcf[rowRanges(g1somaticvcf)$FILTER %in% c(".")]

ggplot(
	data.frame(g1bpgr) %>%
		filter(FILTER==".")) +
	aes(x=QUAL, fill=tp) +
	geom_histogram() +
	scale_x_log10() +
	labs(title="QUAL distribution of somatic passing variants")

ggplot(
	data.frame(
		normal_reads = geno(g1somaticpassvcf)$SR[,"normal"] + geno(g1somaticpassvcf)$RP[,"normal"],
		tumour_reads = geno(g1somaticpassvcf)$SR[,"tumour"] + geno(g1somaticpassvcf)$RP[,"tumour"],
		tp = names(g1somaticpassvcf) %in% g1bpgr[g1bpgr$tp]$sourceId
		)) +
	aes(x=tumour_reads, y=normal_reads, fill=tp) +
	geom_point() +
	scale_x_log10() +
	scale_y_log10()
	

geno_to_dataframe = function(vcf, sample) {
	g = geno(vcf)
	outdf = data.frame(GT=g$GT[,sample])
	for (field in names(g)) {
		outdf[field] = g[field][,sample]
	}
	return(outdf)
}
z = geno_to_dataframe(g1somaticvcf, "normal")

lapply(geno(g1somaticvcf[rowRanges(g1somaticvcf)$QUAL > 5000]), 
