source("libbenchmark.R")
lnx_svs = read_csv(paste0(privatedatadir, "/hartwig/LNX_SVS.csv"))
lnx_links = read_csv(paste0(privatedatadir, "/hartwig/LNX_LINKS.csv"))
lnx_clusters = read_csv(paste0(privatedatadir, "hartwig/LNX_CLUSTERS.csv"))
pcawg_linx_svs = read_tsv(paste0(privatedatadir, "/pcawg/pcawg.linx.vis_sv_data.tsv"))
pcawg_linx_svs$Id = pcawg_linx_svs$SvId
pcawg_lnx_vis_svs = read_tsv(paste0(privatedatadir, "pcawg/all_pcawg_linx_vis_sv_data.tsv"))
pcawg_lnx_clusters = read_tsv(paste0(privatedatadir, "pcawg/all_pcawg_linx.clusters.tsv"))
pcawg_lnx_links = read_tsv(paste0(privatedatadir, "pcawg/all_pcawg_linx.links.tsv"))

#### Figure candidates:
hartwig_clusters_with_many_breakends = lnx_svs %>%
	group_by(SampleId, ClusterId) %>%
	replace_na(list(RepeatClass="", RepeatType="")) %>%
	mutate(
		ncentromeric=sum(RepeatClass == "Satellite/centr" | RepeatType %in% c("(CATTC)n", "(GAATG)n", "HSATII")),
		nbreakends=sum(Type=="SGL"),
		nnoncentromericbreakends=nbreakends-ncentromeric,
		nsvs=n()) %>%
	filter(ResolvedType=="COMPLEX" & ncentromeric > 10)
hartwig_clusters_with_many_breakends %>%
	dplyr::select(SampleId, ClusterId, ncentromeric, nbreakends, nnoncentromericbreakends, nsvs) %>%
	distinct() %>%
	inner_join(lnx_clusters, by=c("SampleId", "ClusterId")) %>%
	arrange(FullyChained, desc(ncentromeric)) %>%
	View()


#### Figure candidates:
lnx_svs %>% 
	dplyr::select(SampleId, Id, ChainId, ChainCount, ChainIndex) %>%
	arrange(SampleId, ChainId, ChainIndex) %>%
	filter(!is.na(ChainIndex)) %>%
	View()

#######
# Candidate fully assembled rearrangements
pcawg_lnx_links %>%
	group_by(sample, clusterId) %>%
	filter(length(unique(chainId)) == 1) %>%
	group_by(sample, clusterId, chainId) %>%
	filter(all(assembled==TRUE)) %>%
	filter(chainCount == 4) %>%
	inner_join(pcawg_lnx_clusters, by=c("sample", "clusterId")) %>%
	filter(resolvedType == "SIMPLE") %>%
	View()




#### Figure candidates:
lnx_links %>%
	dplyr::select(SampleId, ClusterId, ClusterCount, ChainId, ChainCount, ChainConsistent) %>%
ggplot() +
	aes(x=ChainCount) +
	geom_histogram(bins=51) +
	scale_x_continuous(limits=c(0, 50))


lnx_links %>%
	dplyr::select(SampleId, ClusterId, ChainId, ChainCount) %>%
	distinct() %>%
	mutate(chainedLength=ChainCount+1) %>%
	group_by(chainedLength) %>%
	summarise(svs=max(chainedLength)*n()) %>%
	ungroup() %>%
	arrange(desc(chainedLength)) %>%
	mutate(svsAtLeast=cumsum(svs)) %>%
	mutate(pctSvsChainedAtLengthOrMore=svsAtLeast/nrow(lnx_svs)) %>%
	ggplot() +
	aes(x=chainedLength, y=pctSvsChainedAtLengthOrMore) +
	geom_line() +
	scale_x_continuous(limits=c(3,30), expand=c(0,0))

# Inferred FDR

####
# Paper figures
# % of links from assembly
sum(lnx_links$IsAssembled)/nrow(lnx_links)
# % of links from assembly
sum(lnx_links$IsAssembled)/nrow(lnx_svs)

summarise_when_grouped_into_chains = function(df) {
	df %>%
		mutate(variantsInChain=n()) %>%
		group_by(variantsInChain) %>%
		summarise(svs=max(variantsInChain)*n()) %>%
		ungroup() %>%
		arrange(desc(variantsInChain)) %>%
		mutate(svsAtLeast=cumsum(svs)) %>%
		mutate(pctSvsChainedAtLengthOrMore=svsAtLeast/sum(svs))
}
downsample_to_fnr = function(svdf, estimated_fnr, target_fnr) {
	total_svs = nrow(svdf) * 1/(1-estimated_fnr)
	target_svs = nrow(svdf) * (1-target_fnr)
	# minimal data frame with linearised chaining
	svdf = svdf %>% dplyr::select(SampleId, ClusterId, ChainId) %>%
		inner_join(svdf %>%
			dplyr::select(SampleId, ClusterId, ChainId) %>%
			distinct() %>%
			mutate(uid=10000*row_number())) %>%
		dplyr::select(uid) %>%
		group_by(uid) %>%
		mutate(ordinal=row_number())
	# recalculate chaining for target FNR
	subsetdf = svdf %>%
		slice_sample(n=target_svs) %>%
		group_by(uid) %>%
		arrange(ordinal) %>%
		mutate(
			starts_new_chain=!is.na(lag(ordinal)) & lag(ordinal) != ordinal - 1,
			ordinal_of_previous_chain_start=lastTrueOrdinal(starts_new_chain, 0),
			new_uid=uid + ordinal_of_previous_chain_start)
	return(subsetdf %>%
		group_by(new_uid) %>%
		summarise_when_grouped_into_chains() %>%
			mutate(fnr=target_fnr))
}
#' second implmentation that should be much faster
downsample_to_fnr2 = function(svdf, estimated_fnr, target_fnrs) {
	total_svs = nrow(svdf) * 1/(1-estimated_fnr)
	chain_lengths = svdf %>%
		dplyr::select(SampleId, ClusterId, ChainId) %>%
		group_by(SampleId, ClusterId, ChainId) %>%
		summarise(n=n()) %>%
		pull(n) 
	uids = rep(seq_along(chain_lengths), times=chain_lengths)
	ordinals = distanceToLastTrueValue(is.na(lag(uids)) | uids != lag(uids))
	result = data.frame()
	for (target_fnr in target_fnrs) {
		target_svs = nrow(svdf) * (1-target_fnr)
		subset_ordinals = sample.int(length(uids), target_svs)
		subset_ordinals = sort(subset_ordinals)
		new_uids = uids[subset_ordinals]
		new_ordinals = ordinals[subset_ordinals]
		starts_chain = is.na(lag(new_uids)) |
			lag(new_uids != new_uids) | lag(new_ordinals) != new_ordinals -1
		replaced_uids = rep(-1, length(new_uids))
		replaced_uids[starts_chain] = seq_along(replaced_uids[starts_chain])
		replaced_uids=cummax(replaced_uids)
		new_chain_lengths = as.numeric(table(replaced_uids))
		result = bind_rows(
			result,
			data.frame(variantsInChain=new_chain_lengths) %>%
				group_by(variantsInChain) %>%
				summarise(svs=max(variantsInChain)*n()) %>%
				ungroup() %>%
				arrange(desc(variantsInChain)) %>%
				mutate(
					svsAtLeast=cumsum(svs),
					pctSvsChainedAtLengthOrMore=svsAtLeast/sum(svs),
					fnr=target_fnr)
			)
	}
	return(result)
}

hartwig_fdr_downsample_df = downsample_to_fnr2(
		lnx_svs,
		estimated_fnr=0.025,
		target_fnrs=c(0.025, 0.05, 0.10, 0.112, 0.15, 0.20, 0.30)) %>%
	mutate(cohort="Hartwig")

pcawg_fdr_df = downsample_to_fnr2(pcawg_lnx_vis_svs, estimated_fnr=0.112, target_fnrs=c(0.112))


fdr_df = bind_rows(
		hartwig_fdr_downsample_df %>%
			mutate(label=ifelse(fnr==0.025, "Hartwig cohort", paste0("Simulated FNR: ", fnr * 100, "%"))),
		pcawg_fdr_df %>%
			mutate(label="PCAWG cohort")) %>%
	mutate(isSimulated=!str_detect(label, "cohort"))
plot_fdr_df = fdr_df %>% filter(label != "Simulated FNR: 11.2%")
ggplot(plot_fdr_df) +
	aes(
		x=variantsInChain,
		y=pctSvsChainedAtLengthOrMore,
		colour=label) +
	geom_line(data=plot_fdr_df %>% filter(!isSimulated), size=1.5) +
	geom_line(data=plot_fdr_df %>% filter(isSimulated)) +
	geom_point(data=plot_fdr_df %>% filter(!isSimulated), size=3) +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0), labels = scales::percent) +
	coord_cartesian(xlim=c(-0.00000001, 40.000001), ylim=c(0,0.1500001)) +
	scale_color_brewer(palette="Dark2") +
	scale_color_manual(values=c("#1b9e77", "#7570b3", "#636363", "#969696", "#bdbdbd", "#dddddd", "#252525")) +
	labs(y="Portion of called SVs", x="Chain length")
	#labs(y="Portion of called SVs", x="Chained with at least this many other SVs")
	#labs(y="SVs in chains at least length", x="Chain length")
figsave("figure5_chain_length_fnr", width=6, height=4.5)


#####
# Paper results
lnx_links %>%
	group_by(LinkReason) %>%
	summarise(count=n()) %>%
	ungroup() %>%
	mutate(
		total=sum(count),
		pct=count/total)
# grep clusterId DO7196T.linx.links.tsv > all_pcawg_linx.links.tsv
# cat *.linx.links.tsv | grep -v clusterId >> all_pcawg_linx.links.tsv
pcawg_lnx_links = read_tsv(paste0(privatedatadir, "pcawg/all_pcawg_linx.links.tsv"))
pcawg_lnx_links %>%
	summarise(
		assembled=sum(assembled),
		total=n(),
		pct=assembled/total)

lnx_clusters %>% filter(ResolvedType %in% c("DOUBLE_MINUTE", "COMPLEX")) %>%
	group_by(FullyChained) %>%
	summarise(count=n()) %>%
	ungroup() %>%
	mutate(
		total=sum(count),
		pct=count/total)

fdr_df %>% filter(variantsInChain==20)


# Portion of chains fully-resolved by assembly
pcawg_lnx_links %>%
	group_by(sample, clusterId, chainId) %>%
	summarise(links=n(), assembled=sum(assembled)) %>%
	group_by(links) %>%
	summarise(
		count=n(),
		fully_assembled=sum(links==assembled),
		pct_fully_assembled=100*fully_assembled/count)
















