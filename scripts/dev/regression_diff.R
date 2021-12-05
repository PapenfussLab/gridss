library(StructuralVariantAnnotation)
library(argparser)
library(tidyverse)
baseline_filename = "v2.9.3/colo829/somatic.vcf.bgz"
latest_filename = "vdev-c45021b5/colo829_10_externalaligner/somatic.vcf.bgz"

ann_hit_matches = function(hitdf, queryGr, queryVcf, subjectGr, subjectVcf) {
	hitdf = hitdf |>
		as.data.frame() |>
		mutate(
			queryName=queryGr$sourceId[queryHits],
			subjectName=subjectGr$sourceId[subjectHits]) |>
		mutate(
			matchQUAL=rowRanges(vcf_b)$QUAL[queryHits] == rowRanges(vcf_l)$QUAL[subjectHits],
			matchFILTER=rowRanges(vcf_b)$FILTER[queryHits] == rowRanges(vcf_l)$FILTER[subjectHits],
			matchstart=start(vcf_b)[queryHits] == start(vcf_l)[subjectHits],
			matchend=end(vcf_b)[queryHits] == end(vcf_l)[subjectHits],
			matchpos=matchstart & matchend)
	iq = as.data.frame(info(queryVcf))
	is = as.data.frame(info(subjectVcf))
	#"CIPOS","CIRPOS","IHOMPOS","HOMLEN","HOMSEQ",
	for (field in c("AS","ASC","ASQ","ASRP","ASSR","BA","BANRP","BANRPQ","BANSR",
									"BANSRQ","BAQ","BASRP","BASSR","BMQ","BMQN","BMQX","BQ","BSC",
									"BSCQ","BUM","BUMQ","BVF","CAS","CASQ","CQ",
									"IC","IMPRECISE","INSRMP","INSRMRC","INSRMRO","INSRMRT","IQ",
									"MQ","MQN","MQX","RAS","RASQ","REF","REFPAIR","RP","RPQ",
									"SB","SC","SR","SRQ","VF")) {
		qval = iq[[field]][hitdf$queryHits]
		sval = is[[field]][hitdf$subjectHits]
		if (xor(is.null(qval),is.null(sval))) {
			hitdf[[paste0("match", field)]] = TRUE
		} else {
			hitdf[[paste0("match", field)]] = ifelse(is.na(qval) | is.na(sval), is.na(qval) & is.na(sval), qval == sval)
		}
	}
	match_mat = as.matrix(hitdf |> dplyr::select(starts_with("match")))
	name_mat = matrix(ncol=length(colnames(match_mat)), rep(str_remove(colnames(match_mat), "^match"), each=nrow(match_mat)))
	hitdf$mismatches = apply(ifelse(match_mat, "", name_mat), 1, function(x) paste0(x, collapse=",")) |>
		str_replace_all("[,]{2,}", ",") |>
		str_replace_all("^,", "") |>
		str_replace_all(",$", "")
	hitdf |> dplyr::select(-starts_with("match"), -queryHits, -subjectHits)
}
write_subset = function(vcf, record_names, output_prefix, output) {
	record_names = unique(record_names)
	write(record_names, paste0(prefix, output, ".txt"))
	writeVcf(vcf[record_names], paste0(prefix, output, ".vcf"))
}
dir.create("diff")
for (need_PASS in c(TRUE, FALSE)) {
	vcf_b = readVcf(baseline_filename)
	vcf_l = readVcf(latest_filename)
	prefix="diff/all_"
	if (need_PASS) {
		vcf_b = vcf_b[rowRanges(vcf_b)$FILTER %in% c("", ".", "PASS")]
		vcf_l = vcf_l[rowRanges(vcf_b)$FILTER %in% c("", ".", "PASS")]
		prefix="diff/pass_"
	}
	bpgr_b = breakpointRanges(vcf_b)
	bpgr_l = breakpointRanges(vcf_l)
	begr_b = breakendRanges(vcf_b)
	begr_l = breakendRanges(vcf_l)
	hits = bind_rows(
		findBreakpointOverlaps(bpgr_b, bpgr_l) |> ann_hit_matches(bpgr_b, vcf_b, bpgr_l, vcf_l),
		findOverlaps(begr_b, begr_l) |> ann_hit_matches(begr_b, vcf_b, begr_l, vcf_l))
	write_subset(vcf_b, names(vcf_b)[!names(vcf_b) %in% hits$queryName], prefix, "baseline_only")
	write_subset(vcf_l, names(vcf_l)[!names(vcf_l) %in% hits$subjectName], prefix, "latest_only")
	write_subset(vcf_b, hits |> filter(mismatches=="") |> pull(queryName), prefix, "baseline_perfect_match")
	write_subset(vcf_l, hits |> filter(mismatches=="") |> pull(subjectName), prefix, "latest_perfect_match")
	write_tsv(hits |> filter(mismatches!="") |> dplyr::select(baseline=queryName, latest=subjectName) |> as.data.frame(), paste0(prefix, "mismatched_info_columns.tsv"))
}
