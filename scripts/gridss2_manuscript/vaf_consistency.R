source("libbenchmark.R")
library(openxlsx)
library(ggExtra)
isolation_distance = 3000
lnx_svs = read_csv(paste0(privatedatadir, "/hartwig/LNX_SVS.csv"))
svlong = bind_rows(
	lnx_svs %>%
		mutate(isBreakend=ChrEnd==0) %>%
		dplyr::select(sampleId=SampleId, chr=ChrStart, pos=PosStart, cn=CNStart, cn_delta=CNChgStart, sv_ploidy=Ploidy, isBreakend),
	lnx_svs %>%
		mutate(isBreakend=ChrEnd==0) %>%
		filter(!isBreakend) %>%
		dplyr::select(sampleId=SampleId, chr=ChrEnd, pos=PosEnd, cn=CNEnd, cn_delta=CNChgEnd, sv_ploidy=Ploidy, isBreakend)) %>%
	group_by(sampleId, chr) %>%
	arrange(pos) %>%
	mutate(
		isolated_l = is.na(lead(pos)) | abs(pos - lead(pos)) > isolation_distance,
		isolated_r = is.na( lag(pos)) | abs(pos -  lag(pos)) > isolation_distance,
		isolated=isolated_l & isolated_r) %>%
	ungroup()
isolated_svs = svlong %>% filter(isolated)

isolated_non_inferred_svs = isolated_svs %>% filter(cn_delta != sv_ploidy)

cor(abs(isolated_svs$cn_delta), isolated_svs$sv_ploidy)
ggplot(isolated_non_inferred_svs %>% slice_sample(n=25000)) +
	aes(x=abs(cn_delta), y=sv_ploidy) +
	geom_jitter(alpha=0.12, shape=16) +
	geom_abline(slope=1, intercept=0) +
	scale_x_continuous(limits=c(0, 3.5), expand=c(0,0)) +
	scale_y_continuous(limits=c(0, 3.5), expand=c(0,0)) +
	labs(y="Expected change in copy number\ninfered from variant allele fraction", x="Actual change in copy number")
figsave("estimated_fdr_hartwig", width=4, height=4)

# paper figures
fdr_df = svlong %>%
	mutate(
		likely_fp=cn_delta<0.1 & sv_ploidy > 0.25,
		subclonal=sv_ploidy<0.75)

estimated_isolated_fdr = fdr_df %>%
	filter(isolated) %>%
	group_by(isBreakend, subclonal) %>%
	summarise(fdr=sum(likely_fp)/n())

fdr_df %>% inner_join(estimated_isolated_fdr, by=c("isBreakend", "subclonal")) %>%
	summarise(
		fp=sum(likely_fp),
		total=n(),
		fdr=fp/total) %>% pull(fdr)

isolated_non_inferred_svs %>%
	mutate() %>%
	group_by(likely_fp) %>%
	summarise(count=n()) %>%
	ungroup() %>%
	mutate(pct=count/sum(count))

isolated_non_inferred_svs %>%
	mutate(likely_fp=cn_delta<0.1 & sv_ploidy > 0.25) %>%
	filter(likely_fp & sv_ploidy < 0.75) %>%
	summarise(count=n())

