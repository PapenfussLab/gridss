source("libbenchmark.R")
library("tidyverse")
rawprobeResult = read_csv(paste0(figdir, "supptable_probe_validation_results.raw.csv"))
probeResult = rawprobeResult

required_overlap = 25
probeResult = probeResult %>%
	replace_na(list(
		insertSequence="",
		startHomologySequence="",
		endHomologySequence="")) %>%
	mutate(
		isbe=is.na(probeResult$endChromosome),
		insBases=str_length(insertSequence),
		contig=uniqueId,
		contig_break_pos=ifelse(isbe, ifelse(startOrientation == 1, 2000 - insBases, insBases), 1000),
		contig_break_lower_pos=pmax(1, contig_break_pos - str_length(startHomologySequence) - required_overlap - ifelse(isbe, 0, ceiling(insBases / 2))),
		contig_break_upper_pos=pmin(2000, contig_break_pos + str_length(endHomologySequence) + required_overlap + ifelse(isbe, 0, ceiling(insBases / 2))))
required_overlap_gr = with(probeResult, GRanges(
	seqnames=contig,
	ranges=IRanges(start=contig_break_lower_pos, end=contig_break_upper_pos)))
export(con="probe_required_overlap.bed", required_overlap_gr)

# subset to probe contigs
# for f in HMF*A.bam ; do samtools view -@ $(nproc) -b $f $(grep -vE "^([0-9])|X|Y|MT" probe_ref_genome/probes.fasta.fai | cut -f 1) > $f.probecontigs.bam ; done
#for f in HMF*.probecontigs.bam ; do echo "bedtools intersect -f 1 -c -a probe_required_overlap.bed -b $f > $f.probeoverlaps.bed & "; done
#testing: bedtools intersect -f 1 -c -a <(head -10 probe_required_overlap.bed) -b HMF000214A.bam.probecontigs.bam  HMF000379A.bam.probecontigs.bam
alt_hits_df = bind_rows(lapply(unique(probeResult$sampleId), function(sampleId) {
	read_tsv(paste0("protecteddata/validation/", sampleId, ".bam.probecontigs.bam.probeoverlaps.bed"),
					 col_names = c(
					 	"uniqueId",
					 	# (in bed coordinates)
					 	"contig_break_lower_pos",
					 	"contig_break_upper_pos",
					 	"unused1",
					 	"unused2",
					 	"unused3",
					 	"hits"),
					 col_types = "ciiccci") %>%
		dplyr::select(uniqueId, hits, contig_break_lower_pos, contig_break_upper_pos) %>%
		mutate(
			sampleId=sampleId,
			alt_required_span_length=contig_break_upper_pos - contig_break_lower_pos + 1)
	})) %>%
	left_join(dplyr::select(probeResult, sampleId, uniqueId) %>% mutate(is_source=TRUE), by=c("sampleId", "uniqueId")) %>%
	replace_na(list(is_source=FALSE)) %>%
	group_by(uniqueId, alt_required_span_length) %>%
	dplyr::summarise(
		sourceAltDepth=sum(ifelse(is_source, hits, 0)),
		otherAltDepthMax=max(ifelse(is_source, 0, hits)),
		otherAltDepth=sum(ifelse(is_source, 0, hits)),
		otherAltCount=sum(ifelse(is_source | hits == 0, 0, 1))
		) %>%
	mutate(
		p=ppois(sourceAltDepth,otherAltDepthMax,FALSE,FALSE),
		supported=p<0.001 & otherAltDepth<40 & sourceAltDepth>1) %>%
	ungroup()


ggplot(bind_rows(
		alt_hits_df %>% mutate(type="overlap"),
		probeResult %>% mutate(type="depth"))) +
	aes(x=sourceAltDepth + 1, y=otherAltDepth + 1) +
	facet_wrap(~type) +
	scale_x_log10() +
	scale_y_log10() +
	geom_jitter(size=0.1)
ggplot(bind_rows(
	alt_hits_df %>% mutate(type="overlap"),
	probeResult %>% mutate(type="depth"))) +
	aes(x=width, fill=sourceAltDepth > 0) +
	geom_histogram(bins=100) +
	facet_wrap(~type)

write_csv(
	rawprobeResult %>%
		dplyr::select(
			-sourceAltDepth,
			-otherAltDepthMax,
			-otherAltDepth,
			-otherAltCount,
			-supported,
			-p) %>%
		left_join(alt_hits_df, by="uniqueId"),
	paste0(figdir, "supptable_probe_validation_results.csv"))













