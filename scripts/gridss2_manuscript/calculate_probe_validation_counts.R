rawprobeResult = read_csv(paste0(figdir, "supptable_probe_validation_results.csv"))
probeResult = rawprobeResult %>% filter(is.na(exclusion))

required_overlap = 25
probeResult = probeResult %>%
	mutate(
		contig=uniqueId,
		contig_break_pos=ifelse(is.na(probeResult$endChromosome), ifelse(startOrientation == 1, 2000 - str_length(insertSequence), str_length(insertSequence)), 1000),
		contig_break_lower_pos=pmax(1, ((contig_break_pos - str_length(startHomologySequence)) %na% contig_break_pos) - required_overlap),
		contig_break_upper_pos=pmin(2000, ((contig_break_pos + str_length(endHomologySequence)) %na% contig_break_pos) + required_overlap))
required_overlap_gr = with(probeResult, GRanges(
	seqnames=contig,
	ranges=IRanges(start=contig_break_lower_pos, end=contig_break_upper_pos)))
export(con="probe_required_overlap.bed", required_overlap_gr)

# subset to probe contigs
# for f in HMF*A.bam ; do samtools view -@ $(nproc) -b $f $(grep -vE "^([0-9])|X|Y|MT" probe_ref_genome/probes.fasta.fai | cut -f 1) > $f.probecontigs.bam ; done
#for f in HMF*.probecontigs.bam ; do echo "bedtools intersect -f 1 -c -a probe_required_overlap.bed -b $f > $f.probeoverlaps.bed & "; done
#testing: bedtools intersect -f 1 -c -a <(head -10 probe_required_overlap.bed) -b HMF000214A.bam.probecontigs.bam  HMF000379A.bam.probecontigs.bam 