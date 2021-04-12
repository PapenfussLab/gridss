source("libbenchmark.R")
library(rtracklayer)
library(assertthat)

flank_padding=100
flank_size=10000
# centromeric HOR array (chr8:44243868–46323885)
# chromosome 8 centromeric region (chr8:43600000–47200000)
chm13_actual_centro_gr = GRanges(seqnames="chr8", ranges=IRanges(start=43600000, end=47200000))
chm13_centro_gr = resize(chm13_actual_centro_gr, width=width(chm13_actual_centro_gr) * 2, fix="start") 

#samtools faidx chm13.draft_v1.0.fasta chr8 > ../../gridss2_manuscript/chm13.draft_v1.0.chr8.fasta
#gunzip -c chm13.draft_v1.0_plus38Y_repeatmasker.out.gz | grep "chr8 " | gzip -c > chm13.draft_v1.0_chr8_repeatmasker.out.gz
chm13_repeats_gr = import.repeatmasker.fa.out("publicdata/chm13.draft_v1.0_chr8_repeatmasker.out.gz")
liftover_chm13_to_hg38 = import.chain("publicdata/t2t-chm13-v1.0.hg38.over.chain")
liftover_hg38_to_hg19 = import.chain(system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain"))

chm13_ranges = GRanges(seqnames="chr8", ranges=ranges(liftover_chm13_to_hg38$chr8))
chm13_ranges = chm13_ranges[!overlapsAny(chm13_ranges, chm13_centro_gr)]
hg38_ranges = unlist(liftOver(chm13_ranges, liftover_chm13_to_hg38))
assertthat::assert_that(all(elementNROWS(liftOver(chm13_ranges, liftover_chm13_to_hg38))==1))
# Filter to only intervals that cleanly map to hg19
chm13_ranges = chm13_ranges[seqnames(hg38_ranges) == "chr8" & elementNROWS(liftOver(hg38_ranges, liftover_hg38_to_hg19)) == 1 & width(chm13_ranges) > 2 * flank_size + flank_padding]
hg19_ranges = unlist(liftOver(unlist(liftOver(chm13_ranges, liftover_chm13_to_hg38)), liftover_hg38_to_hg19))
# further restrict to full-length 1:1 mappings
chm13_ranges = chm13_ranges[width(chm13_ranges) == width(hg19_ranges)]
hg19_ranges = unlist(liftOver(unlist(liftOver(chm13_ranges, liftover_chm13_to_hg38)), liftover_hg38_to_hg19))
assertthat::are_equal(length(hg19_ranges), length(chm13_ranges))
assertthat::are_equal(width(hg19_ranges), width(chm13_ranges))
assertthat::assert_that(all(seqnames(hg19_ranges) == "chr8"))

chm13_anchor_position_gr = GRanges(seqnames=seqnames(chm13_ranges), ranges=IRanges(start=start(chm13_ranges) + flank_padding, width=flank_size))
chm13_flank_position_gr = resize(chm13_anchor_position_gr, width=1, fix="start")
chm13_break_position_gr = resize(chm13_anchor_position_gr, width=1, fix="end")
hg19_break_position_gr=unlist(liftOver(unlist(liftOver(chm13_break_position_gr, liftover_chm13_to_hg38)), liftover_hg38_to_hg19))
centr_stride = floor((width(chm13_centro_gr) - 2 * flank_size) / length(chm13_break_position_gr))
chm13_centro_position_gr = GRanges(
	seqnames="chr8",
	ranges=IRanges(
		start=seq(from=start(chm13_centro_gr) + flank_size, by=centr_stride, length.out=length(chm13_break_position_gr)),
		width=flank_size))
chm13_centro_break_position_gr = resize(chm13_centro_position_gr, 1, fix="start")
hg38_cbp_grl=liftOver(chm13_centro_break_position_gr, liftover_chm13_to_hg38)
chm13_centro_break_position_gr$hg38pos = NA_integer_
chm13_centro_break_position_gr$hg38chr = NA_character_
chm13_centro_break_position_gr$hg38pos[elementNROWS(hg38_cbp_grl) == 1] = start(unlist(hg38_cbp_grl[elementNROWS(hg38_cbp_grl) == 1]))
chm13_centro_break_position_gr$hg38chr[elementNROWS(hg38_cbp_grl) == 1] = seqnames(unlist(hg38_cbp_grl[elementNROWS(hg38_cbp_grl) == 1]))
cbp_has_38 = !is.na(chm13_centro_break_position_gr$hg38pos)
hg19_cbp_grl = liftOver(
	GRanges(
		seqnames=chm13_centro_break_position_gr$hg38chr[cbp_has_38],
		ranges(IRanges(width=1,start=chm13_centro_break_position_gr$hg38pos[cbp_has_38]))),
	liftover_hg38_to_hg19)
valid_hg19_cbp_grl = hg19_cbp_grl[elementNROWS(hg19_cbp_grl) == 1]
cbp_has_38_19 = elementNROWS(hg19_cbp_grl) == 1
chm13_centro_break_position_gr$hg19pos = NA_integer_
chm13_centro_break_position_gr$hg19chr = NA_character_
chm13_centro_break_position_gr$hg19pos[cbp_has_38][cbp_has_38_19] = start(unlist(valid_hg19_cbp_grl))
chm13_centro_break_position_gr$hg19chr[cbp_has_38][cbp_has_38_19] = seqnames(unlist(valid_hg19_cbp_grl))

chm13_ref = readDNAStringSet("publicdata/chm13.draft_v1.0.chr8.fasta")
truthdf = data.frame(
		chm_anchor_flank_pos = start(chm13_flank_position_gr),
		chm_anchor_break_pos =   end(chm13_break_position_gr),
		hg19_anchor_break_pos= start(hg19_break_position_gr),
		chm_centro_flank_pos =   end(chm13_centro_position_gr),
		chm_centro_break_pos = start(chm13_centro_position_gr),
		hg19_centro_break_pos= chm13_centro_break_position_gr$hg19pos,
		centro_repeat_type   = chm13_repeats_gr$repeatType [findOverlaps(chm13_centro_break_position_gr, chm13_repeats_gr, select="first")],
		centro_repeat_class  = chm13_repeats_gr$repeatClass[findOverlaps(chm13_centro_break_position_gr, chm13_repeats_gr, select="first")],
		anchor_repeat_type   = chm13_repeats_gr$repeatType [findOverlaps(       chm13_break_position_gr, chm13_repeats_gr, select="first")],
		anchor_repeat_class  = chm13_repeats_gr$repeatClass[findOverlaps(       chm13_break_position_gr, chm13_repeats_gr, select="first")],
		in_centro            = overlapsAny(chm13_centro_position_gr, chm13_actual_centro_gr)
			) %>%
	mutate(
		id=row_number(),
		id=paste0("run_", id, "_"))
if (!file.exists("publicdata/sim/all_runs.n.fa")) {
	simfa = DNAStringSet(paste0(as.character(chm13_ref[chm13_anchor_position_gr]), as.character(chm13_ref[chm13_centro_position_gr])))
	names(simfa) = truthdf$id
	writeXStringSet(simfa, "publicdata/sim/all_runs.t.fa")
	# in the normal we include the sequence on both sides of both breakends
	# so we have support for the ref allele in the normal
	simfa = DNAStringSet(c(
		as.character(chm13_ref[resize(chm13_anchor_position_gr, fix="start", width=width(chm13_anchor_position_gr)*2)]),
		as.character(chm13_ref[resize(chm13_centro_position_gr, fix="end", width=width(chm13_centro_position_gr)*2)])))
	names(simfa) = c(paste0("anchor_", truthdf$id), paste0("centro_", truthdf$id))
	writeXStringSet(simfa, "publicdata/sim/all_runs.n.fa")
}
# Now run sim.sh to generate all the results





