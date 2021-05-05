setwd("../")
source("libgridss.R")
setwd("gridss2_manuscript")
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
hg19_repeats_gr = import.repeatmasker.fa.out("publicdata/hg19.fa.out.gz")
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
		runid=row_number(),
		id=paste0("run_", runid, "_"))
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
# Then run sim_aggregate_results.sh
process_caller = function(caller) {
	vcf = VariantAnnotation::readVcf(paste0("publicdata/sim/merged_", caller, "_sim.vcf"))
	info(vcf)$RUNID = as.integer(str_extract(info(vcf)$RUNID, "[0-9]+"))
	row.names(vcf) = paste0("run", info(vcf)$RUNID, "_", row.names(vcf))
	if ("MATEID" %in% colnames(info(vcf))) {
		info(vcf)$MATEID = elementExtract(info(vcf)$MATEID)
		info(vcf)$MATEID[!is.na(info(vcf)$MATEID)] = paste0("run", info(vcf)$RUNID[!is.na(info(vcf)$MATEID)], "_", info(vcf)$MATEID[!is.na(info(vcf)$MATEID)])
	}
	if (caller == "novobreak") {
		# Columns aren't even in the correct order
		VariantAnnotation::fixed(vcf)$REF = DNAStringSet(rep("N", length(vcf)))
		names(vcf) = paste0("placeholder", seq(length(vcf)))
	}
	write(paste(caller, "filters:", unique(rowRanges(vcf)$FILTER)), stderr())
	vcf = vcf[VariantAnnotation::fixed(vcf)$FILTER %in% c("PASS", "minRead")]
	bpgr = breakpointRanges(vcf) #, inferMissingBreakends=TRUE)
	begr = breakendRanges(vcf)
	gr = c(bpgr, begr) 
	mcols(gr) = as.data.frame(mcols(gr)) %>%
		mutate(runid=info(vcf[gr$sourceId])$RUNID) %>%
		inner_join(truthdf, by="runid")
	# TODO: add hg19 repeat annotation
	gr$distance_to_anchor_break = abs(start(gr) - gr$hg19_anchor_break_pos)
	gr$distance_to_centro_break = abs(start(gr) - gr$hg19_centro_break_pos)
	#& strand(gr) == "+" no strand matching requirement as the DNA segment could be inverted in chm13
	gr$is_anchor_hit = as.logical(seqnames(gr) == "chr8") & gr$distance_to_anchor_break < 100
	gr$is_centro_hit = as.logical(seqnames(gr) == "chr8") & gr$distance_to_centro_break < 100
	begr = gr[is.na(gr$partner)]
	bpgr = gr[!is.na(gr$partner)]
	bpgr$treat_as_anchor = ifelse(bpgr$is_anchor_hit == partner(bpgr)$is_anchor_hit, bpgr$distance_to_anchor_break < partner(bpgr)$distance_to_anchor_break, bpgr$is_anchor_hit)
	bpanchorgr = bpgr[bpgr$treat_as_anchor]
	bpcentrogr = partner(bpgr)[bpgr$treat_as_anchor]
	# bealn repeatmasker lookup assumes all single breakends have a bealn value
	assert_that(all(elementNROWS(info(vcf[begr$sourceId])$BEALN) > 0))
	bebealn = str_match(elementExtract(info(vcf[begr$sourceId])$BEALN, 1), "^([^:]+):([0-9]+)[|]([-+])[|]([^|]+)[|]([0-9]+)")
	bebealngr = GRanges(seqnames=bebealn[,2], ranges=IRanges(start=as.integer(bebealn[,3]), width=GenomicAlignments::cigarWidthAlongReferenceSpace(bebealn[,5])))
	callerdf = bind_rows(
		data.frame(
			grid=names(bpanchorgr),
			caller_centro_chr=as.character(seqnames(bpcentrogr)),
			caller_centro_pos=start(bpcentrogr),
			caller_centro_strand=as.character(strand(bpcentrogr)),
			caller_centro_repeat_type =hg19_repeats_gr$repeatType [findOverlaps(bpcentrogr, hg19_repeats_gr, select="first")],
			caller_centro_repeat_class=hg19_repeats_gr$repeatClass[findOverlaps(bpcentrogr, hg19_repeats_gr, select="first")],
			is_centro_hit=bpcentrogr$is_centro_hit),
		data.frame(
			grid=names(begr),
			caller_centro_repeat_type =hg19_repeats_gr$repeatType [findOverlaps(bebealngr, hg19_repeats_gr, select="first")],
			caller_centro_repeat_class=hg19_repeats_gr$repeatClass[findOverlaps(bebealngr, hg19_repeats_gr, select="first")]))
	fgr = gr[callerdf$grid]
	callerdf = callerdf	%>%
		mutate(
			caller=caller,
			sourceid=fgr$sourceId,
			runid=info(vcf[fgr$sourceId])$RUNID,
			filters=VariantAnnotation::fixed(vcf[fgr$sourceId])$FILTER,
			caller_anchor_chr=as.character(seqnames(fgr)),
			caller_anchor_pos=start(fgr),
			caller_anchor_strand=as.character(strand(fgr)),
			caller_anchor_repeat_type =hg19_repeats_gr$repeatType [findOverlaps(fgr, hg19_repeats_gr, select="first")],
			caller_anchor_repeat_class=hg19_repeats_gr$repeatClass[findOverlaps(fgr, hg19_repeats_gr, select="first")],
			is_breakend_call=is.na(fgr$partner)) %>%
		inner_join(truthdf, by="runid") %>%
		replace_na(list(
			is_breakpoint_match=FALSE,
			is_breakend_match=FALSE,
			is_centro_hit=FALSE,
			centro_repeat_type="",
			centro_repeat_class="",
			nchor_repeat_type="",
			anchor_repeat_class="",
			caller_centro_repeat_type="",
			caller_centro_repeat_class="",
			caller_anchor_repeat_type="",
			caller_anchor_repeat_class="")) %>%
		mutate(
			centro_repeat_type_match=caller_centro_repeat_type == centro_repeat_type,
			centro_repeat_class_match=caller_centro_repeat_class == centro_repeat_class,
			is_anchor_hit=fgr$is_anchor_hit,
			is_breakpoint_match=is_anchor_hit & is_centro_hit,
			is_breakend_match=is_anchor_hit) %>%
		group_by(runid) %>%
		arrange(is_breakpoint_match, is_breakend_match, is_breakend_call, centro_repeat_type_match, centro_repeat_class_match) %>%
		mutate(tp=(is_breakpoint_match | is_breakend_match) & cumsum(is_breakpoint_match | is_breakend_match) == 1) %>%
		ungroup() %>%
		mutate(
			status=
				ifelse(!tp, paste0("False Positive"),#, ifelse(is_breakpoint_match | is_breakend_match, " (duplicate call)", "")),
					ifelse(is_breakpoint_match, "Matching breakpoint call",
						ifelse(is_breakend_call, "Matching breakend call",
							"Breakpoint call (one side correct)"))))
}
callsdf = bind_rows(lapply(c("gridss", "manta", "svaba", "novobreak"), process_caller))
resultsdf = bind_rows(
	callsdf,
	# fill in false negatives
	data.frame(runid=unique(truthdf$runid), caller=unique(callsdf$caller)) %>%
		expand(runid, caller) %>%
		anti_join(callsdf %>% filter(tp), by=c("runid", "caller")) %>%
		inner_join(truthdf, by="runid") %>%
		mutate(status="False Negative")) %>%
	mutate(status=factor(status, levels=c(
		"False Positive",
		"False Negative",
		"Breakpoint call (one side correct)",
		"Matching breakend call",
		"Matching breakpoint call")))
		
# Plotting time!
result_by_rc = resultsdf %>%
	mutate(print_centro_repeat_class = summaryRepeatClass(centro_repeat_class)) %>%
	group_by(caller, print_centro_repeat_class) %>%
	mutate(
		total_events=sum(status!="False Positive"),
		total_tp=sum(!str_detect(status, "False"))) %>%
	group_by(caller, print_centro_repeat_class, status, total_tp, total_events) %>%
	summarise(events=n(), portion=n()/ifelse(status=="False Positive", total_tp, total_events)) %>%
	distinct()
ggplot(result_by_rc) +
	aes(x=caller, y=portion, fill=status) +
	geom_col() +
	scale_fill_manual(values=c(
		"#99000d",
		"#e41a1c",
		"#a6cee3",
		"#1f78b4",
		"#4daf4a")) +
	facet_grid(status!="False Positive" ~ print_centro_repeat_class) +
	scale_y_continuous(expand=c(0,0,0,0), labels = scales::percent) + 
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
figsave("sim_centromere", width=7, height=4)

# Figures in manuscript main text
result_by_rc %>% filter(print_centro_repeat_class =="Satellite") %>% arrange(status)
save.image("sim.RData")
