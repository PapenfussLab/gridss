library(StructuralVariantAnnotation)
library(rtracklayer)

vcf = readVcf("gridss.vcf")

# Export breakpoints to BEDPE
bpgr = breakpointRanges(vcf)
# TODO: add your event filtering here. The default GRIDSS output is very verbose/sensitive.

bpgr_first_breakend = bpgr[substring(names(bpgr), nchar(names(bpgr))) == "o"]
bpgr_first_breakend_partner = bpgr[bpgr_first_breakend$partner]
bedpedf = data.frame(
	chrom1=seqnames(bpgr_first_breakend),
	start1=start(bpgr_first_breakend),
	end1=end(bpgr_first_breakend),
	chrom2=seqnames(bpgr_first_breakend_partner),
	start2=start(bpgr_first_breakend_partner),
	end2=end(bpgr_first_breakend_partner),
	name=names(bpgr_first_breakend),
	score=bpgr_first_breakend$QUAL,
	strand1=strand(bpgr_first_breakend),
	strand2=strand(bpgr_first_breakend_partner))
write.table(bedpedf, file="gridss_breakpoints.bedpe", sep="\t", quote=FALSE, col.names=FALSE)
	
# Export single breakends to BED
begr = breakendRanges(vcf)
# TODO: add your event filtering here. The default GRIDSS output is very verbose/sensitive.
begr$score = begr$QUAL
export(begr, con="gridss_single_breakends.bed")
