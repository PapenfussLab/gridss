library(StructuralVariantAnnotation)
library(rtracklayer)

vcf = readVcf("gridss.vcf")

# Export breakpoints to BEDPE
bpgr = breakpointRanges(vcf)
# TODO: add your event filtering here. The default GRIDSS output is very verbose/sensitive.
write.table(breakpointgr2bedpe(bpgr), file="gridss_breakpoints.bedpe", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
	
# Export single breakends to BED
begr = breakendRanges(vcf)
# TODO: add your event filtering here. The default GRIDSS output is very verbose/sensitive.
begr$score = begr$QUAL
export(begr, con="gridss_single_breakends.bed")
