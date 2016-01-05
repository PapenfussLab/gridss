##################
# Processing steps (approximate runtime: 186h wall , 243h CPU)
##################
# - Follow na12878.R processing steps
# - Create symbolic links from ~/i/data.na12878/00000000000000000000000000000000.* to ~/i/data.assembly/00000000000000000000000000000000.*
# - ./call_gridss.sh assembly
# - Run this script

source("libna12878.R")

kvcfs <- NULL
pwd <- getwd()
setwd(paste0(rootdir, "i/data.assembly"))
kmetadata <- LoadMetadata()
kvcfs <- LoadVcfs(kmetadata, existingVcfs=kvcfs)
setwd(pwd)

kvcfs <- lapply(kvcfs, function(vcf) {
  vcf <- vcf[isDeletionLike(vcf, minsize), ]
  if (is.null(attr(vcf, "sourceMetadata")) || is.na(attr(vcf, "sourceMetadata")$CX_CALLER) || str_split(attr(vcf, "sourceMetadata")$CX_CALLER, "/")[[1]][1] !="gridss") return(vcf)
  # assembly calls only
  df <- gridss.vcftodf(vcf)
  score <- df$ASQ + df$RASQ
  rowRanges(vcf)$QUAL <- score
  vcf <- vcf[score > 0,]
  return(vcf)
})

###########
# Gridss assembly approaches
#assdelroc <- bedToROC(longReadBed(gridssassemblyvcfs))
#assdelroc$Filter <- str_extract(assdelroc$caller, "[^$]*$")
#assdelroc$caller <- str_extract(assdelroc$caller, "[^$]*")
#ggplot(assdelroc) + aes(x=fp, y=tp, color=caller, linetype=Filter) + geom_line() +
#coord_cartesian(xlim=c(0, 1000), ylim=c(0, 3000)) + 
# scale_color_brewer(palette="Set2") +
#  labs(x="False Positives", y="True Positives", title="gridss assembly comparison")
#saveplot(paste0("na12878_assembly_tp_fp_pacbiomoleculo", "_error_", maxLengthRelativeDifference, "_", maxPositionDifference_LengthOrEndpoint, ""), width=7, height=5)
#ggplot(assdelroc) + aes(x=tp, y=precision, color=caller, linetype=Filter) + geom_line() + 
#  scale_color_brewer(palette="Set2") +
#  labs(x="True Positives", y="Precision", title="gridss assembly comparison")
#saveplot(paste0("na12878_assembly_prec_pacbiomoleculo", "_error_", maxLengthRelativeDifference, "_", maxPositionDifference_LengthOrEndpoint, ""), width=7, height=5)


###########
# Gridss assembly and kmer length comparison
kdelroc <- bedToROC(longReadBed(kvcfs))
kdelroc$assembly <- as.character(kdelroc$assembly)
kdelroc$assembly[kdelroc$assembly=="Subgraph"] <- "Windowed"
kdelroc$kmer <- as.character(kdelroc$kmer)
ggplot(kdelroc[kdelroc$Filter=="Including Filtered"]) + aes(x=fp, y=tp, color=kmer, linetype=assembly) + geom_line() +
  coord_cartesian(xlim=c(0, 1000), ylim=c(0, 3000)) + 
  scale_color_brewer(palette="Set2") +
  labs(x="False Positives", y="True Positives", title="gridss assembly kmer size comparison")
saveplot(paste0("na12878_gridss_kmer_tp_fp_pacbiomoleculo", "_error_", maxLengthRelativeDifference, "_", maxPositionDifference_LengthOrEndpoint, ""), width=7, height=5)

# only one data point per tpotherwise illustrator crashs
ggplot(kdelroc[kdelroc$Filter=="Including Filtered"][, .SD[!duplicated(.SD$tp, fromLast=TRUE),], by=c("kmer", "assembly")]) +
  aes(x=tp, y=precision, color=kmer, linetype=assembly) + 
  geom_line() + 
  scale_color_brewer(palette="Set2") +
  labs(x="True Positives", y="Precision", title="gridss assembly kmer size comparison")
saveplot(paste0("na12878_gridss_kmer_prec_pacbiomoleculo", "_error_", maxLengthRelativeDifference, "_", maxPositionDifference_LengthOrEndpoint, ""), width=7, height=5)
