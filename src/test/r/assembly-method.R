##################
# Processing steps (approximate runtime: 186h wall , 243h CPU)
##################
# - Follow na12878.R processing steps
# - Create symbolic links from ~/i/data.na12878/00000000000000000000000000000000.* to ~/i/data.assembly/00000000000000000000000000000000.*
# - ./call_gridss.sh assembly
# - Run this script
source("common.R")
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
  localqual <- info(vcf)$ASQ
  localqual[is.na(localqual)] <- 0
  remotequal <- info(vcf)$RASQ
  remotequal[is.na(remotequal)] <- 0
  rowRanges(vcf)$QUAL <- localqual + remotequal
  vcf <- vcf[rowRanges(vcf)$QUAL > 0,]
  return(vcf)
})

###########
# Gridss assembly and kmer length comparison
klrbed <- longReadBed(kvcfs)
kdelroc <- bedToROC(klrbed)
kdelroc$assembly <- as.character(kdelroc$assembly)
kdelroc$assembly[kdelroc$assembly=="Subgraph"] <- "Windowed"
kdelroc$kmer <- as.character(kdelroc$kmer)
ggplot(kdelroc[kdelroc$Filter==text_all_calls]) + aes(x=fp, y=tp, color=kmer, linetype=assembly) + geom_line() +
  coord_cartesian(xlim=c(0, 1000), ylim=c(0, 3000)) + 
  scale_color_brewer(palette="Set2") +
  labs(x="False Positives", y="True Positives", title="gridss assembly kmer size comparison")
saveplot(paste0("na12878_gridss_kmer_tp_fp_pacbiomoleculo", "_error_", maxLengthRelativeDifference, "_", maxPositionDifference_LengthOrEndpoint, ""), width=7, height=5)

# only one data point per tpotherwise illustrator crashs
ggplot(kdelroc[kdelroc$Filter==text_all_calls][, .SD[!duplicated(.SD$tp, fromLast=TRUE),], by=c("kmer", "assembly")]) +
  aes(x=tp, y=precision, color=kmer, linetype=assembly) + 
  geom_line() + 
  scale_color_brewer(palette="Set2") +
  labs(x="True Positives", y="Precision", title="gridss assembly kmer size comparison")
saveplot(paste0("na12878_gridss_kmer_prec_pacbiomoleculo", "_error_", maxLengthRelativeDifference, "_", maxPositionDifference_LengthOrEndpoint, ""), width=7, height=5)

