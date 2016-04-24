##################
# Processing steps
##################
# - Follow na12878.R processing steps
# - Create symbolic links from ~/i/data.na12878/00000000000000000000000000000000.* to ~/i/data.model/00000000000000000000000000000000.*
# - ./call_gridss.sh model
# - Run this script
source("common.R")
source("libna12878.R")

mvcfs <- NULL
pwd <- getwd()
setwd(paste0(rootdir, "i/data.model"))
mmetadata <- LoadMetadata()
mvcfs <- LoadVcfs(mmetadata, existingVcfs=mvcfs)
setwd(pwd)

mvcfs <- lapply(mvcfs, function(vcf) {
  vcf <- vcf[isDeletionLike(vcf, minsize), ]
  return(vcf)
})

mlrbed <- longReadBed(mvcfs)
mdelroc <- bedToROC(mlrbed)
mdelroc$assembly <- as.character(mdelroc$assembly)
mdelroc$kmer <- as.character(mdelroc$kmer)
mdelroc$model <- as.character(mdelroc$model)
mdelroc$exclusion <- as.character(mdelroc$exclusion)
mdelroc$subset <- ifelse(mdelroc$exclusion == "None", "SC RP", ifelse(mdelroc$exclusion == "SC", "RP", "SC"))
mdelroc$subset <- paste0(mdelroc$subset, ifelse(mdelroc$assembly=="Positional", " AS", ""))


ggplot(mdelroc) +
  aes(x=fp, y=tp, color=model, linetype=Filter) +
  facet_wrap(exclusion ~ assembly) +
  geom_line() +
  coord_cartesian(xlim=c(0, 1000), ylim=c(0, 3000)) + 
  scale_color_brewer(palette="Set2") +
  labs(x="False Positives", y="True Positives", title="GRIDSS model comparison")
saveplot(paste0("model_breakdown_tp_fp_pacbiomoleculo", "_error_", maxLengthRelativeDifference, "_", maxPositionDifference_LengthOrEndpoint, ""), width=7, height=5)

# only one data point per tp otherwise illustrator crashs
ggplot(mdelroc) + #[, .SD[!duplicated(.SD$tp, fromLast=TRUE),], by=c("kmer", "assembly", "model", "exclusion", "Filter")]) +
  aes(x=tp, y=precision, color=model, linetype=Filter) + 
  facet_grid(exclusion ~ assembly) +
  geom_line() + 
  scale_color_brewer(palette="Set2") +
  labs(x="True Positives", y="Precision", title="GRIDSS model comparison")
saveplot(paste0("model_breakdown_prec_pacbiomoleculo", "_error_", maxLengthRelativeDifference, "_", maxPositionDifference_LengthOrEndpoint, ""), width=7, height=5)

mcmp <- mdelroc[mdelroc$Filter==text_all_calls & mdelroc$assembly=="Positional" & mdelroc$exclusion=="None",]

ggplot(mcmp) +
  aes(x=fp, y=tp, color=model) +
  geom_line() +
  coord_cartesian(xlim=c(0, 1000), ylim=c(0, 3000)) + 
  scale_color_brewer(palette="Set2") +
  labs(x="False Positives", y="True Positives", title="GRIDSS model comparison")
saveplot(paste0("model_tp_fp_pacbiomoleculo", "_error_", maxLengthRelativeDifference, "_", maxPositionDifference_LengthOrEndpoint, ""), width=7, height=5)

# only one data point per tp otherwise illustrator crashs
ggplot(mcmp[, .SD[!duplicated(.SD$tp, fromLast=TRUE),], by=c("kmer", "assembly", "model", "exclusion")]) +
  aes(x=tp, y=precision, color=model) + 
  geom_line() + 
  scale_color_brewer(palette="Set2") +
  labs(x="True Positives", y="Precision", title="GRIDSS model comparison")
saveplot(paste0("model_prec_pacbiomoleculo", "_error_", maxLengthRelativeDifference, "_", maxPositionDifference_LengthOrEndpoint, ""), width=7, height=5)

maxfp <- 1000
auc <- mcmp[, list(auc=rocauc(.SD, maxfp=maxfp)), by=c("kmer", "assembly", "model", "exclusion")]
lowertriangle <- max(mcmp[mcmp$fp <= maxfp]$tp) * maxfp / 2
auc$relativePerformance <- (auc$auc - lowertriangle) / max(auc$auc - lowertriangle)

# breakdown
breakdown <- mdelroc[mdelroc$Filter==text_all_calls & mdelroc$model=="FastEmpiricalReferenceLikelihood",]
breakdown$Evidence <- factor(ifelse(breakdown$exclusion=="None", "Both", ifelse(breakdown$exclusion=="RP", "Split Read", "Discordant Pairs")))
breakdown$Assembly <- relevel(factor(breakdown$assembly), "Positional")

ggplot(breakdown) +
  aes(x=fp, y=tp, color=Evidence, linetype=Assembly) +
  geom_line() +
  coord_cartesian(xlim=c(0, 1000), ylim=c(0, 3000)) + 
  scale_color_brewer(palette="Dark2") +
  labs(x="False Positives", y="True Positives", title="GRIDSS support breakdown")
saveplot(paste0("breakdown_tp_fp_pacbiomoleculo", "_error_", maxLengthRelativeDifference, "_", maxPositionDifference_LengthOrEndpoint, ""), width=7, height=5)

# only one data point per tp otherwise illustrator crashs
ggplot(breakdown[, .SD[!duplicated(.SD$tp, fromLast=TRUE),], by=c("kmer", "assembly", "model", "exclusion")]) +
  aes(x=tp, y=precision, color=Evidence, linetype=Assembly) + 
  geom_line() + 
  scale_color_brewer(palette="Dark2") +
  labs(x="True Positives", y="Precision", title="GRIDSS support breakdown")
saveplot(paste0("breakdown_prec_pacbiomoleculo", "_error_", maxLengthRelativeDifference, "_", maxPositionDifference_LengthOrEndpoint, ""), width=7, height=5)



#gridssbreakdownvcfs <- unlist(recursive=FALSE, lapply(vcfs, function(vcf) {
#   if (is.null(attr(vcf, "sourceMetadata")) || is.na(attr(vcf, "sourceMetadata")$CX_CALLER) || str_split(attr(vcf, "sourceMetadata")$CX_CALLER, "/")[[1]][1] !="gridss") return(NULL)
#   if (!is.na(attr(vcf, "sourceMetadata")$CX_CALLER_ARGS) && !is.null(attr(vcf, "sourceMetadata")$CX_CALLER_ARGS) && attr(vcf, "sourceMetadata")$CX_CALLER_ARGS == "assembly.method") return (NULL) # exclude subgraph assembly
#   df <- gridss.vcftodf(vcf)
#   return (lapply(list(
#     list(f=function(df) df$QUAL, label="DP SR Assembly (Default)"),
#     list(f=function(df) df$ASQ + df$RASQ, label="Assembly only"),
#     list(f=function(df) df$RPQ + df$SRQ + df$RSRQ, label="DP SR"),
#     list(f=function(df) df$RPQ, label="DP only"),
#     list(f=function(df) df$SRQ + df$RSRQ, label="SR only"),
#     list(f=function(df) df$ASRP + df$ASSR, label="Assembly only$count"),
#     list(f=function(df) df$RP + df$SR + df$RSR, label="DP SR$count"),
#     list(f=function(df) df$RP, label="DP only$count"),
#     list(f=function(df) df$SR + df$RSR, label="SR only$count")),
#     function(scoring) {
#       ovcf <- vcf
#       rowRanges(ovcf)$QUAL <- scoring$f(df)
#       rowRanges(ovcf)$FILTER <- "."
#       attr(ovcf, "sourceMetadata")$CX_CALLER <- scoring$label
#       return(ovcf)
#     }))
# }))

