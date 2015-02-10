source("libneochromosome.R")

sample <- "778"
rp <- getrpcalls("C:/dev/neochromosome/mmc3.xlsx", "778 (DR)")
cgr <- getcgr("C:/dev/neochromosome/mmc4.xlsx", "778_CGRs")
cn <- getcn("C:/dev/neochromosome/mmc4.xlsx", "778_CN")
vcf <- readVcf("W:/778/778.vcf", "hg19")

sample <- "T1000"
rp <- getrpcalls("C:/dev/neochromosome/mmc3.xlsx", "T1000 (DR)")
cgr <- getcgr("C:/dev/neochromosome/mmc4.xlsx", "T1000_CGRs")
cn <- getcn("C:/dev/neochromosome/mmc4.xlsx", "T1000_CN")
vcf <- readVcf("W:/T1000/T1000.vcf", "hg19")

sample <- "GOT3"
rp <- getrpcalls("C:/dev/neochromosome/mmc3.xlsx", "GOT3 (DR)")
cgr <- getcgr("C:/dev/neochromosome/mmc4.xlsx", "GOT3_CGRs")
cn <- getcn("C:/dev/neochromosome/mmc4.xlsx", "GOT3_CN")
vcf <- readVcf("W:/Papenfuss_lab/projects/liposarcoma/data/gridss/GOT3/GOT3.vcf", "hg19")



vcf <- gridss.removeUnpartnerededBreakend(vcf)
vcfdf <- gridss.truthdetails.processvcf.vcftodf(vcf)
rpmaxgap=120
hits <- countVcfGrBreakpointHits(vcf, rp, maxgap=rpmaxgap)
vcfdf$rpHits <- hits$queryHits
vcfdf$hits_all_gap120 <- hits$subjectHits
rp$hits_hc_gap120 <- countVcfGrBreakpointHits(vcf[vcfdf$confidence=="High",], rp, maxgap=rpmaxgap)$subjectHits
rp$hits_mc_gap120 <- countVcfGrBreakpointHits(vcf[vcfdf$confidence=="High" | vcfdf$confidence=="Medium",], rp, maxgap=rpmaxgap)$subjectHits
rp$hits_gap120 <- as.factor((rp$hits_all_gap120 > 0) + (rp$hits_mc_gap120 > 0) + (rp$hits_hc_gap120 > 0))
levels(rp$hits_gap120) <- c("None", "Low", "Medium", "High")

###############
# Read Pair variant calling concordance
###############
#library(pROC) #install.packages("pROC")
dfroc <- NULL
for (confidence in as.numeric(unique(vcfdf$confidence))) {
  subset <- as.numeric(vcfdf$confidence) >= confidence
  dftmp <- toroc(vcfdf$QUAL[subset], vcfdf$rpHits[subset] > 0, totaltp=length(rp))
  dftmp$confidence <- levels(vcfdf$confidence)[confidence]
  dfroc <- rbind(dfroc, dftmp)  
}
ggplot(dfroc, aes(x=log(qual), y=sens, color=confidence)) + geom_line() + scale_x_reverse() + labs(title=paste0("Sensitivity of curated RP call detection - ", sample))
ggsave(paste0(sample, "_rp_roc.png"), width=10, height=7.5)

ggplot(as.data.frame(mcols(rp)), aes(x=nreads, fill=hits_gap120)) +
  geom_histogram() +
  scale_x_log10() + 
  labs(title=paste0("Sensitivity of curated RP call detection - ", sample))
ggsave(paste0(sample, "_rp_hist.png"), width=10, height=7.5)
