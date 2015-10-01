#source("http://bioconductor.org/biocLite.R")
#biocLite("rtracklayer")
library(rtracklayer)
source("libgridss.R")
source("libneochromosome.R")
source("libvcf.R")

#setwd("W:/dev/gridss/src/test/R/")

theme_set(theme_bw())

minqual <- 140 # 20 QUAL * 7 reads
matchmaxgaps <- c(200) #c(1, 10, 50, 100, 200))
bed <- NULL
gridss <- NULL
graphs <- TRUE

sample <- "778"
rp <- getrpcalls("C:/dev/neochromosome/mmc3.xlsx", "778 (DR)")
# no need to explicitly remove these calls as they're not in CGRs
#rp <- rp[!(seqnames(rp) == seqnames(rp[rp$mate,]) & seqnames(rp) %in% c("chr7", "chr22")),] # remove 7 and 22 intrachromosomal arrangements that had really low limits
cgr <- getcgr("C:/dev/neochromosome/mmc4.xlsx", "778_CGRs")
cn <- getcn("C:/dev/neochromosome/mmc4.xlsx", "778_CN")
vcf <- readVcf("W:/Papenfuss_lab/projects/liposarcoma/data/gridss/778/778.vcf", "hg19")
for (matchmaxgap in matchmaxgaps) {
  out <- go(sample, vcf, rp, cgr, minimumEventSize=500, matchmaxgap=matchmaxgap, graphs=graphs, minqual=minqual)
  bed <- rbind(bed, as.data.frame(out$bed))
  gridss <- rbind(gridss, out$gridss)
}

sample <- "T1000"
rp <- getrpcalls("C:/dev/neochromosome/mmc3.xlsx", "T1000 (DR)")
cgr <- getcgr("C:/dev/neochromosome/mmc4.xlsx", "T1000_CGRs")
cn <- getcn("C:/dev/neochromosome/mmc4.xlsx", "T1000_CN")
vcf <- readVcf("Z:/projects/liposarcoma/data/gridss/T1000/T1000.vcf", "hg19")
for (matchmaxgap in matchmaxgaps) {
  out <- go(sample, vcf, rp, cgr, minimumEventSize=500, matchmaxgap=matchmaxgap, graphs=graphs, minqual=minqual)
  bed <- rbind(bed, as.data.frame(out$bed))
  gridss <- rbind(gridss, out$gridss)
}

sample <- "GOT3"
rp <- getrpcalls("C:/dev/neochromosome/mmc3.xlsx", "GOT3 (DR)")
cgr <- getcgr("C:/dev/neochromosome/mmc4.xlsx", "GOT3_CGRs")
cn <- getcn("C:/dev/neochromosome/mmc4.xlsx", "GOT3_CN")
vcf <- readVcf("W:/Papenfuss_lab/projects/liposarcoma/data/gridss/GOT3/bt2/GOT3-bt2.vcf", "hg19")
for (matchmaxgap in matchmaxgaps) {
  out <- go(sample, vcf, rp, cgr, minimumEventSize=500, matchmaxgap=matchmaxgap, graphs=graphs, minqual=minqual)
  bed <- rbind(bed, as.data.frame(out$bed))
  gridss <- rbind(gridss, out$gridss)
}

dtsummarised <- rbind(
  data.table(bed)[, list(count=.N, hit=TRUE), by=list(sample, matchmaxgap, assembly, spanning)],
  data.table(gridss[is.na(gridss$bedid),])[, list(count=.N, hit=FALSE, spanning=NA), by=list(sample, matchmaxgap, assembly)])
dtsummarised$count <- dtsummarised$count / 2 # breakend -> breakpoint conversion
write.csv(dtsummarised, file=paste0("neo_summarised_", matchmaxgap, ".csv"))

# portion less than RP threshold = most
table(gridss[is.na(gridss$bedid),]$RP >= 7, gridss[is.na(gridss$bedid),]$assembly)
table(gridss[is.na(gridss$bedid),]$RP >= 7, gridss[is.na(gridss$bedid),]$sample) # breakdown by sample -> most less than threshold
table(gridss[is.na(gridss$bedid) & gridss$RP >= 7,]$assembly, gridss[is.na(gridss$bedid) & gridss$RP >= 7,]$sample)
data.table(gridss[gridss$assembly=="Both" & !is.na(gridss$bedid),])[, list(count=.N, homology=sum(HOMLEN>0), untemplated=sum(nchar(INSSEQ)>0)), by=list(sample, matchmaxgap, assembly)]

#####################
# realignment rates


#####################
# Working set
#
df <- out$gridss
df$size <- vcftobpgr(out$vcf)$size
ggplot(df) +
  aes(x=distanceToCall, y=QUAL, color=confidence) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap( ~ found)
table(out$bed$hits)
ggplot(as.data.frame(out$bed)) +
  aes(x=distanceToMedHigh, color=is.na(gridssid)) +
  geom_histogram() +
  scale_x_log10() +
  facet_wrap( ~ hits)

out$bed[out$bed$hits=="None",]

table(df[is.na(df$bedid), ]$confidence)
table(df[df$confidence=="High" & is.na(df$bedid), ]$size <= 500, useNA="ifany")

distanceToClosest(rp, rowRanges(vcf))
distanceToClosest(rowRanges(vcf), rp)
rp$distanceToCall
rp[distanceHits]

# Breakdown of small fragment spanning calls by call type
table(out$bed$spanning, out$bed$hits)
print("gridss found ", FindFragments(vcftobpgr(vcf)), "fragments")


ggplot(out$gridss[out$gridss$confidence %in% c("High", "Medium") & is.na(out$gridss$bedid),]) +
  aes(x=RP, color=confidence) + 
  geom_histogram() +
  labs(title="Read Pair support for unmatched gridss calls")
ggsave("rpsupport.png")
  


# Spanning
spanning <- FindFragmentSpanningEvents(vcftobpgr(out$vcf), vcftobpgr(out$vcf))
spanning <- spanning[substr(as.character(spanning$query), 1, nchar(as.character(spanning$query)) - 1) != substr(as.character(spanning$subjectAlt1), 1, nchar(as.character(spanning$subjectAlt1)) - 1) &
                       substr(as.character(spanning$query), 1, nchar(as.character(spanning$query)) - 1) != substr(as.character(spanning$subjectAlt2), 1, nchar(as.character(spanning$subjectAlt2)) - 1),]
table(out$gridss[as.character(spanning$query),]$confidence)

# Fragment size estimation
frag <- FindFragments(vcftobpgr(out$vcf), 100000)
# exclude tandem duplications:
frag[frag$start_local != frag$end_remote & frag$start_local != frag$start_remote & frag$end_local != frag$end_remote & frag$end_local != frag$start_remote,]
frag <- frag[order(frag$size),]
frag <- frag[!duplicated(frag$start_local),]
frag <- frag[!duplicated(frag$end_local),]
ggplot(frag) + aes(x=size) + geom_density()

# emperical distribution:
x <- runif(nrow(out$gridss)/2, min=1, max=100000000)
x <- x[order(x)]
deltax <- c(x, 100000000) - c(0, x)
ggplot() +
  geom_density(data=data.frame(x=deltax), color="red", aes(x=x)) +
  geom_density(data=frag, color="blue", aes(x=size)) +
  scale_x_continuous(limits=c(0, 100000)) +
  labs(x="fragment size")

bed$assembly <- factor(bed$assembly, c("Both", "Single", "No assembly", "None"))
ggplot(as.data.frame(bed), aes(x=sample, color=!spanning, fill=assembly)) + geom_bar()

