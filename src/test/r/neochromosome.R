##################
# Processing steps
##################
# 1) download Socrates supplementary materials spreadsheets # PMID: 25517748
# 2) download fastq http://www.ebi.ac.uk/ena/data/view/ERP004006 # PMID: 25517748
# 3) Align against hg19. 778, T1000 aligned using bowtie2, GOT3 aligned with bwa. One BAM per library
# 4) Remove duplicates using picard tools.
# 5) run gridss with one input file per library. Explicit fragment size threshold for 778 were {500,350,350,750,700,700}
# 6) run this script
# 7) generate excel pivot table from dtsummarised (neo_summarised_200.csv)

#source("http://bioconductor.org/biocLite.R")
#biocLite("rtracklayer")
library(rtracklayer)
#setwd("W:/dev/gridss/src/test/R/") # location of this script
source("libgridss.R")
source("libneochromosome.R")
source("libvcf.R")
source("common.R")

socratesSupplementaryMaterialsLocation <- "C:/dev/neochromosome/"
location <- list(
  sample778="W:/Papenfuss_lab/projects/liposarcoma/data/gridss/778/778.vcf",
  sampleT1000="W:/Papenfuss_lab/projects/liposarcoma/data/gridss/T1000/T1000.vcf",
  sampleGOT3="W:/Papenfuss_lab/projects/liposarcoma/data/gridss/GOT3/bt2/GOT3-bt2.vcf")

minqual <- 250
matchmaxgaps <- c(200) #c(1, 10, 50, 100, 200))
bed <- NULL
gridss <- NULL
graphs <- TRUE

for (sample in c("778", "T1000", "GOT3")) {
  rp <- getrpcalls(paste0(socratesSupplementaryMaterialsLocation, "mmc3.xlsx"), paste0(sample," (DR)")) 
  cgr <- getcgr(paste0(socratesSupplementaryMaterialsLocation, "mmc4.xlsx"), paste0(sample,"_CGRs"))
  cn <- getcn(paste0(socratesSupplementaryMaterialsLocation, "mmc4.xlsx"), paste0(sample,"_CN"))
  vcf <- readVcf(location[[paste0("sample", sample)]], "hg19")
  for (matchmaxgap in matchmaxgaps) {
    out <- go(sample, vcf, rp, cgr, minimumEventSize=500, matchmaxgap=matchmaxgap, graphs=graphs, minqual=minqual)
    bed <- rbind(bed, as.data.frame(out$bed))
    gridss <- rbind(gridss, out$gridss)
  }
}
gridss$lessthan7rp <- gridss$RP < 7
dtcalls <- rbind(data.frame(sample=bed$sample, assembly=ifelse(bed$spanning, "Spanning assembly", as.character(bed$assembly)), hit="Published"),
                 data.frame(sample=gridss[is.na(gridss$bedid),]$sample, assembly=gridss[is.na(gridss$bedid),]$assembly, hit="gridss additional"))
# convert from breakend to breakpoint counts by trimming every second row
dtcalls <- dtcalls[order(dtcalls$sample, dtcalls$assembly, dtcalls$hit),][seq(from=1, to=nrow(dtcalls), by=2),]
ggplot(dtcalls) +
  aes(x=hit, fill=assembly) +
  facet_wrap(~ sample) +
  geom_bar() +
  labs("Concordance with neochromosome", x="", y="Count", fill="Support")
saveplot(paste0("neo_summarised_", matchmaxgap, ""))

dtsummarised <- rbind(
  data.table(bed)[, list(count=.N, hit=TRUE), by=list(sample, matchmaxgap, assembly, spanning)],
  data.table(gridss[is.na(gridss$bedid),])[, list(count=.N, hit=FALSE, spanning=NA), by=list(sample, matchmaxgap, assembly)])
dtsummarised$count <- dtsummarised$count / 2 # breakend -> breakpoint conversion
write.csv(dtsummarised, file=paste0("neo_summarised_", matchmaxgap, ".csv"))
ggplot(dtsummarised) + aes(x=hit, y=count, fill=assembly, color=spanning) + geom_bar()

# portion less than RP threshold = most
table(gridss[is.na(gridss$bedid),]$RP >= 7, gridss[is.na(gridss$bedid),]$assembly)
table(gridss[is.na(gridss$bedid),]$RP >= 7, gridss[is.na(gridss$bedid),]$sample) # breakdown by sample -> most less than threshold
table(gridss[is.na(gridss$bedid),]$assembly, gridss[is.na(gridss$bedid),]$lessthan7rp) / 2
table(gridss[is.na(gridss$bedid) & gridss$RP >= 7,]$assembly, gridss[is.na(gridss$bedid) & gridss$RP >= 7,]$sample)
data.table(gridss[gridss$assembly=="Both" & !is.na(gridss$bedid),])[, list(count=.N, homology=sum(HOMLEN>0), untemplated=sum(nchar(INSSEQ)>0)), by=list(sample, matchmaxgap, assembly)]
length(gridss[nchar(gridss$INSSEQ)>50,]$INSSEQ)

table(gridss[is.na(gridss$bedid),]$RP >= 7 & gridss[is.na(gridss$bedid),]$asWithin1kbp == 1 & gridss[gridss[is.na(gridss$bedid),]$mate,]$asWithin1kbp == 1) / 2
table(gridss[is.na(gridss$bedid),]$RP >= 7 & gridss[is.na(gridss$bedid),]$asWithin1kbp == 1 & gridss[gridss[is.na(gridss$bedid),]$mate,]$asWithin1kbp == 1, gridss[is.na(gridss$bedid),]$assembly) / 2


#####################
# assembly-only identification
gr <- gridss[!gridss$found & gridss$RP<1 & gridss$SR < 1 & gridss$RSR < 1,]$POS
gr <- GRanges(unlist(lapply(str_split(gr, ":"), function (x) x[1])), IRanges(width=1, start=unlist(lapply(str_split(gr, ":"), function (x) as.integer(x[2])))))

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
saveplot("neo_rpsupport")
  


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



# edge Analysis
#dtg <- data.table(gridss[gridss$assembly == "Both",])
#dtbed <- data.table(as.data.frame(bed[bed$assembly=="Both",]))
#dte <- merge(dtg, dtbed, by.x=c("variantid", "sample"), by.y=c("gridssid", "sample"))
# TODO: why are both the gridssid for the same size?
vcf <- gridss.removeUnpartnerededBreakend(vcf)
sdf <- gridss.vcftodf(vcf)
sdf$chr <- seqnames(rowRanges(vcf))
sdf$strand <- ifelse(str_detect(as.character(rowRanges(vcf)$ALT), "[[:alpha:]]+(\\[|]).*(\\[|])"), "+", "-")

sdf[sdf$AS>0 & sdf$RAS>0,]

ggplot(sdf[sdf$AS>0 & sdf$RAS>0,]) + aes(x=size, y=QUAL, color=strand) + geom_point() + scale_x_log10() + scale_y_log10()


