library(ggplot2)
library(RColorBrewer)
library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(stringr)
setwd("C:/dev/idsv/src/main/R/")
source("libgridss.R")

# gridss TODO: split out *_RM into madfing to us, and madfing elsewhere
#vcf <- readVcf("C:/dev/778.vcf", "hg19_random")
bed <- import.bed(con="W:/778/idsv/cgrs-from-table3.bed")
socBed <- import.bed(con="C:/dev/tracks/778_breakpoints_SOCRATES.bed")
rpBed <- import.bed(con="C:/dev/tracks/778_breakpoints_v15.bed")
vcf <- readVcf("W:/778/idsv/778.vcf", "hg19_random")
df <- gridss.truthdetails.processvcf.vcftodf(vcf)
df$cgr <- gridss.overlaps(vcf, bed)
#df$socrates <- gridss.overlaps(vcf, socBed, maxgap=4)
#df$rpBed <- gridss.overlaps(vcf, socBed, maxgap=16)
df$hasSC <- paste("SC", ifelse(df$SC > 0 & df$RSC > 0, "Both", ifelse(df$SC > 0, "local", ifelse(df$RSC > 0, "remote", "znone"))))
df$hasAS <- paste("AS", ifelse(df$AS > 0 & df$RAS > 0, "Both", ifelse(df$AS > 0, "local", ifelse(df$RAS > 0, "remote", "znone"))))
df$hasRP <- ifelse(df$RP > 0, "RP", "RP None")
# set remote breakend columns
df$cgrmate <- df[df$mateid,]$cgr

##################
# Sanity checks
##################
# should all have mates
df[!(df$mateid %in% row.names(df)),]
mdf <- df[df$mateid,]
# breakpoint fields should match on both sides of the breakend
df[df$FILTER != mdf$FILTER,]
df[df$SOMATIC != mdf$SOMATIC,]
df[df$IMPRECISE != mdf$IMPRECISE,]
df[df$SVLEN != mdf$SVLEN,]
df[df$SVTYPE != mdf$SVTYPE,]
df[df$SPV != mdf$SPV,]
df[df$QUAL != mdf$QUAL,]
df[df$CQ != mdf$CQ,]
df[df$AS != mdf$RAS,]
df[df$RP != mdf$RP,]
df[df$SC != mdf$RSC,]
df[abs(df$ASQ - mdf$RASQ) > 0.1,]
df[abs(df$RPQ - mdf$RPQ) > 0.1,]
df[abs(df$SCQ - mdf$RSCQ) > 0.1,]

# remove those that fail mate check so we can continue
vcf <- vcf[as.character(info(vcf)$MATEID) %in% row.names(vcf),]
callPos <- rowData(vcf)
callPos$mate <- as.character(info(vcf)$MATEID)
strand(callPos) <- ifelse(str_detect(as.character(callPos$ALT), "[[:alpha:]]+(\\[|]).*(\\[|])"), "+", "-")
callPosMate <- callPos[callPos$mate,]

###############
# Read Pair variant calling concordance
###############
# chr12:84262289 call position for RP calls incorrect by 111 bases
rpBedMate <- str_match(rpBed$name, "[[:alnum:]]+:[[:digit:]]+\\([+-]\\)-([[:alnum:]]+):([[:digit:]]+)\\(([+-])\\)")
rpBedMate <- GRanges(seqnames=rpBedMate[,2], ranges=IRanges(start=as.numeric(rpBedMate[,3]), width=1, names=rpBedMate[,1]), strand=rpBedMate[,4])
rpBed <- gridss.annotateBreakpointHits(rpBed, rpBedMate, vcf, maxgap=120)

ggplot(as.data.frame(mcols(rpBed)), aes(x=score, fill=called)) +
  geom_histogram() +
  scale_x_log10() + 
  labs(title="Concordance with existing read pair based calls")
ggsave("778_rp_histogram.png", width=10, height=7.5)

ggplot(as.data.frame(mcols(rpBed)), aes(x=score, y=qual, color=called)) +
  geom_point() +
  scale_x_log10() + scale_y_log10() + 
  labs(title="Concordance with existing read pair based calls")
ggsave("778_rp_concordance.png", width=10, height=7.5)

###############
# Socrates variant calling concordance
###############
socBedMate <- str_match(socBed$name, "bp[[:digit:]]+_to_bp[[:digit:]]+_([[:alnum:]]+):([[:digit:]]+)-([[:digit:]]+)")
socBedMate <- GRanges(seqnames=socBedMate[,2], ranges=IRanges(start=as.numeric(socBedMate[,3]), end=as.numeric(socBedMate[,4]), names=socBedMate[,1]), strand="*")
socBed <- gridss.annotateBreakpointHits(socBed, socBedMate, vcf, maxgap=16, ignore.strand=TRUE) # TODO: consistent strand meaning

ggplot(as.data.frame(mcols(socBed)), aes(x=score, fill=called)) +
  geom_histogram() +
  scale_x_log10() + 
  labs(title="Concordance with existing Socrates calls")
ggsave("778_soc_histogram.png", width=10, height=7.5)

ggplot(as.data.frame(mcols(socBed)), aes(x=score, y=qual, color=called)) +
  geom_point() +
  scale_x_log10() + scale_y_log10() + 
  labs(title="Concordance with existing Socrates calls")
ggsave("778_soc_concordance.png", width=10, height=7.5)






###############
# Exploratory plots
###############

# Contribution of local breakend evidence
ggplot(df, aes(x=QUAL, y=BQ, color=factor(pmin(AS, 1)+pmin(RAS, 1)), size=BAS+1)) + geom_point() + scale_x_log10() + scale_y_log10() + facet_grid(cgrmate ~ cgr) +
  labs(title="Breakpoint and breakend quality score distributions according to breakend CGR location")

# Effect of unique read madfing
# TODO: why to points exist where QUAL > CQ ? this should not be possible
ggplot(df, aes(x=QUAL, y=CQ, color=factor(pmin(AS, 1)+pmin(RAS, 1)))) + geom_point() + scale_x_log10() + scale_y_log10()
# histogram of 
ggplot(df, aes(x=QUAL/CQ)) + geom_histogram()# + scale_y_log10()
head(df[df$QUAL > df$CQ & df$CQ > 100,])

# evidence counts
ggplot(df, aes(x=RP, y=SC+RSC, color=factor(pmin(AS, 1)+pmin(RAS, 1)))) + facet_grid(cgrmate ~ cgr) + geom_point() + scale_x_log10() + scale_y_log10() + geom_jitter(position = position_jitter(width = 1, height=1))

table(data.frame(a=df$cgr, b=df$cgrmate)) /2 

# qual dist
# AS remote can be removed due to symmetry
ggplot(df[df$hasAS != "AS remote" & !(df$hasSC == "SC remote" & df$hasAS != "AS local"), ],
    aes(x=QUAL)) +
  facet_grid(hasRP + hasAS ~ hasSC) + 
  geom_histogram() + scale_x_log10() + scale_y_log10() +
  theme_bw() + labs(title="Breakpoint quality distribution by evidence type")
ggplot(df, aes(x=QUAL, y=BQ)) + facet_grid(hasRP + hasSC ~ hasAS) + geom_point() + scale_x_log10() + scale_y_log10()

# local anchor length shouldn't mean much
ggplot(df, aes(x=A_BLRM, y=A_BLLM)) + geom_point()

# should have a linear relationship
ggplot(df, aes(x=RC+1, y=PC+1, color=log(QUAL))) + geom_smooth() + geom_jitter() + geom_point() + scale_y_log10() + scale_x_log10() + scale_colour_gradientn(colours=c("red", "blue"))

ggplot(df, aes(x=SC + RSC, y=REF, color=REFPAIR)) + geom_point()

# assembly length vs breakpoint quality
ggplot(df, aes(x=A_BLRM, y=QUAL)) + geom_point()

# coverage
ggplot(df, aes(x=REF)) + geom_histogram() + scale_x_log10()
ggplot(df, aes(x=REFPAIR)) + geom_histogram() + scale_x_log10()
ggplot(df, aes(x=SC)) + geom_histogram() + scale_x_log10()
ggplot(df, aes(x=RSC)) + geom_histogram() + scale_x_log10()
ggplot(df, aes(x=AS)) + geom_histogram() + scale_x_log10()
ggplot(df, aes(x=RAS)) + geom_histogram() + scale_x_log10()
ggplot(df, aes(x=RP)) + geom_histogram() + scale_x_log10()


# distribution of variant evidence contribution
ggplot(df, aes(x=SC, y=REF, color=1)) + geom_point() + scale_y_log10() + scale_x_log10() + scale_colour_gradientn(colours=rainbow(4)) 
# normal coverage
ggplot(df, aes(x=RCNormal + SCECNormal)) + geom_histogram()

# Effect of local breakend evidence & unique evidence assignment
ggplot(df, aes(x=QUAL, y=CQUAL)) + geom_point() + scale_x_log10() + scale_y_log10()




# QUAL distribution by breakpoint evidence
ggplot(df, aes(x=RP + SC + RSC)) + geom_histogram()

# total evidence distribution
ggplot(df, aes(x=log10(RPEC+A_RP+1), y=log10(SCEC+A_SC+1), color=factor(A_EC))) +  geom_point() + ggtitle("All evidence") + scale_x_continuous(limits=c(0, 3.5)) + scale_y_continuous(limits=c(0, 3.5))
ggsave("evidence_distribution_all.png")
# non-assembly evidence distribution
ggplot(df, aes(x=log10(RPEC+1), y=log10(SCEC+1), color=factor(A_EC))) + geom_point() + ggtitle("Non-assembly evidence") + scale_x_continuous(limits=c(0, 3.5)) + scale_y_continuous(limits=c(0, 3.5))
ggsave("evidence_distribution_notassembled.png")
# distribution of evidence not in assembly
ggplot(df, aes(x=log10(A_RP+1), y=log10(A_SC+1), color=factor(A_EC))) + geom_point() + ggtitle("Assembly evidence") + scale_x_continuous(limits=c(0, 3.5)) + scale_y_continuous(limits=c(0, 3.5))
ggsave("evidence_distribution_assembled.png")


# how many assemblies where we have no madfed evidence?
nrow(df[df$SCRM==0 & df$RPRM==0 & df$A_RM==0,]) # = number of calls in which *ALL* the breakpoint evidence has been removed -> no sudfort for this variant yet we're calling them
# created breakpoint from assembly
ggplot(df[df$SCRM==0 & df$RPRM==0,], aes(x=A_MQT)) + geom_histogram()
# great assemblies in both directions
head(df[df$SCRM==0 & df$RPRM==0 & df$A_MQT>80,])

