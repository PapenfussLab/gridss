library(ggplot2)
library(RColorBrewer)
setwd("C:/dev/idsv/src/main/R/")
source("libgridss.R")

# gridss TODO: split out *_RM into mapping to us, and mapping elsewhere
vcf <- readVcf("C:/dev/778.vcf", "hg19_random")
#vcf <- readVcf("W:/778/idsv/778.vcf", "hg19_random")
df <- gridss.truthdetails.processvcf.vcftodf(vcf)


# Contribution of local breakend evidence
ggplot(df, aes(x=QUAL, y=BQ, color=factor(pmin(AS, 1)+pmin(RAS, 1)))) + geom_point() + scale_x_log10() + scale_y_log10()

# Effect of unique read mapping
# TODO: why to points exist where QUAL > CQ ? this should not be possible
ggplot(df, aes(x=QUAL, y=CQ, color=factor(pmin(AS, 1)+pmin(RAS, 1)))) + geom_point() + scale_x_log10() + scale_y_log10()
head(df[df$QUAL > 100 * df$CQ,])

# local anchor length shouldn't mean much
ggplot(pp, aes(x=A_BLRM, y=A_BLLM)) + geom_point()

# should have a linear relationship
ggplot(pp, aes(x=RC+1, y=PC+1, color=log(QUAL))) + geom_smooth() + geom_jitter() + geom_point() + scale_y_log10() + scale_x_log10() + scale_colour_gradientn(colours=c("red", "blue"))

ggplot(pp, aes(x=A_SC + A_RP, y=RC, color=PC)) + geom_point()

# assembly length vs breakpoint quality
ggplot(pp, aes(x=A_BLRM, y=QUAL)) + geom_point()

# normal coverage
ggplot(pp, aes(x=RCNormal + SCECNormal)) + geom_histogram()
#total support
ggplot(pp, aes(x=A_SC)) + geom_histogram()
ggplot(pp, aes(x=A_RP)) + geom_histogram()
ggplot(pp, aes(x=PC)) + geom_histogram()
ggplot(pp, aes(x=RC)) + geom_histogram()
ggplot(pp, aes(x=A_SCCLM)) + geom_histogram()
ggplot(pp, aes(x=A_MQRM)) + geom_histogram()
ggplot(pp, aes(x=A_MQRT)) + geom_histogram()
ggplot(pp, aes(x=A_EC)) + geom_histogram()

# distribution of variant evidence contribution
ggplot(df, aes(x=SCEC+RPEC, y=RC+PC, color=log(QUAL))) + geom_point() + scale_y_log10() + scale_x_log10() + scale_colour_gradientn(colours=rainbow(4)) 
# normal coverage
ggplot(df, aes(x=RCNormal + SCECNormal)) + geom_histogram()

# Effect of local breakend evidence & unique evidence assignment
ggplot(df, aes(x=QUAL, y=CQUAL)) + geom_point() + scale_x_log10() + scale_y_log10()




# QUAL distribution by breakpoint evidence
# A_EC, A_RP, A_SC
# A_RM # mapped assembly
# RPEC # inc OEA
# RPRM # DP
# SCEC # inc short SC
# SCRM # mapped SC
# mapped vs unmapped
ggplot(df, aes(x=log10(RPRM+1), y=log10(RPEC-RPRM+1), color=log10(QUAL+1))) + geom_point()  + geom_jitter(position = position_jitter(width=.1, height=.1)) + scale_colour_gradientn(colours=rainbow(4)) 
ggplot(df, aes(x=log10(SCRM+1), y=log10(SCEC-SCRM+1), color=log10(QUAL+1))) + geom_point()  + geom_jitter(position = position_jitter(width=.1, height=.1)) + scale_colour_gradientn(colours=rainbow(4)) 
# mapped evidence level by type
ggplot(df[df$QUAL>=0,], aes(x=log10(RPRM+1), y=log10(SCRM+1), color=factor(A_RM), size=log10(QUAL+1))) + geom_point() + geom_jitter(position = position_jitter(width=.1, height=.1))

# Distribution of BP and BE evidence
ggplot(df, aes(x=log10(BPQUAL+1), y=log10(BEQUAL+1), color=factor(ASSCNT))) + geom_point() + facet_wrap( ~ ASSCNT) + geom_jitter(position = position_jitter(width=.1, height=.1))

# total evidence distribution
ggplot(df, aes(x=log10(RPEC+A_RP+1), y=log10(SCEC+A_SC+1), color=factor(A_EC))) +  geom_point() + ggtitle("All evidence") + scale_x_continuous(limits=c(0, 3.5)) + scale_y_continuous(limits=c(0, 3.5))
ggsave("evidence_distribution_all.png")
# non-assembly evidence distribution
ggplot(df, aes(x=log10(RPEC+1), y=log10(SCEC+1), color=factor(A_EC))) + geom_point() + ggtitle("Non-assembly evidence") + scale_x_continuous(limits=c(0, 3.5)) + scale_y_continuous(limits=c(0, 3.5))
ggsave("evidence_distribution_notassembled.png")
# distribution of evidence not in assembly
ggplot(df, aes(x=log10(A_RP+1), y=log10(A_SC+1), color=factor(A_EC))) + geom_point() + ggtitle("Assembly evidence") + scale_x_continuous(limits=c(0, 3.5)) + scale_y_continuous(limits=c(0, 3.5))
ggsave("evidence_distribution_assembled.png")


# how many assemblies where we have no mapped evidence?
nrow(df[df$SCRM==0 & df$RPRM==0 & df$A_RM==0,]) # = number of calls in which *ALL* the breakpoint evidence has been removed -> no support for this variant yet we're calling them
# created breakpoint from assembly
ggplot(df[df$SCRM==0 & df$RPRM==0,], aes(x=A_MQT)) + geom_histogram()
# great assemblies in both directions
head(df[df$SCRM==0 & df$RPRM==0 & df$A_MQT>80,])

