library(ggplot2)
library(RColorBrewer)
source("libgridss.R")

vcf <- readVcf("W:/778/idsv/778.in.vcf", "hg19_random")
df <- gridss.truthdetails.processvcf.vcftodf(vcf)

pp <- df

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

# Effect of local breakend evidence
ggplot(df, aes(x=QUAL, y=CQUAL)) + geom_point() + scale_x_log10() + scale_y_log10()

