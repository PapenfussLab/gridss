library(ggplot2)
pp <- df[
  (df$A_SCNormal + df$A_RPNormal) / (df$A_SC + df$A_RP) < 0.05 & # less than 5% of support is from the normal # grep "A_SC=0," CPCG0100.vcf | grep "A_RP=0," 
    df$RCNormal >= 6 & # min coverage in the normal supporting the reference allele
    df$A_SCTumour / (df$A_SCTumour + df$RCTumour) >= 0.20 & # at least 20% support for variant in the tumour data
    df$A_RM >=2 & # assembly support from both sides
#    df$A_SC >= 2 & # found an exact breakpoint on both sides
    df$A_SC + df$A_RP/2 >= 4 & # at least this many read pairs support the variant
    df$A_MQRT >= 30 & # reasonable mapping
    1==1,]
# local anchor length shouldn't mean much
ggplot(pp, aes(x=A_BLRM, y=A_BLLM)) + geom_point()

# should have a linear relationship
ggplot(pp, aes(x=log(RC), y=log(PC))) + geom_point() + geom_jitter()

ggplot(pp, aes(x=A_SC + A_RP, y=RC, color=PC)) + geom_point()

# assembly length vs breakpoint quality
ggplot(pp, aes(x=A_BLRM, y=A_BLLM, color=A_BLRM)) + geom_point()
+ scale_y_log10() + scale_colour_gradientn(colours=rainbow(4)) 

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
ggplot(pp, aes(x=A_SCTumour / (A_SCTumour + RCTumour))) + geom_histogram(binwidth=0.05)

# distribution of tumour variant evidence contribution
ggplot(df, aes(x=(A_SCTumour + A_RPTumour) / (A_SC + A_RP))) + geom_histogram(binwidth=0.01)
# distribution of variant evidence contribution
ggplot(df, aes(x=A_SC+A_RP, y=RC+PC, color=(A_SC+A_RP)/(A_SC+A_RP+RC+PC))) + geom_point() + scale_y_log10() + scale_x_log10()
# normal coverage
ggplot(df, aes(x=RCNormal + SCECNormal)) + geom_histogram()

