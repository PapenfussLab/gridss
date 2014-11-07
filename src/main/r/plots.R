library(ggplot2)

# local anchor length shouldn't mean much
ggplot(pp, aes(x=A_BLRM, y=A_BLLM)) + geom_point()

# should have a linear relationship
ggplot(pp, aes(x=log(RC), y=log(PC))) + geom_point() + geom_jitter()

ggplot(pp, aes(x=A_SC + A_RP, y=RC, color=PC)) + geom_point()

# assembly length vs breakpoint quality
ggplot(pp, aes(x=A_BLRM, y=A_BLLM, color=A_BLRM)) + geom_point()
+ scale_y_log10() + scale_colour_gradientn(colours=rainbow(4)) 

# normal coverage
ggplot(pp, aes(x=RCNormal + SCECNormal)) + geom_histogram(binwidth=1)
#total support
ggplot(pp, aes(x=A_SC)) + geom_histogram(binwidth=1)
ggplot(pp, aes(x=A_RP)) + geom_histogram(binwidth=1)
ggplot(pp, aes(x=PC)) + geom_histogram(binwidth=1)
ggplot(pp, aes(x=RC)) + geom_histogram(binwidth=1)
ggplot(pp, aes(x=A_SCCLM)) + geom_histogram(binwidth=1)
ggplot(pp, aes(x=A_MQRM)) + geom_histogram(binwidth=1)
ggplot(pp, aes(x=A_MQRT)) + geom_histogram(binwidth=1)
ggplot(pp, aes(x=A_EC)) + geom_histogram(binwidth=1)
ggplot(pp, aes(x=A_SCTumour / (A_SCTumour + RCTumour))) + geom_histogram(binwidth=0.05)

# distribution of tumour variant evidence contribution
ggplot(df, aes(x=(A_SCTumour + A_RPTumour) / (A_SC + A_RP))) + geom_histogram(binwidth=0.01)
# distribution of variant evidence contribution
ggplot(df, aes(x=A_SC+A_RP, y=RC+PC, color=(A_SC+A_RP)/(A_SC+A_RP+RC+PC))) + geom_point() + scale_y_log10() + scale_x_log10()
# normal coverage
ggplot(df, aes(x=RCNormal + SCECNormal)) + geom_histogram(binwidth=1)