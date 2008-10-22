# ANALYSIS "Expression"/"ANOVA" (t-test and ANOVA for differencies between groups.)
# INPUT GENE_EXPRS normalised.txt OUTPUT p-anova.txt

# JTT, 13.7.2005

# Load the needed packages
library(multtest)

# Loads the normalized data
# This needs a parameter "file" that tells which data file should be loaded.
# Datafiles are produced by normalization scripts, and at the moment
# datafile can be any of "rma.txt", "mas5.txt", or "liwong.txt".
file<-c("normalised.txt")
data<-read.table(file, header=T, sep="\t")

# In order to be able to calculate the test, data should be transposed
transposed.data<-t(data)

# Tests need a parameter "groups" that specifies the grouping of the samples
groups<-c(1,1,2,2,3,3)
# Vector can also be read from a file
# groups<-scan("groups", sep=",")

# Calculation of test statistics
# For Affymetrix arrays with 23 000 genes this might take about 15 minutes
# The calculation of test statistic takes only about a minute,
# loop itself take much longer, about 3 minutes and 30 seconds.
p<-0
len<-length(data[,1])
for(i in 1:len) {
p<-c(p, na.omit((anova(lm(transposed.data[,i]~groups)))$Pr)[1])
}
p.raw<-p[-1]


# This is just for testing speed...
# aaa<-date()
# p<-0
# test<-function(x) {
# p<-c(p, na.omit(anova(lm(transposed.data[,1:10]~groups))$Pr)[1])
# na.omit(summary(aov(transposed.data[,1:len]~groups))[[1]]$Pr)[1]
# }
# p<-c(p, apply(transposed.data, MARGIN=1, FUN=test))
# aaaa<-date()
# aaa;aaaa;
# 
# aaa<-date()
# aa<-summary(aov(transposed.data[,1:len]~groups))
# aaaa<-date()
# aaa;aaaa;
# Speed testing ends here


# Adjusting the p-values
procs<-c("Bonferroni","Holm","Hochberg","SidakSS","SidakSD","BH","BY")
p.adjusted<-mt.rawp2adjp(p.raw,procs)

# Saving the p-values with Affymetrix IDs
p.saved<-data.frame(row.names(data), p.adjusted$adjp[order(p.adjusted$index),])
write.table(p.saved, "p-anova.txt", sep="\t", row.names=F, col.names=T)


