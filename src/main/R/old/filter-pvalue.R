# ANALYSIS "Expression"/"P-Value Filtering" (P-value filtering)
# INPUT GENE_EXPRS p-anova.txt OUTPUT p-anova-filter.txt
# PARAMETER adj [rawp, Bonferroni, Holm, Hochberg, SidakSS, SidakSD, BH, BY] DEFAULT BH (which p-values are used for filtering)

# JTT, 13.7.2005

# Loads the normalized data
# This needs a parameter "file" that tells which data file should be loaded.
# Datafile is produced by some statistical tool (currently, only stat-anova.R)
# At the moment datafile can be only be "p-anova.txt"
file<-c("p-anova.txt")
p<-read.table(file, header=T, sep="\t")

# Filter
# This needs a parameter that tells which p-values are used for filtering.
# Options are "rawp","Bonferroni","Holm","Hochberg","SidakSS","SidakSD","BH","BY"             

#adj<-c("BH")
rawp<-paste("p",adj,sep="$")
#rawp<-p$BH
genenames<-p$row.names
pfilter<-ifelse(rawp<=0.05, rawp, NA)
p2<-data.frame(genenames,pfilter)
p3<-na.omit(p2)

# Saving the results
write.table(p3, "p-anova-filter.txt", sep="\t", row.names=F, col.names=T)