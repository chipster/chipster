# ANALYSIS "Normalisation"/"Li-Wong (Affymetrix)" (Li-Wong normalization for Affymetrix data. Uses the same algorithm as dChip, but does not give the same results.)
# INPUT AFFY microarray.cel OUTPUT liwong.txt

# JTT, 12.7.2005

# loads the needed libraries
library(affy)

# Reads data
data<-ReadAffy()

# Li-Wong normalization method
data.li.wong<-data.liwong <- expresso(data, normalize.method="invariantset", bg.correct=FALSE, pmcorrect.method="pmonly",summary.method="liwong")

# Getting expression values
eset.li.wong<-exprs(data.li.wong)

# Saving expression values to the data directory
write.exprs(data.li.wong, file="liwong.txt", col.names=c('PROBESET_ID\tEXPRS'), quote=F)
