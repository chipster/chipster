# ANALYSIS "Normalisation"/"RMA (Affymetrix)" (RMA normalization for Affymetrix data.)
# INPUT AFFY microarray.cel OUTPUT rma.txt

# JTT, 12.7.2005

# loads the needed libraries
library(affy)

# Reads data
data<-ReadAffy()

# RMA normalization
data.rma<-rma(data)

# Getting expression values
eset.rma<-exprs(data.rma)

# Saving expression values to the data directory
write.exprs(data.rma, file="rma.txt", col.names=c('PROBESET_ID\tEXPRS'), quote=F)
