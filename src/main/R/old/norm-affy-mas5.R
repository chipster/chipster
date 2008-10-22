# ANALYSIS "Normalisation"/"MAS5 (Affymetrix)" (MAS5 preprocessing and scaling normalization for Affymetrix data.)
# INPUT AFFY microarray.cel OUTPUT mas5.txt

# JTT, 12.7.2005

# loads the needed libraries
library(affy)

# Reads data
data<-ReadAffy()

# MAS5 normalization
data.mas5<-mas5(data) # MAS5 analyses, produces non-normalized values, scale!
data.mas5.norm<-affy.scalevalue.exprSet(data.mas5, sc = 1000, analysis="absolute")

# Getting expression values
eset.mas5<-exprs(data.mas5.norm)

# Getting A/M/P calls
calls<-mas5calls(data)

# Saving normalized expression values and calls to the data directory
write.exprs(data.mas5.norm, file="mas5.txt", col.names=c('PROBESET_ID\tEXPRS'), quote=F)
write.exprs(calls, file="calls.txt")