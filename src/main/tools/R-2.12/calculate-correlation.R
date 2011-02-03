# ANALYSIS Utilities/"Calculate sample correlations" (Calculates correlations between samples.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT cor.tsv
# PARAMETER correlation.method [pearson, spearman] DEFAULT pearson (Correlation measure)


# Correlogram
# JTT 18.10.2007
# jeejeejeeJeu

# Parameter settings (default) for testing purposes
#correlation.method<-c("pearson")

# JOOO
 
# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Writes a table of correlation values
write.table(as.data.frame(cor(dat2, method=correlation.method)), "cor.tsv", sep="\t", quote=F)
