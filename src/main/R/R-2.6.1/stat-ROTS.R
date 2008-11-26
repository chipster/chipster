# ANALYSIS ROTS/"Two-group test" (Reproducibility-optimized test statistic for comparing the expression levels between two groups; Elo et al. 2008) 
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT ROTS.tsv, ROTSparameters.tsv
# PARAMETER column METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER fdr.threshold DECIMAL FROM 0 TO 1 DEFAULT 0.05 (FDR cut-off for significant results)
# PARAMETER B INTEGER FROM 10 TO 100000 DEFAULT 500 (Number of bootstrap and permutation resamplings)
# PARAMETER K INTEGER FROM 1000 TO 100000 DEFAULT 5000 (Largest top list size considered)



# Reproducibility-optimized two-group test

# LLE 23.9.2008



# Loads the ROTS library.

library(ROTS)



# Loads the normalized data

file <- c("normalized.tsv")
dat <- read.table(file, header=T, sep="\t", row.names=1)



# Separates expression values and flags

calls <- dat[,grep("flag", names(dat))]
data <- dat[,grep("chip", names(dat))]



# Test needs a parameter "groups" that specifies the grouping of the samples

phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
groups <- phenodata[,grep(column, colnames(phenodata))]



# Sanity checks

if(length(unique(groups))==1 | length(unique(groups))>=3) {
   stop("You need to have exactly two groups to run this analysis")
}



# Testing

result <- ROTS(data, groups, B, K) 
  
d <- result$d
FDR <- result$FDR
FDR[FDR > fdr.threshold] <- NA
write.table(na.omit(data.frame(dat, d=round(d, digits=2), FDR=round(FDR, digits=4))), file="ROTS.tsv", sep="\t", row.names=T, col.names=T, quote=F)
  
param<- round(c(result$a1, result$a2, result$k, result$R, result$Z, result$B),2)
names(param)<- c("a1","a2","Top list size","Reproducibility value","Z-score","Number of resamplings")
write.table(param, file="ROTSparameters.tsv", sep="\t", row.names=T, col.names=F, quote=F)







