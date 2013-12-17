# TOOL extract-genes-from-stattest.R: "Extract genes using a p-value" (Extracts genes from a statistical test result. Specify the p-value column you want to use for extracting the genes.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT extract.tsv: extract.tsv 
# PARAMETER p.value.threshold: p.value.threshold TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off for significant results)
# PARAMETER p.value.column: p.value.column TYPE COLUMN_SEL (Column that contains the p-values)


# Extracts genes from the statistical test results
# JTT 8.11.2007

# Loads the data file
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1, check.names=F)

# Sanity checks
if(p.value.column=="empty") {
   stop("You haven't selected a data column for p-value! Tool cannot be executed.")
}
if(sum(grep("p.adjusted", colnames(dat)))==0) {
   stop("You don't have any P-value columns in the dataset! Please run some statistical test first.")
}

# Extracting the genes
dat2<-dat[which(dat[,grep(p.value.column, colnames(dat))]<=p.value.threshold),]

# Writing the data to disk
write.table(data.frame(dat2), file="extract.tsv", sep="\t", row.names=T, col.names=T, quote=F)

