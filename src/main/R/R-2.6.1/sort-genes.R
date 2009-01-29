# ANALYSIS Statistics/"Sort genes" (Sort genes according to the selected column in ascending or descending order.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT sort-genes.tsv
# PARAMETER column COLUMN_SEL DEFAULT EMPTY (Data file column containing the data to sort)
# PARAMETER method [ascending, descending] DEFAULT ascending (Sort the genes from largest to smallest or vice versa)


# Sort genes
# JTT 29.1.2009
# Loads the libraries

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Extracts raw p-values
groups<-dat[,grep(column, colnames(dat))]

# Sorting
if(method="ascending") {
   dat2<-dat2[,order(groups, decreasing=F)]
} 
if(method="descending") {
   dat2<-dat2[,order(groups, decreasing=T)]
} 

# Writing out data
write.table(data.frame(dat, adjp2df), file="sort-genes.tsv", sep="\t", row.names=T, col.names=T, quote=F)

