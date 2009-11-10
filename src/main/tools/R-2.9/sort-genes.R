# ANALYSIS Utilities/"Sort genes" (Sort genes according to the selected column in ascending or descending order.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT sort-genes.tsv
# PARAMETER column COLUMN_SEL DEFAULT EMPTY (Data file column containing the data to sort)
# PARAMETER method [ascending, descending] DEFAULT ascending (Sort the genes from largest to smallest or vice versa)


# Sort genes
# JTT 29.1.2009
# Loads the libraries

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Extracts data to be sort
tosort<-dat[,grep(column, colnames(dat))]

# Sorting
if(method=="ascending") {
   dat<-dat[order(tosort, decreasing=F),]
} 
if(method=="descending") {
   dat<-dat[order(tosort, decreasing=T),]
} 

# Writing out data
write.table(data.frame(dat), file="sort-genes.tsv", sep="\t", row.names=T, col.names=T, quote=F)

