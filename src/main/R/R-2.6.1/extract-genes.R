# ANALYSIS Utilities/"Extract genes" (Extracts a specified number of genes from the top of the data.) 
# INPUT GENE_EXPRS normalized.tsv OUTPUT extract-genes.tsv
# PARAMETER number.of.genes INTEGER FROM 1 TO 1000000 DEFAULT 50 (Number of genes to extract)


# Extracts genes 
# JTT 29.1.2009

# Loads the data file
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Extracting the genes
dat2<-dat[number.of.genes,]

# Writing the data to disk
write.table(data.frame(dat2), file="extract-genes.tsv", sep="\t", row.names=T, col.names=T, quote=F)

