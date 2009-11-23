# ANALYSIS Normalisation/"Normalize to chip average" (Normalizes data so that the average signal of each chip becomes equal to the value of 1.)
# INPUT GENE_EXPRS normalized.tsv, GENE_EXPRS normalized-too.tsv OUTPUT normalized2genes.tsv
# PARAMETER average.method [median, mean] DEFAULT mean (Method for calculating the average)

# Normalize the data on each chip to the average of the chip
# MG 23.11.2009

# Name of the first table
name1<-c("normalized.tsv")

# Name of the second table
name2<-c("normalized-too.tsv")

# Loads the tables
table1<-read.table(file=name1, sep="\t", header=T, row.names=1)
#table2<-read.table(file=name2, sep="\t", header=T, row.names=1)

# Table2 should always hold the genelist specifying the control probes
#if(nrow(table2)>nrow(table1)) {
#   table3<-table1
#   table1<-table2
#   table2<-table3
#}

# Extract the control genes/probes from the whole dataset
#table1c<-table1[rownames(table1),]

# Separates expression values
table1cd<-table1[,grep("chip", names(table1c))]

# Calculating the normalization values
means<-colMeans(table1cd)

# Performing the normalization
cols<-grep("chip", names(table1c))
for(i in 1:length(cols)) {
   table1[,cols[i]]<-table1[,cols[i]]-means[i]
}

# Write out data
write.table(table1, file="normalized2genes.tsv", sep="\t", row.names=T, col.names=T, quote=F)
