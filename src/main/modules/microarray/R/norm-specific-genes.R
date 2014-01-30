# TOOL norm-specific-genes.R: "Normalize to specific genes" (Normalizes data to specific genes.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT normalized-too.tsv: normalized-too.tsv TYPE GENE_EXPRS 
# OUTPUT normalized2genes.tsv: normalized2genes.tsv 


# Normalize the data to specific genes
# JTT 29.1.2009

# Name of the first table
name1<-c("normalized.tsv")

# Name of the second table
name2<-c("normalized-too.tsv")

# Loads the tables
table1<-read.table(file=name1, sep="\t", header=T, row.names=1)
table2<-read.table(file=name2, sep="\t", header=T, row.names=1)

# Table2 should always hold the genelist specifying the control probes
if(nrow(table2)>nrow(table1)) {
   table3<-table1
   table1<-table2
   table2<-table3
}

# Extract the control genes/probes from the whole dataset
table1c<-table1[rownames(table2),]

# Separates expression values
table1cd<-table1c[,grep("chip", names(table1c))]

# Calculating the normalization values
means<-colMeans(table1cd)

# Performing the normalization
cols<-grep("chip", names(table1c))
for(i in 1:length(cols)) {
   table1[,cols[i]]<-table1[,cols[i]]-means[i]
}

# Write out data
write.table(table1, file="normalized2genes.tsv", sep="\t", row.names=T, col.names=T, quote=F)
