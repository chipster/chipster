# ANALYSIS Normalisation/"Normalize to gene average" (Normalizes data so that the average signal for each gene becomes equal to the value of 1.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT normalized2geneaverage.tsv
# PARAMETER average.method [median, mean] DEFAULT mean (Method for calculating the average)

# Normalize the data so that each gene has an average of 1 across the
# chips of one data set
# MG 24.11.2009

# Name of the first table
name1<-c("normalized.tsv")

# Loads the table with data
table1<-read.table(file=name1, sep="\t", header=T, row.names=1)

# Get the indices for the column with expression values
indices <- grep("chip", names(table1))

# Calculate chip means and medians
gene.means <- apply(table1[,indices],FUN=mean,MARGIN=1)
gene.medians <- apply(table1[,indices],FUN=median,MARGIN=1)

# Calculating the normalization values
if (average.method=="mean") {
	for (row in 1:length(gene.means)) {
		table1[row,indices] <- table1[row,indices]/gene.means[row]
	}
}
if (average.method=="median") {
	for (row in 1:length(gene.means)) {
		table1[row,indices] <- table1[row,indices]/gene.medians[row]
	}
}

# Write out data
write.table(table1, file="normalized2geneaverage.tsv", sep="\t", row.names=T, col.names=T, quote=F)
