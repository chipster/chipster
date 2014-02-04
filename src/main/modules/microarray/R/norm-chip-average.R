# TOOL norm-chip-average.R: "Normalize to chip average" (Normalizes data so that the average signal of each chip becomes equal to the value of 1.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT normalized2chipaverage.tsv: normalized2chipaverage.tsv 
# PARAMETER average.method: average.method TYPE [median: median, mean: mean] DEFAULT mean (Method for calculating the average)

# Normalize the data on each chip to the average of the chip
# MG 23.11.2009

# Name of the first table
name1<-c("normalized.tsv")

# Loads the table with data
table1<-read.table(file=name1, sep="\t", header=T, row.names=1)

# Get the indices for the column with expression values
indices <- grep("chip", names(table1))

# Calculate chip means and medians
chip.means <- apply(table1[,indices],FUN=mean,MARGIN=2)
chip.medians <- apply(table1[,indices],FUN=median,MARGIN=2)

# Calculating the normalization values
if (average.method=="mean") {
	for (column in 1:length(indices)) {
		table1[,indices[column]] <- table1[,indices[column]]/chip.means[column]
	}
}
if (average.method=="median") {
	for (column in 1:length(indices)) {
		table1[,indices[column]] <- table1[,indices[column]]/chip.medians[column]
	}
}

# Write out data
write.table(table1, file="normalized2chipaverage.tsv", sep="\t", row.names=T, col.names=T, quote=F)
