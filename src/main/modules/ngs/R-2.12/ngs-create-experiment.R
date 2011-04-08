# TOOL ngs-create-experiment.R: "Define NGS experiment" (This tool creates a phenodata file containing descriptive information about samples and experiment setup. Please fill in, at least, the group column to describe the experiment groups to be compared for example in an analysis of differential gene expression.)
# INPUT counts_sample_{...}.tsv: "Individual sample files with reads and counts" TYPE GENERIC 
# OUTPUT ngs-data-table.tsv: "Data table with read counts"
# OUTPUT phenodata.tsv: "Experiment description file"
# PARAMETER experiment_type: "Type of experiment" TYPE [chip_seq: "ChIP-seq", rna_seq: "RNA-seq", mirna_seq: "miRNA-seq"] DEFAULT rna_seq (The value that is used for a sequence read that is not present in a sample.) 
# PARAMETER impute_with: "Impute missing data" TYPE INTEGER FROM 0 TO 1000000 DEFAULT 0 (The value that is used for a sequence read that is not present in a sample.) 

##############################
#                            #
# MG, 21.3.2011              #
#                            #
# Takes as an input a set of #
# of files with read counts. #
# Additional information     #
# includes:                  #
#                            #
# read sequence              #
# chromosome name            #
# start position             #
# end position               #
# read sequence length       #
#                            #
##############################

# Sanity check the input data
files<-dir()
files<-files[grep("counts_sample_", files)]
number_files <- length(files)
if (number_files < 1) {
	stop("CHIPSTER-NOTE: You need to have 1 data file to run this tool!")
}

# Read in data files
for (count in 1:number_files) {
	assign (paste("data_", count, sep=""), read.table(files[count], header=T, sep="\t")) 
}

# Merge into single data table
# When there are no count data for some samples for a given sequence read impute the number requested by the
# impute_with parameter

# Make a union of all identifiers
identifiers <- data_1[,1]
for (count in 2:number_files) {
	identifiers <- (union (identifiers, get (paste ("data_", count, sep=""))[,1]))
}
identifiers <- as.character(unique(identifiers))

# Create table for all sample counts
annotation_columns <- 5
data_table <- array (dim = c(length(identifiers), number_files+annotation_columns))
data_table <- as.data.frame(data_table)
rownames (data_table) <- identifiers
colnames(data_table) [1:annotation_columns] <- c("chr","start","end","length","sequence")
for (count in 1:number_files) {
	print(count)
	# Fill in data where there is real data
	indices <- get (paste ("data_", count, sep=""))[,1]
	indices <- as.character(indices)
	print(indices)
	print(length(indices))
	data_table[indices, (count+annotation_columns)] <- get (paste ("data_", count, sep=""))[get (paste ("data_", count, sep=""))[,1]==indices,3]
	# Impute data where there isn't any
	indices_empty <- setdiff(identifiers, indices)
	print(indices_empty)
	data_table[indices_empty, (count+annotation_columns)] <- impute_with
	colnames(data_table) [count+annotation_columns] <- paste("Sample_", count, sep="")
	print (data_table [,count+annotation_columns])
}


# Generates the variables necessary to the phenodata file
# sample<-colnames(dat2)
groups <- c(rep("", number_files))
samples <- colnames (data_table) [(annotation_columns+1):(annotation_columns+number_files)]

# Writes out the data and the phenodata table
write.table(data_table, file="ngs-data-table.tsv", sep="\t", row.names=T, col.names=T, quote=F)
write.table(data.frame(sample=samples, chiptype="not applicable", experiment=experiment_type, group=groups), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# EOF



