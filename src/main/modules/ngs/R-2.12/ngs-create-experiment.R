# TOOL ngs-create-experiment.R: "Define NGS experiment" (This tool creates a combined count table for all the samples, and phenodata file which allows you to describe the samples and experiment setup. Please fill in the group column to describe the experiment groups to be compared.
# You have also the option to give the library size for each sample. If you fill in this information, the edgeR and DESeq tools will use these numbers rather than summing the counts from the count table.)
# INPUT sample{...}.tsv: "Individual sample files with counts" TYPE GENERIC 
# OUTPUT ngs-data-table.tsv: "Data table with read counts" 
# OUTPUT META phenodata.tsv: "Experiment description file" 
# PARAMETER experiment_type: "Type of experiment" TYPE [chip_seq: ChIP-seq, rna_seq: RNA-seq, mirna_seq: miRNA-seq] DEFAULT rna_seq (Experiment type.)
# PARAMETER alignment_type: "Does your data contain genomic coordinates" TYPE [genome: yes, other: no] DEFAULT other (Does your data table contain genomic coordinates.)
# PARAMETER impute_with: "Impute missing data" TYPE INTEGER FROM 0 TO 1000000 DEFAULT 0 (The value that is used for a sequence read that is not present in a sample.)
# PARAMETER count_column: "Count column" TYPE COLUMN_SEL DEFAULT EMPTY (Data column containing count data)

# MG 21.3.2011, takes as an input a set of files with read counts. Can also take read sequence, genomic location (chr, start, end) and length 
# MG, modified to deal with data that was aligned to other than a genome
# MK 30.4.2013, added ability to leave out genomic location data

# Sanity check the input data
files <- dir()
files <- files[grep("sample", files)]
number_files <- length(files)
if (number_files < 1) {
	stop("CHIPSTER-NOTE: You need to have at least 1 data file to run this tool!")
}

if (count_column == "EMPTY") {
	stop("CHIPSTER-NOTE: You must define the column containing count data!")
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

if("chr" %in% colnames(data_1)) {
	identifier_table <- data_1[, 1:4]
	for (count in 2:number_files) {
		this.identifier <- identifier_table[,1]
	}

	identifier_table <- identifier_table[!duplicated(identifier_table[,1]), ]	
	rownames(identifier_table) <- identifier_table[,1]
	if(length(identifiers) != nrow(identifier_table)) { stop("CHIPSTER-NOTE: An error occurred in table matching"); }
}

# Extract chromosome, start, end and length from id if aligned against genome
if (alignment_type == "genome") {
	if(!("chr" %in% colnames(data_1))) {
		extract_info <- strsplit(identifiers, "_")
		extract_info <- unlist(extract_info)
		chr_info <-  extract_info [seq(1,length(extract_info),4)]
		start_info <- extract_info [seq(2,length(extract_info),4)]
		end_info <- extract_info [seq(3,length(extract_info),4)]
		length_info <- as.numeric(end_info)-as.numeric(start_info)+1
		sequence_info <-  extract_info [seq(4,length(extract_info),4)]
	} else {
		chr_info <- identifier_table[identifiers, 2] 
		start_info <- identifier_table[identifiers, 3]
		end_info <- identifier_table[identifiers, 4]
		length_info <- as.numeric(identifier_table[identifiers, 4]) - as.numeric(identifier_table[identifiers, 3])
		sequence_info <-  rep(NA, length(identifiers))
	}
}

# Create table for all sample counts
annotation_columns <- 5
data_table <- array (dim = c(length(identifiers), number_files+annotation_columns))
data_table <- as.data.frame(data_table)
rownames (data_table) <- identifiers
colnames(data_table) [1:annotation_columns] <- c("chr","start","end","length","sequence")
if (alignment_type == "genome") {
	data_table$chr <- chr_info
	data_table$start <- start_info
	data_table$end <- end_info
	data_table$length <- length_info
	data_table$sequence <- sequence_info
}
if (alignment_type == "other") {
	data_table$chr <- "NA"
	data_table$start <- "NA"
	data_table$end <- "NA"
	data_table$length <- "NA"
	data_table$sequence <- "NA"
}

for (count in 1:number_files) {
	#print(count)
	# Fill in data where there is real data
	indices <- get (paste ("data_", count, sep=""))[,1]
	indices <- as.character(indices)
	#print(indices)
	#print(length(indices))
	if (alignment_type == "genome") {
		data_table[indices, (count+annotation_columns)] <- get (paste ("data_", count, sep=""))[get (paste ("data_", count, sep=""))[,1]==indices, count_column]
	}
	if (alignment_type == "other") {
		data_table[indices, (count+annotation_columns)] <- get (paste ("data_", count, sep=""))[get (paste ("data_", count, sep=""))[,1]==indices, count_column]
	}
	# Impute data where there isn't any
	indices_empty <- setdiff(identifiers, indices)
	#print(indices_empty)
	data_table[indices_empty, (count+annotation_columns)] <- impute_with
	colnames(data_table) [count+annotation_columns] <- paste("chip.", files[count], sep="")
	print (data_table [,count+annotation_columns])
}

# Generates the variables necessary to the phenodata file
# sample<-colnames(dat2)
groups <- c(rep("", number_files))
# Check whether column names contain "chip." and if so, remove from name
samples <- colnames (data_table) [(annotation_columns+1):(annotation_columns+number_files)]
if (length (grep("chip.", samples)) >= 1) {
	samples<-sub("chip.","",samples)
}

# Force identifiers to lower case if miRNA-seq data has been aligned against other than a reference genome
if (experiment_type == "mirna_seq" && alignment_type == "other") {
	rownames (data_table) <- tolower (rownames(data_table))
}

if(alignment_type == "other") {
	drops <- c("chr", "start", "end", "length", "sequence")
	data_table <- data_table[,!(names(data_table) %in% drops)]
}

# Writes out the data and the phenodata table
write.table(data_table, file="ngs-data-table.tsv", sep="\t", row.names=T, col.names=T, quote=F)
write.table(data.frame(sample=samples, chiptype="not applicable", experiment=experiment_type, group=groups, library_size=""), file="phenodata.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# EOF



