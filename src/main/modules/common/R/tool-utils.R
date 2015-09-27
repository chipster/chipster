read_input_definitions <- function(){
	# Read in the data
	list <- scan("chipster-inputs.tsv", what="", sep="\n", comment.char="#")
	# Separate elements by one or more whitepace
	inputdef <- strsplit(list, "[[:space:]]+")
	# Extract the first vector element and set it as the list element name
	names(inputdef) <- sapply(inputdef, function(list) list[[1]])
	# Remove the first vector element from each list element
	inputdef <- lapply(inputdef, function(list) list[-1])
	
	return(inputdef)
}							

write_output_definitions <- function(output_names){
	write.table(output_names, file = "chipster-outputs.tsv", row.names = F, col.names= F, quote = F, sep = "\t")
}

# Removes user-specified postfix if present
#
remove_postfix <- function(name, postfix){
	pf <- paste(postfix, "$", sep="")
	if (grepl(pf, name)){
		basename <- substr(name, 1, (nchar(name) - nchar(postfix)))
		return(basename)
	}else{
		return(name)
	}
}

# Removes extension, i.e. everything after the last dot (including the dot)
#
remove_extension <- function(name){
	return(sub("\\.[^.]*$", "", name))
}

# Strips common file extensions from a file name
#
strip_name <- function(name){
	known_postfixes <- c(".gz", ".bam", ".sam", ".fa", ".fasta", ".fq", ".fastq", ".gtf", ".tsv", "_trimmed", "_filtered")
	newname <- name
	while (TRUE){
		for (i in known_postfixes){
			newname <- remove_postfix(newname, i)
		}
		if (nchar(newname) == nchar(name)){
			break
		}
		name <- newname
	}
	return(name)
	
}

# If the names look like typical paired-end names: *_1, *_2, remove the ending and return the name.
# If not, return first name as-is.
#
paired_name <- function(name1, name2){
	if (grepl("_1$", name1) && grepl("_2$", name2)){
			return(remove_postfix(name1, "_1"))
	}
	if (grepl("_2$", name1) && grepl("_1$", name2)){
			return(remove_postfix(name1, "_2"))
	}
	return(name1)
}