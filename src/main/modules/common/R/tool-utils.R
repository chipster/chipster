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
	#if (nchar(name) > 0 && nchar(postfix) > 0 ){
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
	known_postfixes <- c(".gz", ".bam", ".fa", ".fasta", ".fq", ".fastq", ".gtf", ".tsv", "_trimmed", "_filtered")
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