# Utility function for adding a vector of string into same length (maximum length found from the vector).
pad.with.zeroes <- function(str.vector) {
	len <- max(nchar(str.vector))
	str.vector.2 <- paste(paste(rep("0", len), collapse=""), str.vector, sep="")
	return( substr(str.vector.2, nchar(str.vector.2)-(len-1), nchar(str.vector.2)) )
}


# First sort by chromosome, then by start coordinate.
# Chromosomes that are numeric are compated numerically. Chromosomes that
# are non-numeric are compared lexically (in their normalised form).
# Numeric names are always considered smaller than non-numeric.
# This is the same sorting principle that Chipster genome browser has, hence
# using this function guarantees that BED files are compatible with it.
#
# R does not have flexible sorting functionalities, so we need to resort to a hackish solution.
# We create alphabetic presentation of the values and use that for sorting.    
#
# Usage:
#  bed <- read.table(file="sortme.bed", skip=1, sep="\t") # assume file has 1 line header
#  colnames(bed)[1:2] <- c("chr", "start")  # these named columns are required for sorting 
#  sorted.bed <- sort.bed(bed)
#  write.table(sorted.bed, file="sorted.bed", sep="\t", row.names=F, col.names=F, quote=F)
# 
sort.bed <- function(bed) {

	# Check inputs
	if (!'chr' %in% colnames(bed)) {
		stop('BED is missing column chr');
	}
	if (!'start' %in% colnames(bed)) {
		stop('BED is missing column start');
	}
	
	# Normalise chromosome names
	chr.without.postfix <- gsub('(.*)\\..*', '\\1', bed$chr)
	chr.normalised <- gsub('chr(.*)', '\\1', chr.without.postfix)

	# Convert all fields to right format
	chr.is.nonnumeric <- ifelse(is.na(as.numeric(chr.normalised)), '1', '0')
	chr.padded <- pad.with.zeroes(chr.normalised)
	start.padded <- pad.with.zeroes(bed$start)
	
	# Create the strings to be sorted
	string.matrix <- rbind(chr.is.nonnumeric, chr.padded, start.padded)
	string.vector <- apply(string.matrix, 2, function(row) paste(row[1], row[2], row[3], sep=""))
	
	# Sort the bed, using the generated strings
	so <- order(string.vector)
	sorted.bed <- bed[so,]
	
	return(sorted.bed)
}

