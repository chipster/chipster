# First sort by chromosome, then by start coordinate.
# Chromosomes that are numeric are compated numerically. Chromosomes that
# are non-numeric are compared lexically (in their normalised form).
# Numeric names are always considered smaller than non-numeric.
# This is the same sorting principle that Chipster genome browser has, hence
# using this function guarantees that sorted files are compatible with it.    
#
sort.gtf <- function(input, output) {
	system(paste("java -cp  '", chipster.module.path, "/../../../shared/lib/*' fi.csc.chipster.tools.ngs.SortGtf ", input, " ", output, sep = ""))
}

# Adds "chr" to the beginning of each line that starts with a number or with X, Y, Z, W or M
#
addChrToGtf <- function(input, output){
	awk.command <- paste("awk '{ if ($1 ~ /^[0-9,X,Y,Z,W,M]/) {print \"chr\"$0} else {print $0} }'", input, ">", output)
	system(awk.command)
}