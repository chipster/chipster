# TOOL bedtools-overlap.R: "Overlap BED" (Computes the amount of overlap (positive values\) or distance (negative values\) between genome features and reports the result at the end of the same line. This tool is based on the BEDTools package.)
# INPUT file.a: "BED file A" TYPE GENERIC
# OUTPUT OPTIONAL overlap.bed
# OUTPUT OPTIONAL error.txt
# PARAMETER cols: "Columns to use" TYPE STRING DEFAULT "2,3,6,7" (Specify the columns (1-based\) for the starts and ends of the features for which you'd like to compute the overlap/distance. The columns must be listed in the following order: start1,end1,start2,end2.)

# AMS 23.4.2012
# AMS 23.9.2013 Improved outout/error file handling

# binary
# Note: overlap binary for now from from version 2.12.0 as it is not present in 2.17.0
binary <- c(file.path(chipster.tools.path, "BEDTools-Version-2.12.0", "bin", "overlap"))

# optional options
options <- paste("")
options <- paste(options, "-cols", cols)

# input files
options <- paste(options,"-i file.a")

# command
command <- paste(binary, options, " > overlap.tmp 2> error.tmp")

# run
system(command)

# Generate output/error message
if (file.info("overlap.tmp")$size > 0) {
	system("mv overlap.tmp overlap.bed")
} else if (file.info("error.tmp")$size > 0) {
	system("mv error.tmp error.txt")
} else{
	system("echo \"# No results found\" > error.txt")
}
