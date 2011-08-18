# TOOL bedtools-overlap.R: "BEDTools overlap" (Computes the amount of overlap (positive values\) or distance (negative values\) between genome features and reports the result at the end of the same line.)
# INPUT file.a: "BED file A" TYPE GENERIC
# OUTPUT overlap.bed 
# PARAMETER cols: "Columns to use" TYPE STRING DEFAULT "2,3,6,7" (Specify the columns (1-based\) for the starts and ends of the features for which you'd like to compute the overlap/distance. The columns must be listed in the following order: start1,end1,start2,end2.)


# binary
binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "overlap"))

# optional options
options <- paste("")
options <- paste(options, "-cols", cols)

# input files
options <- paste(options,"-i file.a")

# command
command <- paste(binary, options, " > overlap.bed")

# run
system(command)
if (file.info("overlap.bed")$size == 0) {system("echo \"No results found\" > overlap.bed")}