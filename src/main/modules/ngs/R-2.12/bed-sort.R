# TOOL bed-sort.R: "Sort BED files" (Sorts a BED file by chromosome and then by start position. This tool is based on the BEDTools package.)
# INPUT regions.bed TYPE GENERIC 
# OUTPUT sorted.bed 



# EK 22.6.2011

# binary
binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "sortBed"))

# command
command <- paste(binary, "-i regions.bed > sorted.bed")

# run
system(command)


