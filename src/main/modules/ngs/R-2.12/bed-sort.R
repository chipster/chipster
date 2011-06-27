# TOOL bed-sort.R: "Sort BED files" (Sort BED files.)
# INPUT regions.bed TYPE GENERIC 
# OUTPUT sorted.bed 



# EK 22.6.2011

# binary
binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "sortBed"))

# command
command <- paste(binary, "-sizeD -i regions.bed > sorted.bed")

# run
system(command)


