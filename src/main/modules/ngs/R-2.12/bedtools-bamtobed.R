# TOOL bedtools-bamtobed.R: "Convert BAM to BED" (Converts BAM alignments to BED6, BED12 or BEDPE format. This tool is based on the BEDTools package.)
# INPUT file.a: "BAM file" TYPE GENERIC
# OUTPUT OPTIONAL bamtobed.bed
# OUTPUT OPTIONAL error.txt
# PARAMETER output.type: "Output type" TYPE [BED6, BED12, BEDPE] DEFAULT "BED6" (Select the output type (BED6, BED12 or BEDPE\).)
# PARAMETER OPTIONAL split: "Report split BAM alignments as separate BED entries" TYPE [yes,no] DEFAULT no (Report each portion of a split BAM (i.e., having an N CIGAR operation\) alignment as a distinct BED intervals.)
# PARAMETER OPTIONAL ed: "Use BAM edit distance for BED score" TYPE [yes,no] DEFAULT no (Use BAM edit distance (NM tag\) for BED score. Default for BED is to use mapping quality. Default for BEDPE is to use the minimum of the two mapping qualities for the pair.)
# PARAMETER OPTIONAL color: "Color string" TYPE STRING DEFAULT "255,0,0" (An R,G,B string for the color used with BED12 format. Default is (255,0,0\).)
# PARAMETER OPTIONAL cigar: "Add the CIGAR string" TYPE [yes,no] DEFAULT no (Add the CIGAR string to the BED entry as a 7th column.)

# AMS 23.4.2012
# AMS 23.9.2013 Improved outout/error file handling

# binary
binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "bamToBed"))

# options
options <- paste("")
if (output.type == "BEDPE") {options <- paste(optons, "-bedpe")}
if (output.type == "BED12") {
	options <- paste(options, "-bed12")
	options <- paste(options, "-color", color)
}
if (split == "yes") {options <- paste(options, "-split")}
if (ed == "yes") {options <- paste(options, "-ed")}
if (cigar == 'yes') {options <- paste(options,"-cigar")}

# input files
options <- paste(options, "-i", "file.a")

# command
command <- paste(binary, options, "> bamtobed.tmp 2> error.tmp")

# run
system(command)

# Generate output/error message
if (file.info("bamtobed.tmp")$size > 0) {
	system("mv bamtobed.tmp bamtobed.bed")
} else if (file.info("error.tmp")$size > 0) {
	system("mv error.tmp error.txt")
} else{
	system("echo \"# No results found\" > error.txt")
}
