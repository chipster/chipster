# TOOL bedtools-subtractbed.R: "Subtract BED" (Removes the portion(s\) of an interval that is overlapped by another feature(s\). This tool is based on the BEDTools package.)
# INPUT file.a: "Input file A" TYPE GENERIC
# INPUT file.b: "Input file B" TYPE GENERIC
# OUTPUT OPTIONAL subtractbed.bed 
# OUTPUT OPTIONAL error.txt
# PARAMETER OPTIONAL f: "Minimum overlap required as a fraction of A" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.000000001 (Minimum overlap required as a fraction of A. Default is 1E-9 (effectively 1bp\))
# PARAMETER OPTIONAL s: "Force strandedness" TYPE [yes, no] DEFAULT no (Only report hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.)

# AMS 23.4.2012
# AMS 23.9.2013 Improved outout/error file handling


# binary
binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "subtractBed"))

# options
options <- paste("")
options <- paste(options, "-f", f)
if (s == "yes") {options <- paste(options, "-s")}


# input files
options <- paste(options, "-a file.a -b file.b")

# command
command <- paste(binary, options, "> subtractbed.tmp 2> error.tmp")

# run
system(command)

# Generate output/error message
if (file.info("subtractbed.tmp")$size > 0) {
	system("mv subtractbed.tmp subtractbed.bed")
} else if (file.info("error.tmp")$size > 0) {
	system("mv error.tmp error.txt")
} else{
	system("echo \"# No results found\" > error.txt")
}
