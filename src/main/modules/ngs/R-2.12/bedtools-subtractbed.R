# TOOL bedtools-subtractbed.R: "BEDTools subtractBed" (Removes the portion(s\) of an interval that is overlapped by another feature(s\).)
# INPUT file.a: "Input file A" TYPE GENERIC
# INPUT file.b: "Input file B" TYPE GENERIC
# OUTPUT subtractbed.bed 
# PARAMETER OPTIONAL f: "Minimum overlap required as a fraction of A" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.000000001 (Minimum overlap required as a fraction of A. Default is 1E-9 (effectively 1bp\))
# PARAMETER OPTIONAL s: "Force strandedness" TYPE [yes, no] DEFAULT no (Only report hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.)

# binary
binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "subtractBed"))

# options
options <- paste("")
options <- paste(options, "-f", f)
if (s == "yes") {options <- paste(options, "-s")}


# input files
options <- paste(options, "-a file.a -b file.b")

# command
command <- paste(binary, options, "> subtractbed.bed")

# run
system(command)
if (file.info("subtractbed.bed")$size == 0) {system("echo \"No results found\" > subtractbed.bed")}