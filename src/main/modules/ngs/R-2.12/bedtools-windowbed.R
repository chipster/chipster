# TOOL bedtools-windowbed.R: "BEDTools windowBed" (Examines a window around each feature in A and reports all features in B that overlap the window. For each overlap the entire entry in A and B are reported.)
# INPUT file.a: "Input file A" TYPE GENERIC
# INPUT file.b: "Input file B" TYPE GENERIC
# OUTPUT windowbed.bed 
# PARAMETER OPTIONAL abam: "The A input file is in BAM format" TYPE [yes,no] DEFAULT no (The A input file is in BAM format. By default output will be BAM as well.) 
# PARAMETER OPTIONAL ubam: "Write uncompressed BAM output" TYPE [yes,no] DEFAULT no (Write uncompressed BAM output. Default is to write compressed BAM.) 
# PARAMETER OPTIONAL bed: "When using BAM input, write output as BED." TYPE [yes,no] DEFAULT no (When using BAM input, write output as BED. The default is to write output in BAM.) 
# PARAMETER OPTIONAL l: "Base pairs added upstream" TYPE INTEGER DEFAULT 1000 (Base pairs added upstream (left of\) of each entry in A when searching for overlaps in B. Allows one to define assymterical windows. Default is 1000 bp.) 
# PARAMETER OPTIONAL r: "Base pairs added downstream" TYPE INTEGER DEFAULT 1000 (Base pairs added downstream (right of\) of each entry in A when searching for overlaps in B. Allows one to define assymterical windows. Default is 1000 bp.) 
# PARAMETER OPTIONAL sw: "Define window based on strand" TYPE [yes,no] DEFAULT no (Define window based on strand.  For example if used, upstream window size 500 for a negative-stranded feature will add 500 bp downstream.) 
# PARAMETER OPTIONAL sm: "Only report hits in B that overlap A on the same strand" TYPE [yes,no] DEFAULT no (Only report hits in B that overlap A on the same strand.) 
# PARAMETER OPTIONAL u: "Write the original A entry once if any overlaps found in B" TYPE [yes,no] DEFAULT no (Write the original A entry once if any overlaps found in B. In other words, just report the fact >=1 hit was found.) 
# PARAMETER OPTIONAL c: "For each entry in A, report the number of overlaps with B" TYPE [yes,no] DEFAULT no (For each entry in A, report the number of overlaps with B. Reports 0 for A entries that have no overlap with B.) 
# PARAMETER OPTIONAL v: "Only report those entries in A that have no overlaps with B" TYPE [yes,no] DEFAULT no (Only report those entries in A that have _no overlaps_ with B) 

# binary
binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "windowBed"))

# optional options
options <- paste("")
if (abam == "yes") {
	options <- paste(options, "-abam")}
	if (ubam == "yes") {options <- paste(options, "-ubam")}
	if (bed == "yes") {options <- paste(options, "-bed")}
}
if (sw == "yes") {options <- paste(options, "-sw")}
if (sm == "yes") {options <- paste(options, "-sm")}
if (u == "yes") {options <- paste(options, "-u")}
if (c == "yes") {options <- paste(options, "-c")}
if (v == "yes") {options <- paste(options, "-v")}

# common options
options <- paste(options, "-l", l, "-r", r, "-a file.a -b file.b")

# command
command <- paste(binary, options, " > windowbed.bed")

# run
system(command)
