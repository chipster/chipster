# TOOL bedtools-pairtopair.R: "BEDTools pairToPair" (pairToPair compares two BEDPE files in search of overlaps where each end of a BEDPE feature in A overlaps with the ends of a feature in B.)
# INPUT file.a: "BEDPE file A" TYPE GENERIC
# INPUT file.b: "BEDPE file B" TYPE GENERIC
# OUTPUT pairtopair.bed 
# PARAMETER OPTIONAL f: "Minimum overlap" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.000000001 (Minimum overlap required as a fraction of A. Default is 1E-9 (effectively 1bp\))
# PARAMETER OPTIONAL type: "Approach to reporting overlaps" TYPE [either, neither, both] DEFAULT both (either: Report overlaps if either ends of A overlap B. neither: Report A if neither end of A overlaps B. both: Report overlaps if both ends of A overlap B.)
# PARAMETER OPTIONAL is: "Ignore strands when searching for overlaps" TYPE [yes, no] DEFAULT no (Ignore strands when searching for overlaps. By default, strands are enforced.)
# PARAMETER OPTIONAL rdn: "Require different names" TYPE [yes, no] DEFAULT no (Require the hits to have different names (i.e. avoid self-hits\). By default, same names are allowed.)
# PARAMETER OPTIONAL addslop: "Add slop" TYPE [yes, no] DEFAULT no (Note: Slop is subtracted from start1 and start2 and added to end1 and end2.)
# PARAMETER OPTIONAL slop: "The amount of slop in bp to be added to each footprint" TYPE INTEGER DEFAULT 0 (The amount of slop in bp to be added to each footprint.)
# PARAMETER OPTIONAL ss: "Add slop based to each BEDPE footprint based on strand" TYPE [yes, no] DEFAULT no (Add slop based to each BEDPE footprint based on strand. If strand is +, slop is only added to the end coordinates. If strand is -, slop is only added to the start coordinates. By default, slop is added in both directions.) 

# binary
binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "pairToPair"))

# optional options
options <- paste("")
if (is == "yes") {options <- paste(options, "-si")}
if (rdn == "yes") {options <- paste(options, "-rdn")}
options <- paste(options, "-f", f) 
options <- paste(options,"-type", type)
if (addslop == "yes"){
	options <- paste(options, "-slop", slop)
	if (ss == "yes") {options <- paste(options, "-ss")}
}
# input files
options <- paste(options,"-a file.a -b file.b")

# command
command <- paste(binary, options, " > pairtopair.bed")

# run
system(command)
if (file.info("pairtopair.bed")$size == 0) {system("echo \"No results found\" > pairtopair.bed")}
