# TOOL bedtools-pairtopair.R: "BEDTools pairToPair" (pairToPair compares two BEDPE files in search of overlaps where each end of a BEDPE feature in A overlaps with the ends of a feature in B.)
# INPUT a: "BEDPE file A" TYPE GENERIC
# INPUT b: "BEDPE file B" TYPE GENERIC
# OUTPUT result.txt 
# PARAMETER OPTIONAL f: "Minimum overlap" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.000000001 (Minimum overlap required as a fraction of A. Default is 1E-9 (effectively 1bp\))
# PARAMETER OPTIONAL type: "Approach to reporting overlaps" TYPE [either, neither, both] DEFAULT both (either: Report overlaps if either ends of A overlap B. neither: Report A if neither end of A overlaps B. both: Report overlaps if both ends of A overlap B.)
# PARAMETER OPTIONAL is: "Force strandedness" TYPE [yes, no] DEFAULT no (Report only hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.)
# PARAMETER OPTIONAL rdn: "Require different names" TYPE [yes, no] DEFAULT no (Require the hits to have different names (i.e. avoid self-hits\). By default, same names are allowed.)

# binary
binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "pairToPair"))

# command
command <- paste(binary, "-f", f, "-type", type, "-is", is, "-rdn", rdn, "-a", a, "-b", b)

# run
system(command)
