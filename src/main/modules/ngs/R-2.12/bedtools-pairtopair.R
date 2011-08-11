# TOOL bedtools-pairtopair.R: "BEDTools pairToPair" (pairToPair compares two BEDPE files in search of overlaps where each end of a BEDPE feature in A
# overlaps with the ends of a feature in B. For example, using pairToPair, one could screen for the exact
# same discordant paired-end alignment in two files. This could suggest (among other things\) that the
# discordant pair suggests the same structural variation in each file/sample.)
# INPUT a "BEDPE file A. Each feature in A is compared to B in search of overlaps." TYPE GENERIC
# INPUT b "BEDPE file B." TYPE GENERIC
# OUTPUT result.txt 
# PARAMETER OPTIONAL f: "Minimum overlap" TYPE DECIMAL FROM 0 to 1 DEFAULT 0.000000001 (Minimum overlap required as a fraction of A. Default is 1E-9 (effectively 1bp\))
# PARAMETER OPTIONAL type: "Approach to reporting overlaps" TYPE [either, neither, both] DEFAULT both (either: Report overlaps if either ends of A overlap B. neither: Report A if neither end of A overlaps B. both: Report overlaps if both ends of A overlap B.)
# PARAMETER OPTIONAL is: "Force strandedness" TYPE [yes, no] DEFAULT no (Report only hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.)
# PARAMETER OPTIONAL slop:
# PARAMETER OPTIONAL ss: "Add slop based on strand" TYPE [yes, no] DEFAULT no (Add slop based to each BEDPE footprint based on strand. If strand is "+", slop is only added to the end coordinates. If strand is "-", slop is only added to the start coordinates. By default, slop is added in both directions.)
# PARAMETER OPTIONAL rdn: "Require different names" TYPE [yes, no] DEFAULT no (Require the hits to have different names (i.e. avoid self-hits\). By default, same names are allowed.)
# AMS 11.8.2011

# binary
binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "pairToPair"))

# command
command <- paste(binary, "-f", f, "-type", type, "-is", is, -"slop", slop, "-ss", ss, "-a", a, "-b", b)

# run
system(command)
