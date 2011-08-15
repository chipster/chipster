# TOOL bedtools-closestbed.R: "BEDTools closestBed" (For each feature in A, finds the closest feature (upstream or downstream\) in B.)
# INPUT file.a: "BED file A" TYPE GENERIC
# INPUT file.b: "BED file B" TYPE GENERIC
# OUTPUT closestbed.bed 
# PARAMETER OPTIONAL s: "Force strandedness" TYPE [yes, no] DEFAULT no (Force strandedness. That is, find the closest feature in B that overlaps A on the same strand. By default, overlaps are reported without respect to strand.)
# PARAMETER OPTIONAL d: "Report distance" TYPE [yes, no] DEFAULT no (In addition to the closest feature in B, report its distance to A as an extra column. The reported distance for overlapping features will be 0.)
# PARAMETER OPTIONAL t: "Approach to reporting multiple overlaps" TYPE [all, first, last] DEFAULT all (How ties for closest feature are handled. This occurs when two features in B have exactly the same overlap with A. By default, all such features in B are reported. The options are: all (Report all ties\), first (Report the first tie that occurred in the B file\) and  last (Report the last tie that occurred in the B file\).)


# binary
binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "closestBed"))

# optional options
options <- paste("")
if (s == "yes") {options <- paste(options, "-s")}
if (d == "yes") {options <- paste(options, "-d")}

# common options
options <- paste(options, "-t", t, "-a file.a -b file.b")

# command
command <- paste(binary, options, " >closestbed.bed")

# run
system(command)
