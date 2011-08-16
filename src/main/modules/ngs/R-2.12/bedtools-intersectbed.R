# TOOL bedtools-intersectbed.R: "BEDTools intersectBed" (Report overlaps between two feature files.)
# INPUT file.a: "Input file A" TYPE GENERIC
# INPUT file.b: "Input file B" TYPE GENERIC
# OUTPUT intersectbed.bed 
# PARAMETER abam: "File A is BAM format" TYPE [yes, no] DEFAULT no (Select yes if file A is BAM format.)
# PARAMETER OPTIONAL ubam: "Write uncompressed BAM output" TYPE [yes, no] DEFAULT no (Write uncompressed BAM output. Default is to write compressed BAM.)
# PARAMETER OPTIONAL bed: "Write output as BED" TYPE [yes, no] DEFAULT no (When using BAM input, write output as BED. The default is to write output in BAM when using BAM input.)
# PARAMETER OPTIONAL wa: "Write the original entry in A for each overlap" TYPE [yes, no] DEFAULT no (Report the original A feature when an overlap is found. The entire A feature is reported, not just the portion that overlaps with the B feature.)
# PARAMETER OPTIONAL wb: "Write the original entry in B for each overlap" TYPE [yes, no] DEFAULT no (Report the original B feature when an overlap is found. The entire B feature is reported, not just the portion that overlaps with the A feature.)
# PARAMETER OPTIONAL wo: "Write the original A and B entries" TYPE [yes, no] DEFAULT no (Write the original A and B entries plus the number of base pairs of overlap between the two feature. Only A features with overlap are reported.)
# PARAMETER OPTIONAL wao: "Write the original A and B entries" TYPE [yes, no] DEFAULT no (Write the original A and B entries plus the number of base pairs of overlap between the two features. However, A features w/o overlap are also reported with a NULL B feature and overlap = 0.)
# PARAMETER OPTIONAL u: "Write the original A entry once if any overlaps found in B" TYPE [yes, no] DEFAULT no (Write the original A entry once if any overlaps found in B)
# PARAMETER OPTIONAL c: "For each entry in A, report the number of overlaps with B" TYPE [yes, no] DEFAULT no (For each entry in A, report the number of overlaps with B)
# PARAMETER OPTIONAL v: "Only report those entries in A that have no overlaps with B" TYPE [yes, no] DEFAULT no (Only report those entries in A that have no overlaps with B)
# PARAMETER OPTIONAL f: "Minimum overlap required as a fraction of A" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.000000001 (Minimum overlap required as a fraction of A. Default is 1E-9 (effectively 1bp\))
# PARAMETER OPTIONAL r: "Require that the fraction overlap be reciprocal for A and B" TYPE [yes, no] DEFAULT no (Require that the fraction overlap be reciprocal for A and B. In other words, if minumum overlap is 0.90 and this option is selected, this requires that B overlap 90% of A and A also overlaps 90% of B.)
# PARAMETER OPTIONAL s: "Force strandedness" TYPE [yes, no] DEFAULT no (Only report hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.)
# PARAMETER OPTIONAL split: "Treat split BAM or BED12 entries as distinct BED intervals" TYPE [yes, no] DEFAULT no (Treat "split" BAM (i.e., having an “N” CIGAR operation\) or BED12 entries as distinct BED intervals.)

# binary
binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "intersectBed"))

# options
options <- paste("")
if (abam == "yes") {
	if (ubam == "yes") {options <- paste(options, "-ubam")}
	if (bed == "yes") {options <- paste(options, "-bed")}
}
if (wa == "yes") {options <- paste(options,"-wa")}
if (wb == "yes") {options <- paste(options,"-wb")}
if (wo == "yes") {options <- paste(options,"-wo")}
if (wao == "yes") {options <- paste(options,"-wao")}
if (u == "yes") {options <- paste(options, "-u")}
if (c == "yes") {options <- paste(options, "-c")}
if (v == "yes") {options <- paste(options, "-v")}
if (f == "yes") {options <- paste(options, "-f", f)}
if (r == "yes") {options <- paste(options, "-r")}
if (s == "yes") {options <- paste(options, "-s")}
if (split == "yes") {options <- paste(options, "-split")}

# input files
if (abam == "yes") {options <- paste(options, "-abam file.a -b file.b")}
if (abam == "no") {options <- paste(options, "-a file.a -b file.b")}

# command
command <- paste(binary, options, "> intersectbed.bed")

# run
system(command)