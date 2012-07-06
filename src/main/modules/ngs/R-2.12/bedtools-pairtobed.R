# TOOL bedtools-pairtobed.R: "Overlap BEDPE with BED" (Reports overlaps between a BEDPE file and a BED/GFF/VCF file. This tool is based on the BEDTools package.)
# INPUT file.a: "BEDPE file" TYPE GENERIC
# INPUT file.b: "BED file" TYPE GENERIC
# OUTPUT pairtobed.txt
# PARAMETER abam: "File A is BAM format" TYPE [yes, no] DEFAULT no (Select yes if file A is BAM format.)
# PARAMETER OPTIONAL ubam: "Write uncompressed BAM output" TYPE [yes,no] DEFAULT no (Write uncompressed BAM output. Default is to write compressed BAM.) 
# PARAMETER OPTIONAL bedpe: "When using BAM input, write output as BEDPE." TYPE [yes,no] DEFAULT no (When using BAM input, write output as BEDPE. The default is to write output in BAM.) 
# PARAMETER OPTIONAL ed: "Use BAM total edit distance for BEDPE score" TYPE [yes, no] DEFAULT no (Use BAM total edit distance (NM tag\) for BEDPE score. Default for BEDPE is to use the minimum of of the two mapping qualities for the pair. When this option is used the total edit distance from the two mates is reported as the score.)
# PARAMETER OPTIONAL f: "Minimum overlap" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.000000001 (Minimum overlap required as a fraction of A (e.g. 0.05\). Default is 1E-9 (effectively 1bp\))
# PARAMETER OPTIONAL s: "Force strandedness" TYPE [yes, no] DEFAULT no (Enforce strandedness when finding overlaps. Default is to ignore stand. Not applicable with  output type inspan or outspan.)
# PARAMETER type: "Approach to reporting overlaps" TYPE [either, neither, both, xor, notboth, ispan, ospan, notispan, notospan] DEFAULT either (Approach to reporting overlaps between BEDPE and BED.
# either: Report overlaps if either end of A overlaps B.
# neither: Report A if neither end of A overlaps B.
# both: Report overlaps if both ends of A overlap  B.
# xor: Report overlaps if one and only one end of A overlaps B.
# notboth: Report overlaps if neither end or one and only one end of A overlap B.  That is, xor + neither.
# ispan: Report overlaps between [end1, start2] of A and B. Note: If chrom1 <> chrom2, entry is ignored.
# ospan: Report overlaps between [start1, end2] of A and B. Note: If chrom1 <> chrom2, entry is ignored.
# notispan: Report A if ispan of A doesn't overlap B. Note: If chrom1 <> chrom2, entry is ignored.
# notospan: Report A if ospan of A doesn't overlap B. Note: If chrom1 <> chrom2, entry is ignored.)

# AMS 23.4.2012

# binary
binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "pairToBed"))

# optional options
options <- paste("")
options <- paste("")
if (abam == "yes") {
	if (ubam == "yes") {options <- paste(options, "-ubam")}
	if (bedpe == "yes") {options <- paste(options, "-bedpe")}
}
if (ed == "yes") {options <- paste(options, "-ed")}
options <- paste(options, "-f", f) 
if (s == "yes") {options <- paste(options, "-s")}
options <- paste(options, "-type", type) 

# input files
# input files
if (abam == "yes") {options <- paste(options, "-abam file.a -b file.b")}
if (abam == "no") {options <- paste(options, "-a file.a -b file.b")}


# command
command <- paste(binary, options, " > pairtobed.bed")

# run
system(command)
if (file.info("pairtobed.bed")$size == 0) {system("echo \"No results found\" > pairtobed.bed")}