# TOOL bedtools-intersectbed.R: "Intersect BED" (Report overlaps between two feature files such as BED, VCF and GTF. One of the files can also be BAM. Note that when A and B files are compared, the B file is loaded into memory. Therefore to minimize memory usage, one should set the smaller of the two files as the B file. This tool is based on the BEDTools package.)
# INPUT file.a: "Input file A" TYPE GENERIC
# INPUT file.b: "Input file B, the smaller file" TYPE GENERIC
# OUTPUT OPTIONAL intersectbed.bed 
# OUTPUT OPTIONAL intersectbed.bam
# OUTPUT OPTIONAL intersectbed.bam.bai
# OUTPUT OPTIONAL error.txt
# PARAMETER wa: "Write the original entry in A for each overlap" TYPE [yes, no] DEFAULT no (Report the original A feature when an overlap is found. The entire A feature is reported, not just the portion that overlaps with the B feature.)
# PARAMETER u: "Write the original A entry once if any overlaps found in B" TYPE [yes, no] DEFAULT no (Write the original A entry once if any overlaps found in B)
# PARAMETER v: "Only report those entries in A that have no overlaps with B" TYPE [yes, no] DEFAULT no (Only report those entries in A that have no overlaps with B)
# PARAMETER wb: "Write the original entry in B for each overlap" TYPE [yes, no] DEFAULT no (Report the original B feature when an overlap is found. The entire B feature is reported, not just the portion that overlaps with the A feature.)
# PARAMETER s: "Force strandedness" TYPE [yes, no] DEFAULT no (Only report hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.)
# PARAMETER abam: "File A is BAM format" TYPE [yes, no] DEFAULT no (Select yes if file A is a BAM file.)
# PARAMETER OPTIONAL wo: "Write the original A and B entries for overlapped A features" TYPE [yes, no] DEFAULT no (Write the original A and B entries plus the number of base pairs of overlap between the two feature. Only A features with overlap are reported.)
# PARAMETER OPTIONAL wao: "Write the original A and B entries for all A features" TYPE [yes, no] DEFAULT no (Write the original A and B entries plus the number of base pairs of overlap between the two features. However, A features w/o overlap are also reported with a NULL B feature and overlap = 0.)
# PARAMETER OPTIONAL c: "For each entry in A, report the number of overlaps with B" TYPE [yes, no] DEFAULT no (For each entry in A, report the number of overlaps with B)
# PARAMETER OPTIONAL f: "Minimum overlap required as a fraction of A" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.000000001 (Minimum overlap required as a fraction of A. Default is 1E-9 (effectively 1bp\))
# PARAMETER OPTIONAL r: "Require that the fraction overlap be reciprocal for A and B" TYPE [yes, no] DEFAULT no (Require that the fraction overlap be reciprocal for A and B. In other words, if minumum overlap is 0.90 and this option is selected, this requires that B overlap 90% of A and A also overlaps 90% of B.)
# PARAMETER OPTIONAL ubam: "Write uncompressed BAM output" TYPE [yes, no] DEFAULT no (Write uncompressed BAM output. Default is to write compressed BAM.)
# PARAMETER OPTIONAL bed: "When using BAM input, write output as BED" TYPE [yes, no] DEFAULT no (When using BAM input, the default is to write output in BAM.)
# PARAMETER OPTIONAL split: "Treat split BAM or BED12 entries as distinct BED intervals" TYPE [yes, no] DEFAULT no (Treat "split" BAM (i.e., having an “N” CIGAR operation\) or BED12 entries as distinct BED intervals.)

# AMS 23.4.2012
# AMS 11.10.2012 Fixed BAM file support
# AMS 23.9.2013 Improved output/error file handling

# binary
binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "intersectBed"))
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# options
options <- paste("")
outputfiletype <- "bed"
if (abam == "yes") {
	outputfiletype <- "bam"
	if (ubam == "yes") {
		options <- paste(options, "-ubam")		
	}
	if (bed == "yes") {
		outputfiletype <- "bed"
		options <- paste(options, "-bed")		
	}
}
if (wa == "yes") {options <- paste(options,"-wa")}
if (wb == "yes") {options <- paste(options,"-wb")}
if (wo == "yes") {options <- paste(options,"-wo")}
if (wao == "yes") {options <- paste(options,"-wao")}
if (u == "yes") {options <- paste(options, "-u")}
if (c == "yes") {options <- paste(options, "-c")}
if (v == "yes") {options <- paste(options, "-v")}
options <- paste(options, "-f", f)
if (r == "yes") {options <- paste(options, "-r")}
if (s == "yes") {options <- paste(options, "-s")}
if (split == "yes") {options <- paste(options, "-split")}

# input files
if (abam == "yes") {options <- paste(options, "-abam file.a -b file.b")}
if (abam == "no") {options <- paste(options, "-a file.a -b file.b")}

# command
command <- paste(binary, options, "> intersectbed.tmp 2> error.tmp")

#stop(paste('CHIPSTER-NOTE: ', command))

# run
system(command)

# Generate output/error message
if (file.info("intersectbed.tmp")$size > 0) {
	if (outputfiletype == "bed"){
		system("mv intersectbed.tmp intersectbed.bed")
	}
	if (outputfiletype == "bam"){
	system("mv intersectbed.tmp intersectbed.bam")
	}	
} else if (file.info("error.tmp")$size > 0) {
	system("mv error.tmp error.txt")
} else{
	system("echo \"# No results found\" > error.txt")
}

# Index bam
if (file.exists("intersectbed.bam")){
	system(paste(samtools.binary, "index intersectbed.bam > intersectbed.bam.bai"))
}

# Handle output names
source(file.path(chipster.common.path, "tool-utils.R"))

# read input names
inputnames <- read_input_definitions()

basename  <- strip_name(inputnames$file.a)

# Make a matrix of output names
outputnames <- matrix(NA, nrow=2, ncol=2)
outputnames[1,] <- c("intersectbed.bam", paste(basename, "_intersect.bam", sep =""))
outputnames[2,] <- c("intersectbed.bam.bai", paste(basename, "_intersect.bam.bai", sep =""))

# Write output definitions file
write_output_definitions(outputnames)