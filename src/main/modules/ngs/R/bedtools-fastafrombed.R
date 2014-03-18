# TOOL bedtools-fastafrombed.R: "FASTA from BED" (Extracts sequences from a FASTA file for each of the intervals defined in a BED file. The headers in the input FASTA file must exactly match the chromosome column in the BED file. This tool is based on the BEDTools package.)
# INPUT f.file: "FASTA file" TYPE GENERIC
# INPUT b.file: "BED/GFF/VCF file" TYPE GENERIC
# OUTPUT fastafrombed.txt
# PARAMETER OPTIONAL name: "Use the name field for the FASTA header" TYPE [yes, no] DEFAULT no (Use the name field for the FASTA header)
# PARAMETER OPTIONAL tab: "Write output in TAB delimited format" TYPE [yes, no] DEFAULT no (Write output in TAB delimited format. Default is FASTA format.)
# PARAMETER OPTIONAL s: "Force strandedness" TYPE [yes, no] DEFAULT no (Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented. By default, strand information is ignored.) 

# AMS 23.4.2012

# binary
binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "fastaFromBed"))

# options
options <- paste("")
if (name == "yes") {options <- paste(options, "-name")}
if (tab == "yes") {options <- paste(options, "-tab")}
if (s == "yes") {options <- paste(options, "-s")}

# input files
options <- paste(options, "-fi f.file", "-bed b.file", "-fo fastafrombed.txt")

# command
command <- paste(binary, options)

# run
system("dos2unix f.file")
system(command)
