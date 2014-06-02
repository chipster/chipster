# TOOL convert-RNA-bam-to-edger.R: "Map aligned reads to genes with coverageBed using own BED" (Calculates how many reads in a BAM file map to each gene. You have to provide the gene locations in the BED format. Please note that the chromosome names have to be same in the BED and BAM files. This tool is based on the coverageBed tool of the BEDTools package. In order to use the output in edgeR, you need to select all samples and run the tool \"Utilities - Define NGS experiment\".)
# INPUT file.a: "BAM file" TYPE GENERIC
# INPUT file.b: "BED file" TYPE GENERIC
# OUTPUT edgeR-input.tsv
# PARAMETER OPTIONAL s: "Force strandedness" TYPE [yes, no] DEFAULT no (Force strandedness. That is, only count reads in BAM that overlap genes on the same strand in the BED file. By default, hits are counted regardless of strand.)
# PARAMETER OPTIONAL split: "Treat split BAM or BED12 entries as distinct BED intervals" TYPE [yes, no] DEFAULT no (Treat "split" BAM or BED12 entries as distinct BED intervals when computing coverage. For BAM files, this uses the CIGAR N and D operations to infer the blocks for computing coverage. For BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds fields (i.e., columns 10,11,12\).)

# MG 22.8.2011
# modified from the coveragebed tool by AMS
# EK 10.1.2012

# setup parameters
abam <- "yes"

# binary
binary <- c(file.path(chipster.tools.path, "bedtools", "bin", "coverageBed"))

# optional options
options <- paste("")
if (s == "yes") {options <- paste(options, "-s")}
# if (hist == "yes") {options <- paste(options, "-hist")}
# if (d == "yes") {options <- paste(options, "-d")}
if (split == "yes") {options <- paste(options, "-split")}

# input files
if (abam == "yes") {options <- paste(options, "-abam file.a -b file.b")}
if (abam == "no") {options <- paste(options, "-a file.a -b file.b")}

# command
command <- paste(binary, options, " > coveragebed.tsv")

# run
system(command)

# bring in file to R environment for formating
file <- c("coveragebed.tsv")
dat <- read.table(file, header=F, sep="\t")
id_list <- paste(dat$V1,dat$V2,dat$V3,dat$V4, sep="_")
length_list <- dat$V3-dat$V2+1
results_table <- data.frame	(
		id=id_list,
		sequence=dat$V4,
		chr=dat$V1,
		start=dat$V2,
		end=dat$V3,
		length=length_list,
		count=dat$V5
)

# write result table to output
write.table(results_table, file="edgeR-input.tsv", col.names=T, quote=F, sep="\t", row.names=F)

# EOF
