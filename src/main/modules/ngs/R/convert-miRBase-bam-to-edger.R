# TOOL convert-miRBase-bam-to-edger.R: "Convert miRBased BAM file to count table" (This tool takes a BAM file as an input, calculates the number of times each miRNA is identified, and removes the ones for which the count is under the user defined threshold.)
# INPUT bam_file.bam: "Alignment against miRBase in BAM format" TYPE GENERIC
# OUTPUT miRNA-counts.tsv: "A count file suitable for differential expression analysis"
# PARAMETER count_limit: "Count limit" TYPE INTEGER FROM 0 TO 1000 DEFAULT 10 (Keep miRNAs which have more reads than this.)

# EK 07.07.2011
# MK 13.05.2013, fix bug in header formats
# MK 12.05.2014, added check for file size 

# Convert to SAM, grep miRNA name (removes unaligned reads), count, filter on tag number
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))
counts.command <- paste(samtools.binary, "view bam_file.bam | awk '{print $3}' | grep -v '^*' | sort -k1 | uniq -c | awk '{if($1>", count_limit, ")print $2\"\t\"$1}'> counts.tsv")
system(counts.command)

input.file <- "counts.tsv"
if(file.info(input.file)$size > 0) {
	dat <- read.table(input.file, header=F, sep="\t", row.names=NULL)
	colnames(dat) <- c("id", "count")
	write.table(data.frame(dat), file="miRNA-counts.tsv", col.names=T, quote=F, sep="\t", row.names=F)
} else {
	stop("CHIPSTER-NOTE: No miRNA in your BAM file had the minimum number of reads required.")	
}

# Add column headers
#headers <- paste("id\t","count", sep="")
#input.file <- "counts.tsv"
#header.file <- "header_file"
#system(paste("echo \"", headers, "\"", ">", header.file))
#merge.command <- paste("cat", header.file, input.file, "> miRNA-counts.tsv")
#system(merge.command)

#headers <- paste("id\t","count", sep="")
#system(paste("echo \"", headers, "\"", "> header.file"))
#merge.command <- paste("cat", header.file, "counts.tsv > edgeR-input.tsv")
#system(merge.command)

# EOF

