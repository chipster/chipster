# TOOL reseqc.R: "RNA-seq quality metrics with RseQC" (Given an RNA-seq BAM file and gene and exon locations in a BED file, this tool reports several quality metrics such as coverage uniformity, gene and junction saturation, junction annotation and alignment statistics. This tool is based on the RSeQC package.)
# INPUT alignment_file: "BAM file" TYPE GENERIC
# INPUT reference_file: "BED file" TYPE GENERIC
# OUTPUT OPTIONAL RSeQC.geneBodyCoverage.pdf
# OUTPUT OPTIONAL RSeQC.junctionSaturation_plot.pdf
# OUTPUT OPTIONAL RSeQC.splice_events.pdf
# OUTPUT OPTIONAL RSeQC.splice_junction.pdf
# OUTPUT OPTIONAL RSeQC.RPKM_saturation.pdf
# OUTPUT OPTIONAL RSeQC.bamStat.txt
# OUTPUT OPTIONAL RSeQC.inner_distance_plot.pdf
# PARAMETER OPTIONAL paired: "Generate inner distance plot" TYPE [yes, no] DEFAULT no (Calculate the inner distance (or insert size\) between two paired RNA reads. The distance is the mRNA length between two paired fragments.)

# AMS 09.01.2014
# AMS 23.05.2014 added inner distance plot
# AMS 03.12.2014 improved error handling for the plots

# geneBody_coverage
binary <- c(file.path(chipster.tools.path, "RSeQC", "scripts", "geneBody_coverage.py"))
command <- paste(binary, "-i alignment_file -r reference_file -o RSeQC")
system(command)
try(source("RSeQC.geneBodyCoverage_plot.r"), silent=TRUE)

# junction_saturation
binary <- c(file.path(chipster.tools.path, "RSeQC", "scripts", "junction_saturation.py"))
command <- paste(binary, "-i alignment_file -r reference_file -o RSeQC")
system(command)
try(source("RSeQC.junctionSaturation_plot.r"), silent=TRUE)

# junction_annotation
binary <- c(file.path(chipster.tools.path, "RSeQC", "scripts", "junction_annotation.py"))
command <- paste(binary, "-i alignment_file -r reference_file -o RSeQC")
system(command)
try(source("RSeQC.junction_plot.r"), silent=TRUE)

#RPKM_saturation
binary <- c(file.path(chipster.tools.path, "RSeQC", "scripts", "RPKM_saturation.py"))
command <- paste(binary, "-i alignment_file -r reference_file -o RSeQC")
system(command)
try(source("RSeQC.saturation.r"), silent=TRUE)
system("mv RSeQC.saturation.pdf RSeQC.RPKM_saturation.pdf")

# bam_stat
binary <- c(file.path(chipster.tools.path, "RSeQC", "scripts", "bam_stat.py"))
command <- paste(binary, "-i alignment_file 2> RSeQC.bamStat.txt")
system(command)

# inner_distance
if (paired == "yes"){
	binary <- c(file.path(chipster.tools.path, "RSeQC", "scripts", "inner_distance.py"))
	command <- paste(binary, "-i alignment_file -r reference_file -o RSeQC")
	system(command)
	try(source("RSeQC.inner_distance_plot.r"), silent=TRUE)	
}
