# TOOL reseqc.R: "RNA-seq quality metrics with RseQC" (Given an RNA-seq BAM file and gene and exon locations in a BED file, this tool reports several quality metrics such as coverage uniformity, gene and junction saturation, junction annotation and alignment statistics. This tool is based on the RSeQC package.)
# INPUT alignment_file: "BAM file" TYPE GENERIC
# INPUT reference_file: "BED file" TYPE GENERIC
# OUTPUT OPTIONAL RSeQC.geneBodyCoverage.pdf
# OUTPUT OPTIONAL RSeQC.junctionSaturation_plot.pdf
# OUTPUT OPTIONAL RSeQC.splice_events.pdf
# OUTPUT OPTIONAL RSeQC.splice_junction.pdf
# OUTPUT OPTIONAL RSeQC.RKPM_saturation.pdf
# OUTPUT OPTIONAL RSeQC.bamStat.txt

# AMS 09.01.2014

# geneBody_coverage
binary <- c(file.path(chipster.tools.path, "RSeQC", "scripts", "geneBody_coverage.py"))
command <- paste(binary, "-i alignment_file -r reference_file -o RSeQC")
system(command)
source("RSeQC.geneBodyCoverage_plot.r")

# junction_saturation
binary <- c(file.path(chipster.tools.path, "RSeQC", "scripts", "junction_saturation.py"))
command <- paste(binary, "-i alignment_file -r reference_file -o RSeQC")
system(command)
source("RSeQC.junctionSaturation_plot.r")

# junction_annotation
binary <- c(file.path(chipster.tools.path, "RSeQC", "scripts", "junction_annotation.py"))
command <- paste(binary, "-i alignment_file -r reference_file -o RSeQC")
system(command)
source("RSeQC.junction_plot.r")

#RPKM_saturation
binary <- c(file.path(chipster.tools.path, "RSeQC", "scripts", "RPKM_saturation.py"))
command <- paste(binary, "-i alignment_file -r reference_file -o RSeQC")
system(command)
source("RSeQC.saturation.r")

# bam_stat
binary <- c(file.path(chipster.tools.path, "RSeQC", "scripts", "bam_stat.py"))
command <- paste(binary, "-i alignment_file 2> RSeQC.bamStat.txt")
system(command)

#system("ls > results.txt")
system("mv RSeQC.saturation.pdf RSeQC.RKPM_saturation.pdf")