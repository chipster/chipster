# TOOL cufflinks2.R: "Assemble reads into transcripts using Cufflinks" (Given aligned RNA-seq reads as a BAM file, Cufflinks assembles the alignments into a parsimonious set of transcripts. It then estimates the relative abundances of these transcripts based on how many reads support each one. It is recommended to create the input BAM files using the TopHat aligner. You can view the resulting GTF file in Chipster genome browser to explore the structure of the genes. If you have GTF files from several samples, you can merge them using the Cuffmerge tool, and use the merged GTF file in differential expression analysis using Cuffdiff.)
# INPUT alignment.bam TYPE BAM
# OUTPUT OPTIONAL genes.fpkm_tracking.tsv  
# OUTPUT OPTIONAL isoforms.fpkm_tracking.tsv
# OUTPUT OPTIONAL skipped.gtf
# OUTPUT OPTIONAL transcripts.gtf
# PARAMETER OPTIONAL normalize: "Upper-quartile normalization " TYPE [yes, no] DEFAULT yes (Upper quartile normalization can improve robustness of differential expression calls for less abundant genes and transcripts. It excludes very abundant genes when normalizing expression values for the number of reads in each sample by using the upper quartile of the number of fragments mapping to individual loci.)
# PARAMETER OPTIONAL mmread: "Enable multi-mapped read correction" TYPE [yes, no] DEFAULT no (By default, Cufflinks will uniformly divide each multi-mapped read to all of the positions it maps to. If multi-mapped read correction is enabled, Cufflinks will re-estimate the transcript abundances dividing each multi-mapped read probabilistically based on the initial abundance estimation, the inferred fragment length and fragment bias, if bias correction is enabled.)

# AMS 21.11.2012

# binary
cufflinks.binary <- c(file.path(chipster.tools.path, "cufflinks-2.0.2.Linux_x86_64", "cufflinks"))

# options
cufflinks.options <- ""
if (normalize == "yes") {
	cufflinks.options <- paste(cufflinks.options, "-N")
}
if (mmread == "yes") {
	cufflinks.options <- paste(cufflinks.options, "-u")
}

# command
command <- paste(cufflinks.binary, cufflinks.options, "-q", "-o tmp", "alignment.bam")

# run
system(command)

# Rename files
if (file.exists("tmp/genes.fpkm_tracking") && file.info("tmp/genes.fpkm_tracking")$size > 0) {
	system("mv tmp/genes.fpkm_tracking genes.fpkm_tracking.tsv")
}
if (file.exists("tmp/isoforms.fpkm_tracking") && file.info("tmp/isoforms.fpkm_tracking")$size > 0) {
	system("mv tmp/isoforms.fpkm_tracking isoforms.fpkm_tracking.tsv")
}
if (file.exists("tmp/skipped.gtf") && file.info("tmp/skipped.gtf")$size > 0) {
	system("mv tmp/skipped.gtf skipped.gtf")
}
if (file.exists("tmp/transcripts.gtf") && file.info("tmp/transcripts.gtf")$size > 0) {
	system("mv tmp/transcripts.gtf transcripts.gtf")
}

