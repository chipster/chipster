# TOOL cufflinks2.R: "Assemble reads into transcripts using Cufflinks" (Given aligned RNA-Seq reads, Cufflinks assembles the alignments into a parsimonious set of transcripts. It then estimates the relative abundances of these transcripts based on how many reads support each one, taking into account biases in library preparation protocols. It is recommended to use the TopHat aligner to map your reads to the reference genome. You can view the resulting GTF file in Chipster genome browser to explore the structure of the genes. If you have GTF files from several samples, you can merge them to one using the Cuffmerge tool, and use it in differential expression analysis using Cuffdiff.)
# INPUT alignment.bam TYPE GENERIC
# OUTPUT OPTIONAL genes.fpkm_tracking.tsv  
# OUTPUT OPTIONAL isoforms.fpkm_tracking.tsv
# OUTPUT OPTIONAL skipped.gtf
# OUTPUT OPTIONAL transcripts.gtf

# AMS 21.11.2012

# binary
cufflinks.binary <- c(file.path(chipster.tools.path, "cufflinks-2.0.2.Linux_x86_64", "cufflinks"))

# command
command <- paste(cufflinks.binary, "-q", "-o tmp", "alignment.bam")

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

