# TOOL cufflinks2.R: "Cufflinks2" (Cufflinks2)
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
if (file.info("tmp/genes.fpkm_tracking")$size > 0) {
	system("mv tmp/genes.fpkm_tracking genes.fpkm_tracking.tsv")
}
if (file.info("tmp/isoforms.fpkm_tracking")$size > 0) {
	system("mv tmp/isoforms.fpkm_tracking isoforms.fpkm_tracking.tsv")
}
if (file.info("tmp/skipped.gtf")$size > 0) {
	system("mv tmp/skipped.gtf skipped.gtf")
}
if (file.info("tmp/transcripts.gtf")$size > 0) {
	system("mv tmp/transcripts.gtf transcripts.gtf")
}