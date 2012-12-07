# TOOL cuffmerge2.R: "Cuffmerge2" (Cuffmerge2)
# INPUT sequence.fa: "FASTA sequence" TYPE GENERIC
# INPUT transcripts1.gtf: "GTF file" TYPE GENERIC
# INPUT transcripts2.gtf: "GTF file" TYPE GENERIC
# OUTPUT OPTIONAL genes.fpkm_tracking.tsv  
# OUTPUT OPTIONAL isoforms.fpkm_tracking.tsv  
# OUTPUT OPTIONAL merged.gtf  
# OUTPUT OPTIONAL skipped.gtf  
# OUTPUT OPTIONAL transcripts.gtf
# OUTPUT OPTIONAL assemblies.txt

# AMS 21.11.2012

# binary
cuffmerge.binary <- c(file.path(chipster.tools.path, "cufflinks-2.0.2.Linux_x86_64", "cuffmerge"))

# Make gtf list file
system("ls *.gtf > assemblies.txt")

# command
command <- paste(cuffmerge.binary, "-s", "sequence.fa", "assemblies.txt")

# run
#stop(paste('CHIPSTER-NOTE: ', command))
system(command)

system("ls -l > skipped.gtf")
system("ls -l merged_asm > merged.gtf")

# run
system(command)

# Rename files
if (file.info("merged_asm/genes.fpkm_tracking")$size > 0) {
	system("mv merged_asm/genes.fpkm_tracking genes.fpkm_tracking.tsv")
}
if (file.info("merged_asm/isoforms.fpkm_tracking")$size > 0) {
	system("mv merged_asm/isoforms.fpkm_tracking isoforms.fpkm_tracking.tsv")
}
if (file.info("merged_asm/merged.gtf")$size > 0) {
	system("mv merged_asm/merged.gtf merged.gtf")
}
if (file.info("merged_asm/skipped.gtf")$size > 0) {
	system("mv merged_asm/skipped.gtf skipped.gtf")
}
if (file.info("merged_asm/transcripts.gtf")$size > 0) {
	system("mv merged_asm/transcripts.gtf transcripts.gtf")
}
