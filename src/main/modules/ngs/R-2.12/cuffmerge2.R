# TOOL cuffmerge2.R: "Cuffmerge2" (Cuffmerge2)
# INPUT annotation{...}.gtf: "GTF files" TYPE GENERIC
# INPUT OPTIONAL ownref.fa: "Reference sequence FASTA" TYPE FASTA
# OUTPUT OPTIONAL ls.txt
# OUTPUT OPTIONAL genes.fpkm_tracking.tsv  
# OUTPUT OPTIONAL isoforms.fpkm_tracking.tsv  
# OUTPUT OPTIONAL merged.gtf  
# OUTPUT OPTIONAL skipped.gtf  
# OUTPUT OPTIONAL transcripts.gtf
# PARAMETER genome: "Genome" TYPE [hg19: "Human genome (hg19\)", mm9: "Mouse genome (mm9\)", mm10: "Mouse genome (mm10\)", rn4: "Rat genome (rn4\)"] DEFAULT hg19 (Genome used for bias correction.)

# AMS 21.11.2012

# binary
cuffmerge.binary <- c(file.path(chipster.tools.path, "cufflinks-2.0.2.Linux_x86_64", "cuffmerge"))

# Make gtf list file
system("ls *.gtf > assemblies.txt")

# Reference sequence
cuffmerge.options <- ""
if (file.exists("ownref.fa")){
	genomefile <- "ownref.fa"
}else{
	if (genome == "hg19"){
		genomefile <- "hg19.fa"
	}
	if (genome == "mm9"){
		genomefile <- "mm9.fa"
	}
	if (genome == "mm10"){
		genomefile <- "mm10.fa"
	}
	if (genome == "rn4"){
		genomefile <- "rn4.fa"
	}
	genomefile <- c(file.path(chipster.tools.path, "genomes", "fasta", "nochr", genomefile))	
}


# command
command <- paste(cuffmerge.binary, "-s", genomefile, "assemblies.txt")

# run
#stop(paste('CHIPSTER-NOTE: ', command))
#system(command)
 system("ls -l > ls.txt")
 
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
