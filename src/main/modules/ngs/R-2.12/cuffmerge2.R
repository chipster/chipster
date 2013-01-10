# TOOL cuffmerge2.R: "Merge transcript assemblies with Cuffmerge" (This tool allows you to merge transcript GTF files obtained by Cufflinks from several samples into one, so that you can use it in differential expression analysis using Cuffdiff.)
# INPUT annotation{...}.gtf: "GTF files" TYPE GTF
# INPUT OPTIONAL ownref.fa: "Reference sequence FASTA" TYPE FASTA  
# OUTPUT OPTIONAL merged.gtf  
# OUTPUT OPTIONAL skipped.gtf  
# OUTPUT OPTIONAL transcripts.gtf
# PARAMETER genome: "Genome" TYPE [hg19: "Human genome (hg19\)", mm9: "Mouse genome (mm9\)", mm10: "Mouse genome (mm10\)", rn4: "Rat genome (rn4\)"] DEFAULT hg19 (Genome used for bias correction.)

# AMS 21.11.2012
# LUE TÄMÄ poistin ne tsv:t outputeista, varmaan pitää poistaa skriptistäkin? mitä nää skipped ja transcripts gtf:t on? Genomi ei ole tässä bias correctionia varten, joten sen voisi kommentoida ulos.

# binary
cuffmerge.binary <- c(file.path(chipster.tools.path, "cufflinks-2.0.2.Linux_x86_64", "cuffmerge"))

# Set PATH so cuffmerge can find gtf_to_sam
cufflinkspath <- c(file.path(chipster.tools.path, "cufflinks-2.0.2.Linux_x86_64"))
setpathc1 <- c("export PATH=$PATH:", cufflinkspath, ";")
setpathcommand <- paste(setpathc1 , collapse = '')


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
command <- paste(setpathcommand, cuffmerge.binary, "-s", genomefile, "assemblies.txt")

# run
#stop(paste('CHIPSTER-NOTE: ', command))
system(command)

# Rename files
if (file.exists("merged_asm/genes.fpkm_tracking") && file.info("merged_asm/genes.fpkm_tracking")$size > 0) {
	system("mv merged_asm/genes.fpkm_tracking genes.fpkm_tracking.tsv")
}
if (file.exists("merged_asm/isoforms.fpkm_tracking") && file.info("merged_asm/isoforms.fpkm_tracking")$size > 0) {
	system("mv merged_asm/isoforms.fpkm_tracking isoforms.fpkm_tracking.tsv")
}
if (file.exists("merged_asm/merged.gtf") && file.info("merged_asm/merged.gtf")$size > 0) {
	system("mv merged_asm/merged.gtf merged.gtf")
}
if (file.exists("merged_asm/skipped.gtf") && file.info("merged_asm/skipped.gtf")$size > 0) {
	system("mv merged_asm/skipped.gtf skipped.gtf")
}
if (file.exists("merged_asm/transcripts.gtf") && file.info("merged_asm/transcripts.gtf")$size > 0) {
	system("mv merged_asm/transcripts.gtf transcripts.gtf")
}
 