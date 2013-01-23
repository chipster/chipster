# TOOL cuffmerge2.R: "Merge transcript assemblies with Cuffmerge" (Given several transcript GTF files obtained by Cufflinks, Cuffmerge merges them into one. The merged GTF file can be used in differential expression analysis with Cuffdiff.)
# INPUT annotation{...}.gtf: "GTF files" TYPE GTF
# OUTPUT OPTIONAL merged.gtf  

# AMS 21.11.2012
# AMS 11.01.2013 Removed unnecessary outputs
# EK 21.1.2013

# binary
cuffmerge.binary <- c(file.path(chipster.tools.path, "cufflinks2", "cuffmerge"))

# Set PATH so cuffmerge can find gtf_to_sam
cufflinkspath <- c(file.path(chipster.tools.path, "cufflinks2"))
setpathc1 <- c("export PATH=$PATH:", cufflinkspath, ";")
setpathcommand <- paste(setpathc1 , collapse = '')

# Make gtf list file
system("ls *.gtf > assemblies.txt")

# command
command <- paste(setpathcommand, cuffmerge.binary, "assemblies.txt")

# run
#stop(paste('CHIPSTER-NOTE: ', command))
system(command)

# Rename files
if (file.exists("merged_asm/merged.gtf") && file.info("merged_asm/merged.gtf")$size > 0) {
	system("mv merged_asm/merged.gtf merged.gtf")
}
