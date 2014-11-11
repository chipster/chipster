# TOOL bwa-mem-paired-end-with-index-building.R: "BWA MEM for paired-end reads and own genome" (BWA aligns reads to genomes and transcriptomes using BWA MEM algorithm. Results are sorted and indexed bam files, which are ready for viewing in the Chipster genome browser. 
# Note that this BWA MEM tool requires that you have imported the reference genome to Chipster in fasta format. If you would like to align paired-end reads against publicly available genomes, please use the tool \"BWA MEM for paired-end reads\".)
# INPUT reads1.txt: "Reads to align" TYPE GENERIC 
# INPUT reads2.txt: "Reads to align" TYPE GENERIC 
# INPUT genome.txt: "Reference genome" TYPE GENERIC
# OUTPUT bwa.bam 
# OUTPUT bwa.bam.bai 
# OUTPUT bwa.log 
# PARAMETER mode: "Data source" TYPE [ normal: " Illumina, 454, IonTorrent reads longer than 70 base pairs", pacbio: "PacBio subreads"] DEFAULT normal (Defining the type of reads will instruct the tool to use a predefined set of parameters optimized for that read type.)

# KM 11.11.2014

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads1.txt")
unzipIfGZipFile("reads2.txt")
unzipIfGZipFile("genome.txt")

# bwa
bwa.binary <- file.path(chipster.tools.path, "bwa", "bwa mem")
bwa.index.binary <- file.path(chipster.module.path, "shell", "check_bwa_index.sh")
command.start <- paste("bash -c '", bwa.binary)

# Do indexing
print("Indexing the genome...")
system("echo Indexing the genome... > bwa.log")
check.command <- paste ( bwa.index.binary, "genome.txt| tail -1 ")
genome.dir <- system(check.command, intern = TRUE)
bwa.genome <- file.path( genome.dir , "genome.txt")

mode.parameters <- ifelse(mode == "pacbio", "-x pacbio", "")

# command ending
command.end <- paste(bwa.genome, "reads1.txt reads2.txt 1> alignment.sam 2>> bwa.log'")

# run bwa alignment
bwa.command <- paste(command.start, mode.parameters, command.end)

echo.command <- paste("echo '", bwa.binary , mode.parameters, bwa.genome, "reads.txt ' > bwa.log" )
#stop(paste('CHIPSTER-NOTE: ', bwa.command))
system(echo.command)
system(bwa.command)
		
# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# convert sam to bam
system(paste(samtools.binary, "view -bS alignment.sam -o alignment.bam"))

# sort bam
system(paste(samtools.binary, "sort alignment.bam alignment.sorted"))

# index bam
system(paste(samtools.binary, "index alignment.sorted.bam"))

# rename result files
system("mv alignment.sorted.bam bwa.bam")
system("mv alignment.sorted.bam.bai bwa.bam.bai")
