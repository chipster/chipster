# TOOL bwa-mem-with-index-building.R: "BWA MEM for single or paired end reads and own genome" (Aligns reads to genomes using the BWA MEM algorithm. If just one read file is given, then a single end analysis is run. If two read files are given, then mapping is done in paired end mode. Results are sorted and indexed BAM files, which are ready for viewing in the Chipster genome browser. 
# Note that this BWA tool requires that you have imported the reference genome to Chipster in fasta format. If you would like to align reads against publicly available genomes, please use the tool \"BWA MEM for single or paired end reads\".)
# INPUT genome.txt: "Reference genome" TYPE GENERIC
# INPUT reads.fastq: "Read file 1" TYPE GENERIC 
# INPUT OPTIONAL reads2.fastq: "Read file 2 for paired end reads" TYPE GENERIC 
# OUTPUT bwa.bam 
# OUTPUT bwa.bam.bai 
# OUTPUT bwa.log 
# PARAMETER OPTIONAL minseedlen: "Minimum seed length" TYPE INTEGER DEFAULT 19 (Matches shorter than this will be missed when looking for maximal exact matches or MEMs in the first alignment phase.)
# PARAMETER OPTIONAL bandwith: "Maximum gap length" TYPE INTEGER DEFAULT 100 (Gaps longer than this will not be found. Note also scoring matrix and hit length affect the maximum gap length, in addition to this band width parameter.)
# PARAMETER OPTIONAL matchscore: "Match score" TYPE INTEGER DEFAULT 1 (Score for a matching base.)
# PARAMETER OPTIONAL mismatchscore: "Mismatch penalty" TYPE INTEGER DEFAULT 4 (Penalty for a mismatching base.)
# PARAMETER OPTIONAL gapopen: "Gap opening penalty" TYPE INTEGER DEFAULT 6 (Gap opening penalty.)
# PARAMETER OPTIONAL gapextension: "Gap extension penalty" TYPE INTEGER DEFAULT 1 (Gap extension penalty.)
# PARAMETER OPTIONAL clippenalty: "Penalty for end clipping" TYPE INTEGER DEFAULT 5 (Penalty for 5\'- and 3\'-end clipping. When performing the Smith-Waterman extension of the seed alignments, BWA-MEM keeps track of the best score reaching the end of the read. If this score is larger than the best SW score minus the clipping penalty, clipping will not be applied.) 
# PARAMETER OPTIONAL rgid: "Read group dentifier" TYPE STRING (Read group identifier. If you want to add the read group line in the BAM file, you have to give this information.)
# PARAMETER OPTIONAL rgsm: "Sample name for read group" TYPE STRING (The name of the sample sequenced in this read group. Note that you have to fill in also the read group identifier parameter for the read group information to appear in the BAM file.)
# PARAMETER OPTIONAL rgpl: "Platform for read group" TYPE [ none: "Not defined", ILLUMINA, SOLID, LS454, HELICOS, PACBIO] DEFAULT none (Platform\/technology used to produce the read. Note that you have to fill in also the read group identifier parameter for the read group information to appear in the BAM file.)
# PARAMETER OPTIONAL rglb: "Library identifier for read group" TYPE STRING (DNA preparation library identifier. The Mark Duplicates tool uses this field to determine which read groups might contain molecular duplicates, in case the same DNA library was sequenced on multiple lanes. Note that you have to fill in also the read group identifier parameter for the read group information to appear in the BAM file.)

# KM 1.9.2015


# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads.fastq")
unzipIfGZipFile("genome.txt")

pe_check.command <- ("ls reads2.fastq 2>/dev/null | wc -l")

paired.end <- system(pe_check.command, intern = TRUE )

if ( paired.end == 1 ){
	unzipIfGZipFile("reads2.fastq")
}


# bwa
bwa.binary <- file.path(chipster.tools.path, "bwa", "bwa mem")
bwa.index.binary <- file.path(chipster.module.path, "shell", "check_bwa_index.sh")

# Do indexing
print("Indexing the genome...")
system("echo Indexing the genome... > bwa.log")
check.command <- paste ( bwa.index.binary, "genome.txt| tail -1 ")
genome.dir <- system(check.command, intern = TRUE)
bwa.genome <- file.path( genome.dir , "genome.txt")

# bwa
command.start <-(bwa.binary)

bwa.parameters <- paste("-M", "-k", minseedlen, "-w", bandwith, "-A", matchscore, "-B", mismatchscore, "-O", gapopen, "-E", gapextension, "-L",  clippenalty )

#Read group definition
if ( nchar(rgid) > 0 ){
	rg.string <-("'@RG")
	rg.string <- paste(rg.string, "\\tID:", rgid, sep="")
	if ( nchar(rgsm) > 0 ){
		rg.string <- paste(rg.string, "\\tSM:", rgsm, sep="")
	}
	if ( rgpl != "none" ){
		rg.string <- paste(rg.string, "\\tPL:", rgpl, sep="")
	}
	if ( nchar(rglb) > 0 ){
		rg.string <- paste(rg.string, "\\tLB:", rglb, sep="")
	}
	rg.string <- paste(rg.string, "'", sep="")
	bwa.parameters <- paste(bwa.parameters,  "-R", rg.string )
}

# command ending
command.end <- paste(bwa.genome, "reads.fastq 1> alignment.sam 2>> bwa.log")
if ( paired.end == 1 ){
	command.end <- paste(bwa.genome, "reads.fastq reads2.fastq 1> alignment.sam 2>> bwa.log")
}	

bwa.command <- paste(command.start, bwa.parameters, command.end)

echo.command <- paste("echo '", bwa.binary, bwa.parameters, bwa.genome, "reads.fastq ' > bwa.log" )
if ( paired.end == 1 ){
	echo.command <- paste("echo '", bwa.binary, bwa.parameters, bwa.genome, "reads.fastq reads2.fastq ' > bwa.log" )
}

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

# Handle output names
#
source(file.path(chipster.common.path, "tool-utils.R"))

# read input names
inputnames <- read_input_definitions()

# Determine base name
basename <- strip_name(inputnames$reads.fastq)

# Make a matrix of output names
outputnames <- matrix(NA, nrow=2, ncol=2)
outputnames[1,] <- c("bwa.bam", paste(basename, ".bam", sep =""))
outputnames[2,] <- c("bwa.bam.bai", paste(basename, ".bam.bai", sep =""))

# Write output definitions file
write_output_definitions(outputnames)
