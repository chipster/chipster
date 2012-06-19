# TOOL bwa-with-index-building.R: "BWA for single end reads and own genome" (BWA aligns reads to genomes and transcriptomes. Unlike Bowtie, it supports gapped alignments. Results are sorted and indexed bam files. 
# Note that this BWA tool requires that you have imported the reference genome to Chipster in fasta format. If you would like to align reads against publicly available genomes, please use the tool \"BWA for single end reads\".)
# INPUT reads.txt: "Reads to align" TYPE GENERIC 
# INPUT genome.txt: "Reference genome" TYPE GENERIC
# OUTPUT bwa.bam 
# OUTPUT bwa.bam.bai 
# OUTPUT bwa.log 
# PARAMETER seed.length: "Length of the seed region" TYPE INTEGER DEFAULT 32 (How many bases of the left, good quality part of the read should be used as the seed region. If the seed length is longer than the reads, the seeding will be disabled. Corresponds to the command line parameter -l.) 
# PARAMETER seed.edit: "Maximum number of of differences in the seed region" TYPE INTEGER DEFAULT 2 (Maximum number of differences such as mismatches or indels in the seed region. Corresponds to the command line parameter -k.)
# PARAMETER total.edit: "Maximum edit distance for the whole read" TYPE DECIMAL DEFAULT 0.04 ( Maximum edit distance if the value is more than one. If the value is between 1 and 0 then it defines the fraction of missing alignments given 2% uniform base error rate. In the latter case, the maximum edit distance is automatically chosen for different read lengths. Corresponds to the command line parameter -n.)
# PARAMETER quality.format: "Quality value format used" TYPE [solexa1_3: "Illumina GA v1.3 or later", sanger: Sanger] DEFAULT sanger (Note that this parameter is taken into account only if you chose to apply the mismatch limit to the seed region. Are the quality values in the Sanger format (ASCII characters equal to the Phred quality plus 33\) or in the Illumina Genome Analyzer Pipeline v1.3 or later format (ASCII characters equal to the Phred quality plus 64\)? Please see the manual for details.)
# PARAMETER OPTIONAL num.gaps: "Maximum number of gaps" TYPE INTEGER DEFAULT 1 (Maximum number of gap openings for one read. Corresponds to the command line parameter -o.)
# PARAMETER OPTIONAL num.extensions: "Maximum number of gap extensions" TYPE INTEGER DEFAULT -1 (Maximum number of gap extensions, -1 for disabling long gaps. Corresponds to the command line parameter -e.)
# PARAMETER OPTIONAL gap.opening: "Gap opening penalty " TYPE INTEGER DEFAULT 11 (Gap opening penalty. Corresponds to the command line parameter -O.)
# PARAMETER OPTIONAL gap.extension: "Gap extension penalty " TYPE INTEGER DEFAULT 4 (Gap extension penalty. Corresponds to the command line parameter -E.)
# PARAMETER OPTIONAL mismatch.penalty: "Mismatch penalty threshold" TYPE INTEGER DEFAULT 4 (BWA will not search for suboptimal hits with a score lower than the alignment score minus this. Corresponds to the command line parameter -M.)
# PARAMETER OPTIONAL disallow.gaps: "Maximum occurrences for extending a long deletion"  TYPE INTEGER DEFAULT 10 (Maximum occurrences for extending a long deletion. Corresponds to the command line parameter -d.)
# PARAMETER OPTIONAL disallow.indel: "Disallow an indel within the given number of bp towards the ends"  TYPE INTEGER DEFAULT 5 (Do not put an indel within the defined value of bp towards the ends. Corresponds to the command line parameter -i.)
# PARAMETER OPTIONAL trim.threshold: "Quality trimming threshold" TYPE INTEGER DEFAULT 0 (Quality threshold for read trimming down to 35bp. Corresponds to the command line parameter -q.)
# PARAMETER OPTIONAL barcode.length: "Barcode length"  TYPE INTEGER DEFAULT 0 (Length of barcode starting from the 5 prime end. The barcode of each read will be trimmed before mapping. Corresponds to the command line parameter -B.)
# PARAMETER OPTIONAL alignment.no: "How many valid alignments are reported per read" TYPE  INTEGER DEFAULT 3 (Maximum number of alignments to report. Corresponds to the samse command line parameter -n.)

# KM 24.8.2011
# AMS 19.6.2012 Added unzipping

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads.txt")

# bwa settings
bwa.binary <- file.path(chipster.tools.path, "bwa", "bwa")
bwa.index.binary <- file.path(chipster.tools.path, "bwa", "check_bwa_index.sh")
command.start <- paste("bash -c '", bwa.binary)


# Do indexing
print("Indexing the genome...")
system("echo Indexing the genome... > bwa.log")
check.command <- paste ( bwa.index.binary, "genome.txt| tail -1 ")
genome.dir <- system(check.command, intern = TRUE)
bwa.genome <- file.path( genome.dir , "genome.txt")

# mode specific parameters
if (total.edit >= 1) {
	total.edit <- round(total.edit)
}
quality.parameter <- ifelse(quality.format == "solexa1_3", "-I", "")
mode.parameters <- paste("aln -t 2 -o", num.gaps, "-e", num.extensions, "-d", disallow.gaps, "-i" , disallow.indel , "-l" , seed.length , "-k" , seed.edit , "-O" , gap.opening , "-E" , gap.extension , "-q" , trim.threshold, "-B" , barcode.length, "-M" , mismatch.penalty , "-n" , total.edit , quality.parameter)


# command ending
command.end <- paste( bwa.genome , "reads.txt 1> alignment.sai 2>> bwa.log'")
system("echo Running the alignment with command: >> bwa.log")

# run bwa alignment
bwa.command <- paste(command.start, mode.parameters, command.end)
echo.command <- paste("echo '",bwa.command ,"'  >> bwa.log")
system(echo.command)
#stop(paste('CHIPSTER-NOTE: ', bwa.command))
system(bwa.command)


# sai to sam conversion
samse.parameters <- paste("samse -n", alignment.no )
samse.end <- paste(bwa.genome, "alignment.sai reads.txt 1> alignment.sam 2>> bwa.log'" )
samse.command <- paste( command.start, samse.parameters , samse.end )
system(samse.command)

		
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