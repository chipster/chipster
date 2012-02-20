# TOOL bwasw.R: "BWA-SW for single-end reads" (BWA aligns reads to genomes and transcriptomes. Unlike Bowtie, it allows gapped alignments. This BWA version uses the BWA-SW algorithm that is designed for long, over 300pb, low-quality reads. The alignment results are sorted and indexed bam files, which are ready for viewing in the Chipster genome browser. 
# Note that this BWA tool uses publicly available genomes. If you would like to align reads against your own datasets, please use the tool \"BWA against own genomes\".)
# INPUT reads.txt: "Reads to align" TYPE GENERIC 
# OUTPUT bwa.bam 
# OUTPUT bwa.bam.bai 
# OUTPUT bwa.log 
# PARAMETER genome: "Genome or transcriptome" TYPE [hg19.fa: "Human genome (hg19\)", mm9.fa: "Mouse genome (mm9\)", rn4.fa: "Rat genome (rn4\)"] DEFAULT mm9.fa (Genome or transcriptome that you would like to align your reads against.)
# PARAMETER OPTIONAL match.score: "Score of a match" TYPE INTEGER DEFAULT 1 (Score of a matching base. Corresponds to the command line parameter -a)
# PARAMETER OPTIONAL gap.opening: "Gap opening penalty " TYPE INTEGER DEFAULT 11 (Gap opening penalty. Corresponds to the command line parameter -q)
# PARAMETER OPTIONAL gap.extension: "Gap extension penalty " TYPE INTEGER DEFAULT 4 (Gap extension penalty. Corresponds to the command line parameter -r)
# PARAMETER OPTIONAL mismatch.penalty: "Mismatch penalty" TYPE INTEGER DEFAULT 3 (Mismatch penalty. Corresponds to the command line parameter -b)
# PARAMETER OPTIONAL band.width: "Band width in the banded alignment"  TYPE INTEGER DEFAULT 33 (Band width in the banded alignment. Corresponds to the command line parameter -w .)
# PARAMETER OPTIONAL min.score: "Minimum score threshold divided by the given value"  TYPE INTEGER DEFAULT 37 (Minimum score threshold divided by the given value. Corresponds to the command line parameter -T.)
# PARAMETER OPTIONAL threshold.coeff: "Coefficient for threshold adjustment" TYPE DECIMAL DEFAULT 5.5 (Coefficient for threshold adjustment according to query length. Given an l-long query, the threshold for a hit to be retained is a*max\{T,c*log\(l\)\}. Corresponds to the command line parameter -c )
# PARAMETER OPTIONAL z.best: "Z-best heuristics" TYPE INTEGER DEFAULT 1 (Z-best heuristics. Higher -z increases accuracy at the cost of speed. Corresponds to the command line parameter -z.)
# PARAMETER OPTIONAL sa.interval: "Maximum SA interval size" TYPE INTEGER DEFAULT 3 (Maximum SA interval size for initiating a seed. Higher value increases accuracy at the cost of speed. Corresponds to the command line parameter -s.)
# PARAMETER OPTIONAL min.support: "Reverse alignment limit" TYPE INTEGER DEFAULT 5 (Minimum number of seeds supporting the resultant alignment to skip reverse alignment. Corresponds to the command line parameter -N)
# PARAMETER OPTIONAL alignment.no: "How many valid alignments are reported per read" TYPE  INTEGER DEFAULT 3 (Maximum number of alignments to report. Corresponds to the command line parameter bwa samse -n )

# KM 24.8.2011

# bwa
bwa.binary <- file.path(chipster.tools.path, "bwa", "bwa")
bwa.indexes <- file.path(chipster.tools.path, "bwa_indexes")
bwa.genome <- file.path(bwa.indexes, genome)
command.start <- paste("bash -c '", bwa.binary)


# algorithm parameters
mode.parameters <- paste("bwasw -t 2 -b", mismatch.penalty , "-q" , gap.opening , "-r" , gap.extension , "-a" , match.score , "-w" , band.width , "-T" , min.score , "-c" , threshold.coeff , "-z" , z.best , "-s" , sa.interval , "-N" , min.support)

# command ending
echo.command <- paste("echo '", bwa.binary , mode.parameters, bwa.genome, "reads.txt ' > bwa.log" )
system(echo.command)

command.end <- paste(bwa.genome, "reads.txt 1> alignment.sai 2>> bwa.log'")

# run bwa alignment
bwa.command <- paste(command.start, mode.parameters, command.end)

#stop(paste('CHIPSTER-NOTE: ', bwa.command))
system(bwa.command)

#system ("pwd")
#system ("ls -l >> bwa.log")
# sai to sam conversion
samse.parameters <- paste("samse -n", alignment.no )
samse.end <- paste(bwa.genome, "alignment.sai reads.txt 1> alignment.sam 2>>bwa.log'" )
samse.command <- paste( command.start, samse.parameters , samse.end )

echo.command <- paste("echo '", bwa.binary , samse.parameters, bwa.genome, "alignment.sai reads.txt ' >> bwa.log" )
system(echo.command)

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