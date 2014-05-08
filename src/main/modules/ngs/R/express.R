# TOOL express.R: "Count reads per transcripts using eXpress" (Counts reads per transcripts using the eXpress package. This tool takes as input either paired end or single end reads in fastq files, and transcripts in a multi-fasta file. It aligns the reads to the transcripts using Bowtie2, and you should allow as many multimappings as possible. Note that Bowtie2 can be very slow if many alignments are allowed per read. Both the read and transcript file can be compressed.)
# INPUT transcripts.fasta: "Fasta file containing transcript sequences" TYPE GENERIC (Fasta file containing transcript sequences.)
# INPUT reads1.fq: "Read file 1" TYPE GENERIC (Reads in fastq format)
# INPUT OPTIONAL reads2.fq: "Read file 2 with mates in matching order" TYPE GENERIC (Read file 2 with mates in matching order.)
# OUTPUT OPTIONAL effective-counts-express.tsv: "Effective read counts per transcript from eXpress" (This file contains effective read counts per transcript from eXpress.)
# OUTPUT OPTIONAL full-express-output.tsv: "Express output table containing all the columns." (Express output table containing all the columns.)
# PARAMETER phred.scale: "Quality scale used in the fastq file" TYPE [phred33: "phred + 33", phred64: "phred + 64"] DEFAULT phred33 (Quality scale used in the fastq file.)
# PARAMETER OPTIONAL alignment.no: "How many valid alignments are reported per read" TYPE [20: "20", 100: "100", 1000: "1000", 6: "all alignments"] DEFAULT 20 (How many alignments should Bowtie report per read.)
# PARAMETER OPTIONAL fragmentlength: "Mean fragment length" TYPE INTEGER FROM 10 TO 10000 DEFAULT 200 (While the empirical distribution is estimated from paired-end reads on-the-fly, this value paramaterizes the prior distribution. If only single-end reads are available, this prior distribution is also used to determine the effective length.)
# PARAMETER OPTIONAL fragmentlengthstdev: "Fragment length standard deviation" TYPE INTEGER FROM 5 TO 200 DEFAULT 60 (While the empirical distribution is estimated from paired-end reads on-the-fly, this value paramaterizes the prior distribution. If only single-end reads are available, this prior distribution is also used to determine the effective length.)

# EK 28.4.2014

# check out if the fastq files are compressed and if so unzip them
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads1.fq")
unzipIfGZipFile("reads2.fq")
unzipIfGZipFile("transcripts.fasta")

# set up Bowtie2 tools and SAMtools
bowtie.binary <- c(file.path(chipster.tools.path, "bowtie2", "bowtie2"))
bowtie.build.binary <- c(file.path(chipster.tools.path, "bowtie2", "bowtie2-build"))
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# make Bowtie index for the transcripts
bowtie.build.command <- paste(bowtie.build.binary, "-offrate=1 -f transcripts.fasta transcriptome")
system(bowtie.build.command)

# collect Bowtie paramters
bowtie.params <- paste("-a -X 800 -p", chipster.threads.max)
if(phred.scale == "phred64"){
	bowtie.params <- paste(bowtie.params, "--phred64")
}
 # how many alignments should Bowtie report
if ( alignment.no==6){
	bowtie.params <- paste(bowtie.params, "-a")
}
if ( alignment.no>6){
	bowtie.params <- paste(bowtie.params, "-k", alignment.no )
}


if (file.exists("reads2.fq")){
	# Paired end reads
	bowtie.command <- paste(bowtie.binary, bowtie.params, "-x transcriptome -1 reads1.fq -2 reads2.fq |", samtools.binary, "view -Sb - > alignment.bam")
	
} else{
	# Single end reads
	bowtie.command <- paste(bowtie.binary, bowtie.params, "-x transcriptome -U reads1.fq |", samtools.binary, "view -Sb - > alignment.bam")
}

# run Bowtie2
system(bowtie.command)

# for testing
# stop(paste('CHIPSTER-NOTE: ', bowtie.command))

# set up eXpress
express.binary <- c(file.path(chipster.tools.path, "express", "express"))

# command
command <- paste(express.binary, "-m", fragmentlength, "-s", fragmentlengthstdev, "transcripts.fasta alignment.bam")

# run
system(command)

# rename result file
system("cp results.xprs full-express-output.tsv")

# make a second output file containing only the effective counts
results <- read.table(file="results.xprs", sep="\t", header=T, quote="")
effcounts <- data.frame(target_id = results$target_id, eff_counts=round(results$eff_counts, digits=0))
write.table(effcounts, file="effective-counts-express.tsv", sep="\t", row.names=F, quote=F)




