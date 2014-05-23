# TOOL express.R: "Count reads per transcripts using eXpress" (Counts reads per transcripts using the eXpress package. This tool takes as input either paired end or single end reads in FASTQ files, and transcripts in a multi-fasta file. It aligns the reads to the transcripts using Bowtie2, and you should allow as many multimappings as possible. Note that Bowtie2 can be very slow if many alignments are allowed per read. Both the read and transcript file can be compressed.)
# INPUT transcripts.fasta: "Fasta file containing transcript sequences" TYPE GENERIC (Fasta file containing transcript sequences.)
# INPUT reads1.fq: "Read file 1" TYPE GENERIC (Reads in fastq format)
# INPUT OPTIONAL reads2.fq: "Read file 2 with mates in matching order" TYPE GENERIC (Read file 2 with mates in matching order.)
# OUTPUT OPTIONAL effective-counts-express.tsv: "Effective read counts per transcript from eXpress" (This file contains effective read counts per transcript from eXpress.)
# OUTPUT OPTIONAL full-express-output.tsv: "Express output table containing all the columns." (Express output table containing all the columns.)
# PARAMETER phred.scale: "Quality scale used in the fastq file" TYPE [phred33: "phred + 33", phred64: "phred + 64"] DEFAULT phred33 (Quality scale used in the fastq file.)
# PARAMETER OPTIONAL alignment.no: "How many valid alignments should Bowtie2 search" TYPE [20: "20", 100: "100", 1000: "1000"] DEFAULT 20 (How many alignments should Bowtie2 search per read.)
# PARAMETER OPTIONAL fragmentlength: "Mean fragment length" TYPE INTEGER FROM 10 TO 10000 DEFAULT 200 (While the empirical distribution is estimated from paired-end reads on-the-fly, this value paramaterizes the prior distribution. If only single-end reads are available, this prior distribution is also used to determine the effective length.)
# PARAMETER OPTIONAL fragmentlengthstdev: "Fragment length standard deviation" TYPE INTEGER FROM 5 TO 200 DEFAULT 60 (While the empirical distribution is estimated from paired-end reads on-the-fly, this value paramaterizes the prior distribution. If only single-end reads are available, this prior distribution is also used to determine the effective length.)

# EK 28.4.2014
# EK 23.5.2014 updated Bowtie2 command, added identifier sorting and piping to eXpress

# check out if the fastq files are compressed and if so unzip them
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads1.fq")
unzipIfGZipFile("reads2.fq")
unzipIfGZipFile("transcripts.fasta")

# set up Bowtie2 tools and eXpress
bowtie.binary <- c(file.path(chipster.tools.path, "bowtie2", "bowtie2"))
bowtie.build.binary <- c(file.path(chipster.tools.path, "bowtie2", "bowtie2-build"))
express.binary <- c(file.path(chipster.tools.path, "express", "express"))

# make Bowtie index for the transcripts
bowtie.build.command <- paste(bowtie.build.binary, "-offrate=1 -f transcripts.fasta transcriptome")
system(bowtie.build.command)

# collect Bowtie paramters
bowtie.params <- paste("--rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4 -p", chipster.threads.max, "-k", alignment.no)
if(phred.scale == "phred64"){
	bowtie.params <- paste(bowtie.params, "--phred64")
}

if (file.exists("reads2.fq")){
	# Paired end reads
	command <- paste(bowtie.binary, bowtie.params, "--no-discordant --no-mixed -x transcriptome -1 reads1.fq -2 reads2.fq |", express.binary, "-m", fragmentlength, "-s", fragmentlengthstdev, "transcripts.fasta")
} else{
	# Single end reads
	command <- paste(bowtie.binary, bowtie.params, "-x transcriptome -U reads1.fq |", express.binary, "-m", fragmentlength, "-s", fragmentlengthstdev, "transcripts.fasta")
}

# for testing
# stop(paste('CHIPSTER-NOTE: ', bowtie.command))

# run
system(command)

# rename result file
system("cp results.xprs full-express-output.tsv")

# take the title row and put it to a new file
system("head -1 results.xprs > id_sorted")

# sort the remaining rows by the identifier column and add it to the file. Sorting by the identifier is necessary so that count files for individual samples can be combined to a count table later.
system("tail -n +2 results.xprs | sort -k 2,2 >> id_sorted")

# make a second output file containing only the effective counts
results <- read.table(file="id_sorted", sep="\t", header=T, quote="")
effcounts <- data.frame(target_id = results$target_id, eff_counts=round(results$eff_counts, digits=0))
write.table(effcounts, file="effective-counts-express.tsv", sep="\t", row.names=F, quote=F)




