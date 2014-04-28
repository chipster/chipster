# TOOL express.R: "Count reads per transcripts using eXpress" (Counts reads per transcript using the eXpress package. You can give as input fastq files, which will be aligned to transcriptome first with Bowtie2. Alternatively you can give as input a BAM file of aligned reads. Note that he BAM file needs to be sorted by read name, and exactly the same set of transcripts should be used for the alignment and eXpress. You should allow as many multimappings as possible. You can also allow mismatches, as eXpress performs error correction.)
# INPUT OPTIONAL alignment.bam: "BAM file containing transcript alignments, sorted by read name" TYPE BAM (BAM file containing reads aligned to transcripts.)
# INPUT OPTIONAL reads1.fq: "Read file 1" TYPE GENERIC (Reads in fastq format)
# INPUT OPTIONAL reads2.fq: "Read file 2 with mates in matching order" TYPE GENERIC (Read file 2 with mates in matching order.)
# OUTPUT OPTIONAL effective-counts-express.tsv: "Effective read counts per transcript from eXpress" (This file contains effective read counts per transcript from eXpress.)
# OUTPUT OPTIONAL full-express-output.tsv: "Express output table containing all the columns." (Express output table containing all the columns.)
# OUTPUT OPTIONAL transcriptome-alignment.bam: "Reads aligned to transcriptome" (Reads aligned to transcriptome using Bowtie2.)
# PARAMETER transcriptome: Transcriptome TYPE [humanrefseq: "Human RefSeq", mouserefseq: "Mouse RefSeq", ratrefseq: "Rat RefSeq", transcripts: transcripts] DEFAULT transcripts (Which transcriptome should be used for counting reads? Note that you should use the same set of transcripts as was used for the alignment.)
# PARAMETER fragmentlength: "Mean fragment length" TYPE INTEGER FROM 10 TO 10000 DEFAULT 200 (While the empirical distribution is estimated from paired-end reads on-the-fly, this value paramaterizes the prior distribution. If only single-end reads are available, this prior distribution is also used to determine the effective length.)


# check out if the fastq files are compressed and if so unzip them
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads1.fq")
unzipIfGZipFile("reads2.fq")

# set up Bowtie2 and SAMtools
bowtie.binary <- c(file.path(chipster.tools.path, "bowtie2", "bowtie2"))
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# read the transcriptome selection
transcriptome <- file.path(chipster.tools.path, "express_files", transcriptome)

# build the Bowtie2 command to align reads to the selected transcriptome (if starting from fastq files)
bowtie.command <- paste(bowtie.binary, "-a -X 800 -x", transcriptome, "-1 reads1.fq -2 reads2.fq |", samtools.binary, "view -Sb - > transcriptome-alignment.bam")

# run Bowtie2
system(bowtie.command)

# rename resulting alignment file for eXpress, return the original to the user
system("cp transcriptome-alignment.bam alignment.bam")

# set up eXpress
express.binary <- c(file.path(chipster.tools.path, "express", "express"))

# add the fasta ending to the transcriptome name (Bowtie2 index has the same basename)
transcriptome <- paste(transcriptome, ".fasta", sep="")

# command
command <- paste(express.binary, transcriptome, "alignment.bam")

# run
system(command)

# rename result file
system("cp results.xprs full-express-output.tsv")

# make a second output file containing only the effective counts
results <- read.table(file="results.xprs", sep="\t", header=T, quote="")
effcounts <- data.frame(target_id = results$target_id, eff_counts=round(results$eff_counts, digits=0))
write.table(effcounts, file="effective-counts-express.tsv", sep="\t", row.names=F, quote=F)




