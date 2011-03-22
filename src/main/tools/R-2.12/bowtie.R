# TOOL "NGS" / bowtie.sadl: Bowtie (Alignment of reads to genome, transcriptome, known miRNAs, etc. Bowtie can be run with two modes: mismatches are either considered all along the alignment and quality values are not taken into account, or mismatches are considered only in the user-defined left part of the read and quality values are taken into account. Chipster's Bowtie tool gives the alignment results as sorted and indexed bam files, which are ready for viewing the Chipster genome browser. Note that this Bowtie tool uses genome indeces provided by us. If you would like to align reads against your own datasets, please use the tool "Bowtie with self-made indeces" which will apeear soon. ) 
# INPUT reads.txt: "Reads to align" TYPE GENERIC
# OUTPUT bowtie.bam
# OUTPUT bowtie.bam.bai
# OUTPUT bowtie.log
# OUTPUT OPTIONAL unaligned-reads.fastq
# OUTPUT OPTIONAL multireads.fastq
# PARAMETER genome: "Genome or transcriptome" TYPE [m_musculus_ncbi37: "m_musculus_ncbi37", hg19: "hg19", mm9: "mm9", rn4: "rn4"] DEFAULT m_musculus_ncbi37 (Genome or transcriptome that you would to align your reads against.)
# PARAMETER max.mismatches: "Number of mismatches allowed" TYPE [0, 1, 2, 3] DEFAULT 2 (How many mismatches are the alignments allowed to have?)
# PARAMETER limit.to.seed: "Consider mismatches only in the seed region" TYPE [yes, no] DEFAULT no (Should the mismatch limit be applied only to the left, good quality part of the read? You can define the lenght of this seed region with the next parameter.)
# PARAMETER seed: "Allowed total of mismatch qualities" TYPE INTEGER FROM 10 TO 100 DEFAULT 70 (What is the maximum permitted total of quality values at ALL mismatch positions throughoutthe alignment, not just in the seed region. Note that this parameter is taken into account only if you have chosen to apply the mismatch limit to the seed region.)
# PARAMETER quality: "Lenght of the seed region" TYPE INTEGER FROM 5 TO 50 DEFAULT 28 (If you have chosen to apply the mismatch limit only to the left, good quality part of the read, how many bases should be considered? The minimum length of seed region is 5.)
# PARAMETER quality.format: "Quality value format used" TYPE [solexa1_3: "Illumina GA v1.3 or later", sanger: "Sanger"] FROM 1 TO 2 DEFAULT sanger (Are the quality values in the Sanger format (ASCII characters equal to the Phred quality plus 33) or in the Illumina Genome Analyzer Pipeline v1.3 or later format (ASCII characters equal to the Phred quality plus 64)).
# PARAMETER multiread: "How many places is a read allowed to map to" TYPE [1, 2, 3] DEFAULT 1 (If you want to have alignments only for uniquely mapping reads, select 1.)
# PARAMETER alignment.no: "How many valid alignments can be reported per read" TYPE [1, 2, 3] DEFAULT 1 (If there are several, equally good alignments, how many should be reported?)
# PARAMETER unaligned.file: "Put unaligned reads to a separate file" TYPE [yes, no] DEFAULT no (Would you like to store unaligned reads to a new fastq file?)
# PARAMETER multiread.file: "Put reads mapping to several positions to a separate file" TYPE [yes, no] DEFAULT no (If you chose not to have alignments for reads which map to multiple positions, would you like to store these reads to a separate fastq file?)


# run bowtie
bowtie.binary <- c(file.path(chipster.tools.path, "bowtie", "bowtie"))
bowtie.parameters <- "-q --best -S --strata"
#if limit.to.seed == yes,
bowtie.command <- paste("bash -c '", bowtie.binary, bowtie.parameters, "-m", multiread, "-k", alignment.no, "-n", max.mismatches, -"l", seed, genome, "reads.txt 1> alignment.sam 2> bowtie.log'")

#if limit.to.seed == no,
#bowtie.command <- paste("bash -c '", bowtie.binary, bowtie.parameters, "-m", multiread, "-k", alignment.no, "-v", max.mismatches, genome, "reads.txt 1> alignment.sam 2> bowtie.log'")

#unaligned.file = "yes"
#ifelse(par == "yes", "--un unaligned-reads.fastq", "")
#[1] "--un filename"
#unaligned.file = "no"
#ifelse(par == "yes", "--un unaligned-reads.fastq", "")

#multiread.file = "yes"
#ifelse(par == "yes", "--max multireads.fastq", "")
#[1] "--un filename"
#multiread.file = "no"
#ifelse(par == "yes", "--max multireads.fastq", "")

system(bowtie.command)

# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# convert sam to bam
system(paste(samtools.binary, "view -bS alignment.sam -o alignment.bam"))

# sort bam
system(paste(samtools.binary, "sort alignment.bam alignment.sorted"))

# index bam
system(paste(samtools.binary, "index alignment.sorted.bam"))

# rename result files
system("mv alignment.sorted.bam alignment-bowtie.bam")
system("mv alignment.sorted.bam.bai alignment-bowtie.bam.bai")