# TOOL bowtie.R: "Bowtie for single end reads" (Bowtie aligns reads to genome, transcriptome, known miRNAs, etc. There are two modes: mismatches are considered either throughout the read, or only in the user-defined left part of the read. In the latter case also quality values are taken into account. Results are sorted and indexed bam files, which are ready for viewing in the Chipster genome browser. 
# Note that this Bowtie tool uses publicly available genomes. If you would like to align reads against your own datasets, please use the tool \"Bowtie for single end reads and own genome\".)
# INPUT reads.txt: "Reads to align" TYPE GENERIC 
# OUTPUT bowtie.bam 
# OUTPUT bowtie.bam.bai 
# OUTPUT bowtie.log 
# OUTPUT OPTIONAL unaligned-reads.fastq 
# OUTPUT OPTIONAL multireads.fastq 
# PARAMETER genome: "Genome or transcriptome" TYPE [hg19: "Human genome (hg19\)", mm9: "Mouse genome (mm9\)", mm10: "Mouse genome (mm10\)", Rattus_norvegicus.Rnor_5.0.70.dna.toplevel: "Rat genome (rn5\)", rn4: "Rat genome (rn4\)", Halorubrum_lacusprofundi_ATCC_49239: "Halorubrum lacusprofundi ATCC 49239 genome", Canis_familiaris.CanFam3.1.71.dna.toplevel: "Dog genome (Ensembl canFam3\)", Gasterosteus_aculeatus.BROADS1.71.dna.toplevel: "Gasterosteus aculeatus genome (BROADS1.71\)", Drosophila_melanogaster.BDGP5.73.dna.toplevel: "Drosophila_melanogaster (BDGP5.73\)", athaliana.TAIR10: "A. thaliana genome (TAIR10\)", ovis_aries_texel: "Sheep genome (oar3.1\)", miRBase19_Rattus_norvegicus: "Rat miRBase 19", miRBase19_Mus_musculus: "Mouse miRBase 19", miRBase19_Homo_sapiens: "Human miRBase19", Escherichia_coli_n1.GCA_000303635.1.18.dna.toplevel: "E. coli gemone (GCA_000303635.1.18\)"] DEFAULT hg19 (Genome or transcriptome that you would like to align your reads against.)
# PARAMETER max.mismatches: "Number of mismatches allowed" TYPE [0, 1, 2, 3] DEFAULT 2 (How many mismatches are the alignments allowed to have?)
# PARAMETER limit.to.seed: "Consider mismatches only in the seed region" TYPE [yes, no] DEFAULT no (Should the mismatch limit be applied only to the left, good quality part of the read? You can define the length of this seed region with the next parameter.)
# PARAMETER seed: "Length of the seed region" TYPE INTEGER FROM 5 TO 50 DEFAULT 28 (If you have chosen to apply the mismatch limit only to the left, good quality part of the read, how many bases should be considered? The minimum length of seed region is 5.)
# PARAMETER quality: "Allowed total of mismatch qualities" TYPE INTEGER FROM 10 TO 100 DEFAULT 70 (What is the maximum permitted total of quality values of ALL mismatch positions throughout the read (not just in the seed region\)? Note that this parameter is taken into account only if you have chosen to apply the mismatch limit to the seed region.)
# PARAMETER quality.format: "Quality value format used" TYPE [solexa1_3: "Illumina GA v1.3-1.5", sanger: Sanger] DEFAULT sanger (Note that this parameter is taken into account only if you chose to apply the mismatch limit to the seed region. Are the quality values in the Sanger format (ASCII characters equal to the Phred quality plus 33\) or in the Illumina Genome Analyzer Pipeline v1.3 or later format (ASCII characters equal to the Phred quality plus 64\)? Please see the manual for details.)
# PARAMETER OPTIONAL multiread: "How many best category hits is the read allowed to have" TYPE [1, 2, 1000000: "no limit"] DEFAULT 1000000 (If you select 1, only those reads which have a unique best category hit will be reported. Note that the read can still have other, lower category hits even if they are not reported. The hit categories are based on the number of mismatches.)
# PARAMETER OPTIONAL alignment.no: "How many valid alignments are reported per read" TYPE [1, 2, 3] DEFAULT 1 (If there are several, equally good alignments, how many should be reported?)
# PARAMETER OPTIONAL multiread.file: "Put multireads to a separate file" TYPE [yes, no] DEFAULT no (If you chose not to have alignments for reads which map to multiple positions, would you like to store these reads to a separate fastq file?)
# PARAMETER OPTIONAL unaligned.file: "Put unaligned reads to a separate file" TYPE [yes, no] DEFAULT no (Would you like to store unaligned reads to a new fastq file? Note that also multireads will be added to this file, unless you asked them to be put to a separate file.)

# EK 22.3.2011
# AMS 19.6.2012 Added unzipping

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads.txt")

# bowtie
bowtie.binary <- c(file.path(chipster.tools.path, "bowtie", "bowtie"))
command.start <- paste("bash -c '", bowtie.binary)

# common parameters
common.parameters <- paste("-q --best -S --strata", "-m", multiread, "-k", alignment.no)

# mode specific parameters
quality.parameter <- ifelse(quality.format == "solexa1_3", "--solexa1.3-quals", "")
n.mode.parameters <- paste("-n", max.mismatches, "-l", seed, "-e", quality, quality.parameter)
v.mode.parameters <- paste("-v", max.mismatches)
mode.parameters <- ifelse(limit.to.seed == "yes", n.mode.parameters, v.mode.parameters)

# output parameters
unaligned.output <- ifelse(unaligned.file == "yes", "--un unaligned-reads.fastq", "")
multiread.output <- ifelse(multiread.file == "yes", "--max multireads.fastq", "")
output.parameters <- paste(unaligned.output, multiread.output)

# command ending
command.end <- paste(genome, "reads.txt 1> alignment.sam 2> bowtie.log'")

# run bowtie
bowtie.command <- paste(command.start, common.parameters, quality.parameter, mode.parameters, output.parameters, command.end)
#stop(paste('CHIPSTER-NOTE: ', bowtie.command))
system(bowtie.command)


# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# convert sam to bam, removes unaligned reads
system(paste(samtools.binary, "view -bS -q 1 alignment.sam -o alignment.bam"))

# sort bam
system(paste(samtools.binary, "sort alignment.bam alignment.sorted"))

# index bam
system(paste(samtools.binary, "index alignment.sorted.bam"))

# rename result files
system("mv alignment.sorted.bam bowtie.bam")
system("mv alignment.sorted.bam.bai bowtie.bam.bai")
