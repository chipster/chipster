# TOOL bwa.R: "BWA for single end reads" (BWA aligns reads to genomes and transcriptomes. Results are sorted and indexed bam files, which are ready for viewing in the Chipster genome browser. 
# Note that this BWA tool uses publicly available genomes. If you would like to align reads against your own datasets, please use the tool \"BWA for single end reads and own genome\".)
# INPUT reads.txt: "Reads to align" TYPE GENERIC 
# OUTPUT bwa.bam 
# OUTPUT bwa.bam.bai 
# OUTPUT bwa.log 
# PARAMETER genome: "Genome or transcriptome" TYPE [hg19.fa: "Human genome (hg19\)", mm10.fa: "Mouse genome (mm10\)", mm9.fa: "Mouse genome (mm9\)", Rattus_norvegicus.Rnor_5.0.70.dna.toplevel.fa: "Rat genome (rn5\)", rn4.fa: "Rat genome (rn4\)", Canis_familiaris.CanFam3.1.71.dna.toplevel.fa: "Dog genome (Ensembl canFam3\)", Gasterosteus_aculeatus.BROADS1.71.dna.toplevel: "Gasterosteus aculeatus genome (BROADS1.71\)",  ovis_aries_texel.fa: "Sheep (oar3.1\)", Sus_scrofa.Sscrofa10.2.69.dna.toplevel.fa: "Pig (10.2\)", Drosophila_melanogaster.BDGP5.73.dna.toplevel.fa: "Drosophila_melanogaster (BDGP5.73\)", athaliana.TAIR10.fa: "A. thaliana genome (TAIR10\)", Arabidopsis_lyrata.v.1.0.16.fa: "A. Lyrata (1.0.16\)", Escherichia_coli_n1.GCA_000303635.1.18.dna.toplevel.fa: "E. coli gemone (GCA_000303635.1.18\)" ] DEFAULT hg19.fa (Genome or transcriptome that you would like to align your reads against.)
# PARAMETER seed.length: "Length of the seed region" TYPE INTEGER DEFAULT 32 (How many bases of the left, good quality part of the read should be used as the seed region. If the seed length is longer than the reads, the seeding will be disabled.) 
# PARAMETER seed.edit: "Maximum number of differences in the seed region" TYPE INTEGER DEFAULT 2 (Maximum number of differences such as mismatches or indels in the seed region.)
# PARAMETER total.edit: "Maximum edit distance for the whole read" TYPE DECIMAL DEFAULT 0.04 ( Maximum edit distance if the value is more than one. If the value is between 1 and 0 then it defines the fraction of missing alignments given 2% uniform base error rate. In the latter case, the maximum edit distance is automatically chosen for different read lengths. Corresponds to the command line parameter -n.)
# PARAMETER quality.format: "Quality value format used" TYPE [solexa1_3: "Illumina GA v1.3-1.5", sanger: Sanger] DEFAULT sanger (Note that this parameter is taken into account only if you chose to apply the mismatch limit to the seed region. Are the quality values in the Sanger format (ASCII characters equal to the Phred quality plus 33\) or in the Illumina Genome Analyzer Pipeline v1.3 or later format (ASCII characters equal to the Phred quality plus 64\)? Please see the manual for details. Corresponds to the command line parameter -I.)
# PARAMETER OPTIONAL num.gaps: "Maximum number of gaps" TYPE INTEGER DEFAULT 1 (Maximum number of gap openings for one read. Corresponds to the command line parameter -o)
# PARAMETER OPTIONAL num.extensions: "Maximum number of gap extensions" TYPE INTEGER DEFAULT -1 (Maximum number of gap extensions, -1 for disabling long gaps. Corresponds to the command line parameter -e)
# PARAMETER OPTIONAL gap.opening: "Gap opening penalty" TYPE INTEGER DEFAULT 11 (Gap opening penalty. Corresponds to the command line parameter -O )
# PARAMETER OPTIONAL gap.extension: "Gap extension penalty" TYPE INTEGER DEFAULT 4 (Gap extension penalty. Corresponds to the command line parameter -E)
# PARAMETER OPTIONAL mismatch.penalty: "Mismatch penalty threshold" TYPE INTEGER DEFAULT 3 (BWA will not search for suboptimal hits with a score lower than the alignment score minus this. Corresponds to the command line parameter -M)
# PARAMETER OPTIONAL disallow.gaps: "Disallow gaps in region"  TYPE INTEGER DEFAULT 16 (Disallow a long deletion within the given number of bp towards the 3\'-end. Corresponds to the command line parameter -d )
# PARAMETER OPTIONAL disallow.indel: "Disallow an indel within the given number of bp towards the ends"  TYPE INTEGER DEFAULT 5 (Do not put an indel within the defined value of bp towards the ends. Corresponds to the command line parameter -i)
# PARAMETER OPTIONAL trim.threshold: "Quality trimming threshold" TYPE INTEGER DEFAULT 0 (Quality threshold for read trimming down to 35bp. Corresponds to the command line parameter -q)
# PARAMETER OPTIONAL barcode.length: "Barcode length"  TYPE INTEGER DEFAULT 0 (Length of barcode starting from the 5 prime-end. The barcode of each read will be trimmed before mapping. Corresponds to the command line parameter -B)
# PARAMETER OPTIONAL alignment.no: "How many valid alignments are reported per read" TYPE  INTEGER DEFAULT 3 (Maximum number of alignments to report. Corresponds to the command line parameter bwa samse -n )

# KM 24.8.2011
# AMS 19.6.2012 Added unzippingcd 

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads.txt")

# bwa
bwa.binary <- file.path(chipster.tools.path, "bwa", "bwa")
bwa.indexes <- file.path(chipster.tools.path, "bwa_indexes")
bwa.genome <- file.path(bwa.indexes, genome)
command.start <- paste("bash -c '", bwa.binary)

# mode specific parameters
if (total.edit >= 1) {
	total.edit <- round(total.edit)
}

quality.parameter <- ifelse(quality.format == "solexa1_3", "-I", "")
mode.parameters <- paste("aln -t 2 -o", num.gaps, "-e", num.extensions, "-d", disallow.gaps, "-i" , disallow.indel , "-l" , seed.length , "-k" , seed.edit , "-O" , gap.opening , "-E" , gap.extension , "-q" , trim.threshold, "-B" , barcode.length , "-M" , mismatch.penalty , "-n" , total.edit , quality.parameter)

# command ending
command.end <- paste(bwa.genome, "reads.txt 1> alignment.sai 2>> bwa.log'")

# run bwa alignment
bwa.command <- paste(command.start, mode.parameters, command.end)

echo.command <- paste("echo '", bwa.binary , mode.parameters, bwa.genome, "reads.txt ' > bwa.log" )
#stop(paste('CHIPSTER-NOTE: ', bwa.command))
system(echo.command)
system(bwa.command)

#system ("pwd")
#system ("ls -l >> bwa.log")
# sai to sam conversion
samse.parameters <- paste("samse -n", alignment.no )
samse.end <- paste(bwa.genome, "alignment.sai reads.txt 1> alignment.sam 2>>bwa.log'" )
samse.command <- paste( command.start, samse.parameters , samse.end )
paste('CHIPSTER-NOTE: ', samse.command)
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
