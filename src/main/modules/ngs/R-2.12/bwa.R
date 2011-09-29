# TOOL bwa.R: "BWA for single end reads" (BWA aligns reads to genome, transcriptome, known miRNAs, etc.
# Results are sorted and indexed bam files, which are ready for viewing in the Chipster genome browser. 
# Note that this BWA tool uses public genomes provided by Chipster. If you would like to align reads against your own datasets, please use the tool \"BWA against own genomes\".)
# INPUT reads.txt: "Reads to align" TYPE GENERIC 
# OUTPUT bwa.bam 
# OUTPUT bwa.bam.bai 
# OUTPUT bwa.log 
# PARAMETER genome: "Genome or transcriptome" TYPE [hg19.fa: "Human genome (hg19\)", mm9.fa: "Mouse genome (mm9\)", rn4.fa: "Rat genome (rn4\)", mmu_miRB17mature.fa: "Mouse miRBase17"] DEFAULT mm9.fa (Genome or transcriptome that you would like to align your reads against.)
# PARAMETER algorithm: "Aligning algorithm type" TYPE [aln: "aln", bwasw: "BWA-SW"] DEFAULT aln (Alingment algorithm to be used. The default algorthm is aln. BWA-SW is designmed for long, over 300pb,  single-end reads with more errors)
# PARAMETER seed.length: "Seed length" TYPE INTEGER DEFAULT 32 ( Number of first nucleotides to be used as seed. If the seed length is longer than query sequences, then seeding will be disabled) 
# PARAMETER seed.edit:"maximum of differences in the seed" TYPE INTEGER DEFAULT 2 ( maximum differences in the seed )
# PARAMETER quality.format: "Quality value format used" TYPE [solexa1_3: "Illumina GA v1.3 or later", sanger: Sanger] DEFAULT sanger (Note that this parameter is taken into account only if you chose to apply the mismatch limit to the seed region. Are the quality values in the Sanger format (ASCII characters equal to the Phred quality plus 33\) or in the Illumina Genome Analyzer Pipeline v1.3 or later format (ASCII characters equal to the Phred quality plus 64\)? Please see the manual for details.)
# PARAMETER OPTIONAL num.gaps: "Maximum number of gap opens" TYPE INTEGER DEFAULT 1 ( Maximum number of gap opens for one read )
# PARAMETER OPTIONAL num.extensions: "Maximum number of gap extensions" TYPE INTEGER DEFAULT -1 (Maximum number of gap extensions, -1 for disabling long gaps )
# PARAMETER OPTIONAL gap.opening: "Gap open penalty " TYPE INTEGER DEFAULT 11 (Gap open penalty)
# PARAMETER OPTIONAL gap.extension: "Gap extension penalty " TYPE INTEGER DEFAULT 4 ( Gap extension penalty)
# PARAMETER OPTIONAL mismatch.penalty: "Mismatch penalty threshold" TYPE INTEGER DEFAULT 4 ( BWA will not search for suboptimal hits with a score lower than defined. )
# PARAMETER OPTIONAL disallow.gaps: "Maximum occurrences for extending a long deletion"  TYPE INTEGER DEFAULT 10 ( Maximum occurrences for extending a long deletion )
# PARAMETER OPTIONAL disallow.indel: "Disallow an indel within the given number of pb towards the ends"  TYPE INTEGER DEFAULT 5 ( do not put an indel within the defined value of bp towards the ends )
# PARAMETER OPTIONAL trim.threshold: "quality trimming threshold" TYPE INTEGER DEFAULT 0 ( quality threshold for read trimming down to 35bp )
# PARAMETER OPTIONAL barcode.length: "Barcode length"  TYPE INTEGER DEFAULT 0 ( Length of barcode starting from the 5 pime-end. The barcode of each read will be trimmed before mapping.)
# PARAMETER OPTIONAL alignment.no: "How many valid alignments are reported per read" TYPE  INTEGER DEFAULT 3 ( Maximum number of alignments to output in the XA tag for reads paired properly. If a read has more than the given amount of  hits, the XA tag will not be written)


# KM 24.8.2011

# bwa
bwa.binary <- file.path(chipster.tools.path, "bin", "bwa")
bwa.indexes <- file.path(chipster.tools.path, "bwa_indexes")
bwa.genome <- file.path(bwa.indexes, genome)
command.start <- paste("bash -c '", bwa.binary)


# mode specific parameters
quality.parameter <- ifelse(quality.format == "solexa1_3", "-I", "")
aln.mode.parameters <- paste("aln -t 2 -o", num.gaps, "-e", num.extensions, "-d", disallow.gaps, "-i" , disallow.indel , "-l" , seed.length , "-k" , seed.edit , "-O" , gap.opening , "-E" , gap.extension , "-q" , trim.threshold, "-B" , barcode.length , "-M" , mismatch.penalty , quality.parameter)
bwasw.mode.parameters <- paste("bwasw -t 2 -b", mismatch.penalty , "-q" , gap.opening , "-r" ,  gap.extension )
mode.parameters <- ifelse(algorithm == "aln", aln.mode.parameters, bwasw.mode.parameters)

# command ending
command.end <- paste(bwa.genome, "reads.txt 1> alignment.sai 2> bwa.log'")

# run bwa alignment
bwa.command <- paste(command.start, mode.parameters, command.end)
#stop(paste('CHIPSTER-NOTE: ', bwa.command))
system(bwa.command)


# sai to sam conversion
samse.parameters <- paste("samse -n", alignment.no )
samse.end <- paste(genome, "alignment.sai reads.txt 1> alignment.sam 2>samse.log'" )
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