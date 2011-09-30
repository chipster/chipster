# TOOL bwa-with-index-building.R: "BWA for single end reads and own genome" (BWA aligns reads to genome, transcriptome, known miRNAs, etc.
# Results are sorted and indexed bam files, which are ready for viewing in the Chipster genome browser. 
# Note that this BWA tool requires that you have impoted the refernce genome to chipster in fasta format. If you would like to align reads against your own datasets, please use the tool \"BWA against own genomes\".)
# INPUT reads.txt: "Reads to align" TYPE GENERIC 
# INPUT genome.txt: "Refence genome" TYPE GENERIC
# OUTPUT bwa.bam 
# OUTPUT bwa.bam.bai 
# OUTPUT bwa.log 
# PARAMETER algorithm: "Aligning algorithm type" TYPE [aln: "aln", bwasw: "BWA-SW"] DEFAULT aln (Alingment algorithm to be used. The default algorthm is aln. BWA-SW is designmed for long, over 300pb,  single-end reads with more errors)
# PARAMETER seed.length: "Seed length" TYPE INTEGER DEFAULT 32 ( -l : Number of first nucleotides to be used as seed. If the seed length is longer than query sequences, then seeding will be disabled) 
# PARAMETER seed.edit:"maximum of differences in the seed" TYPE INTEGER DEFAULT 2 ( -k: maximum differences in the seed )
# PARAMETER quality.format: "Quality value format used" TYPE [solexa1_3: "Illumina GA v1.3 or later", sanger: Sanger] DEFAULT sanger (Note that this parameter is taken into account only if you chose to apply the mismatch limit to the seed region. Are the quality values in the Sanger format (ASCII characters equal to the Phred quality plus 33\) or in the Illumina Genome Analyzer Pipeline v1.3 or later format (ASCII characters equal to the Phred quality plus 64\)? Please see the manual for details.)
# PARAMETER OPTIONAL num.gaps: "Maximum number of gap opens" TYPE INTEGER DEFAULT 1 ( -o : Maximum number of gap opens for one read )
# PARAMETER OPTIONAL num.extensions: "Maximum number of gap extensions" TYPE INTEGER DEFAULT -1 ( -e : Maximum number of gap extensions, -1 for disabling long gaps )
# PARAMETER OPTIONAL gap.opening: "gap open penalty " TYPE INTEGER DEFAULT 11 ( -O: Gap open penalty)
# PARAMETER OPTIONAL gap.extension: "gap extension penalty " TYPE INTEGER DEFAULT 4 ( -E: Gap extension penalty)
# PARAMETER OPTIONAL mismatch.penalty: "Mismatch penalty threshold" TYPE INTEGER DEFAULT 4 ( -M: BWA will not search for suboptimal hits with a score lower than defined. )
# PARAMETER OPTIONAL disallow.gaps: "Maximum occurrences for extending a long deletion"  TYPE INTEGER DEFAULT 10 ( -d : "Maximum occurrences for extending a long deletion )
# PARAMETER OPTIONAL disallow.indel: "Maximum occurrences for extending a long deletion"  TYPE INTEGER DEFAULT 5 ( -i : "do not put an indel within the defined value of bp towards the ends )
# PARAMETER OPTIONAL trim.threshold: "quality trimming threshold" TYPE INTEGER DEFAULT 0 (-q : "quality threshold for read trimming down to 35bp")
# PARAMETER OPTIONAL barcode.length: "Barcode length"  TYPE INTEGER DEFAULT 0 ( -B : Length of barcode starting from the 5’-end. The barcode of each read will be trimmed before mapping.)
# PARAMETER OPTIONAL alignment.no: "How many valid alignments are reported per read" TYPE  INTEGER DEFAULT 3 ( smase-n :Maximum number of alignments to output in the XA tag for reads paired properly. If a read has more than the given amount of  hits, the XA tag will not be written)

# KM 24.8.2011

# bwa
bwa.binary <- file.path(chipster.tools.path, "bin", "bwa")
bwa.index.binary <- file.path(chipster.tools.path, "bin", "check_bwa_index.csh")
bwa.indexes <- file.path(chipster.tools.path, "bwa_indexes")
command.start <- paste("bash -c '", bwa.binary)

check.command <- paste ( bwa.index.binary, "genome.txt| tail -1 ")
genome.dir <- system(check.command, intern = TRUE)
bwa.genome <- file.path( genome.dir , "genome.txt")
#stop(paste('CHIPSTER-NOTE: ', files))
# common parameters

# mode specific parameters
quality.parameter <- ifelse(quality.format == "solexa1_3", "-I", "")
aln.mode.parameters <- paste("aln -t 2 -o", num.gaps, "-e", num.extensions, "-d", disallow.gaps, "-i" , disallow.indel , "-l" , seed.length , "-k" , seed.edit , "-O" , gap.opening , "-E" , gap.extension , "-q" , trim.threshold, "-B" , barcode.length, "-M" , mismatch.penalty , quality.parameter)
bwasw.mode.parameters <- paste("bwasw -t 2 -b", mismatch.penalty , "-q" , gap.opening , "-r" ,  gap.extension )
mode.parameters <- ifelse(algorithm == "aln", aln.mode.parameters, bwasw.mode.parameters)

# command ending
command.end <- paste( bwa.genome , "reads.txt 1> alignment.sai 2>> bwa.log'")

# run bwa alignment
bwa.command <- paste(command.start, mode.parameters, command.end)
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