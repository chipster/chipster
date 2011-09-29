# TOOL bwa-paired-end-with-index-building.R: "BWA for paired-end reads and own genome" (BWA aligns reads to genome, transcriptome, known miRNAs, etc.
# Results are sorted and indexed bam files, which are ready for viewing in the Chipster genome browser. 
# Note that this BWA tool uses public genomes provided by Chipster. If you would like to align reads against your own datasets, please use the tool \"BWA against own genomes\".)
# INPUT reads1.txt: "Paired-end read set 1 to align" TYPE GENERIC 
# INPUT reads2.txt: "Paired-end read set 2 to align" TYPE GENERIC 
# INPUT genome.txt: "Refence genome" TYPE GENERIC
# OUTPUT bwa.bam 
# OUTPUT bwa.bam.bai 
# OUTPUT bwa.log 
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
# PARAMETER OPTIONAL alignment.no: "maximum hits to output for paired reads" TYPE INTEGER DEFAULT 3 ( Maximum number of alignments to output in the XA tag for reads paired properly. If a read has more than the given amount of  hits, the XA tag will not be written)
# PARAMETER OPTIONAL max.discordant: "maximum hits to output for discordant pairs" TYPE INTEGER DEFAULT 10 ( Maximum number of alignments to output in the XA tag for disconcordant read pairs, excluding singletons. If a read has more than INT hits, the XA tag will not be written.) 
# PARAMETER OPTIONAL max.insert: "maximum insert size" TYPE INTEGER DEFAULT 500 (Maximum insert size for a read pair to be considered being mapped properly. This option is only used when there are not enough good alignment to infer the distribution of insert sizes.)
# PARAMETER OPTIONAL max.occurrence: "maximum occurrences for one end" TYPE INTEGER DEFAULT 100000 (Maximum occurrences of a read for pairing. A read with more occurrneces will be treated as a single-end read. Reducing this parameter helps faster pairing. The Deafaul vaule is 100000. For reads shorter than 30bp, applying a smaller velue is recommended to get a sensible speed at the cost of pairing accuracy.)


# KM 26.8.2011

# bwa
bwa.binary <- file.path(chipster.tools.path, "bin", "bwa")
bwa.index.binary <- file.path(chipster.tools.path, "bin", "check_bwa_index.csh")
bwa.indexes <- file.path(chipster.tools.path, "bwa_indexes/tmp")

command.start <- paste("bash -c '", bwa.binary)

check.command <- paste ( bwa.index.binary, "genome.txt| tail -1 ")
genome.dir <- system(check.command, intern = TRUE)
bwa.genome <- file.path( genome.dir , "genome.txt")

#stop(paste('CHIPSTER-NOTE: ', bwa.genome))
# common parameters

# mode specific parameters
quality.parameter <- ifelse(quality.format == "solexa1_3", "-I", "")
mode.parameters <- paste("aln -t 2 -o", num.gaps, "-e", num.extensions, "-d", disallow.gaps, "-i" , disallow.indel , "-l" , seed.length , "-k" , seed.edit , "-O" , gap.opening , "-E" , gap.extension , "-q" , trim.threshold, "-B" , barcode.length , "-M" , mismatch.penalty , quality.parameter)

# command ending
command.end <- paste(bwa.genome, "reads1.txt 1> alignment1.sai 2> bwa.log'")
command.end <- paste(bwa.genome, "reads2.txt 1> alignment2.sai 2> bwa.log'")


# run bwa alignment
bwa.command <- paste(command.start, mode.parameters, command.end)
#stop(paste('CHIPSTER-NOTE: ', bwa.command))
system(bwa.command)


# sai to sam conversion
sampe.parameters <- paste("sampe -n", alignment.no, "-a", max.insert, "-o" , max.occurrence , "-N" , max.discordant )
sampe.end <- paste(bwa.genome, "alignment1.sai alignment2.sai reads1.txt reads1.txt 1> alignment.sam 2>samse.log'" )
sampe.command <- paste( command.start, sampe.parameters , sampe.end )
system(sampe.command)

		
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