# TOOL tophat.R: "TopHat" (TopHat is a program that aligns RNA-Seq reads to a genome in order to identify exon-exon splice junctions.
# TopHat was designed to work with reads produced by the Illumina Genome Analyzer, although users have been successful in using TopHat with reads from other technologies. In TopHat 1.1.0, we began supporting Applied Biosystems' Colorspace format. The software is optimized for reads 75bp or longer.
# Currently, TopHat does not allow short (fewer than a few nucleotides) insertions and deletions in the alignments it reports. Support for insertions and deletions will eventually be added.
# Finally, mixing paired- and single- end reads together is not supported.(
# INPUT reads.txt: "Reads to align" TYPE GENERIC 
# OUTPUT tophat.bam 
# OUTPUT tophat.bam.bai 
# OUTPUT tophat.log 

############################################################
#                                                          #
# Alignment tool that allows identification of exon-exon   #
# splice junctions using the TopHat program. TopHat uses   #
# the same index files as Bowtie and will look int the     #
# .../bowtie/indexes/ folder for those.                    #
#                                                          #
# MG, 23.6.2011                                            #
# development version                                      #
#                                                          #
############################################################


# Setting up TopHat
tophat.binary <- c(file.path(chipster.tools.path, "tophat", "tophat"))
command.start <- paste("bash -c '", tophat.binary)

# TopHat usage
Usage: tophat [options]* <index_base> <reads1_1[,...,readsN_1]> [reads1_2,...readsN_2] 

# Arguments:
<ebwt_base> 	The basename of the index to be searched. The basename is the name of any of the five index files up to but not including the first period. bowtie first looks in the current directory for the index files, then looks in the indexes subdirectory under the directory where the currently-running bowtie executable is located, then looks in the directory specified in the BOWTIE_INDEXES environment variable.
<reads1_1[,...,readsN_1]> 	A comma-separated list of files containing reads in FASTQ or FASTA format. When running TopHat with paired-end reads, this should be the *_1 ("left") set of files.
<[reads1_2,...readsN_2]> 	A comma-separated list of files containing reads in FASTA or FASTA format. Only used when running TopHat with paired end reads, and contains the *_2 ("right") set of files. The *_2 files MUST appear in the same order as the *_1 files.

# Options: 	
-h/--help 	Prints the help message and exits
-v/--version 	Prints the TopHat version number and exits
-o/--output-dir <string> 	Sets the name of the directory in which TopHat will write all of its output. The default is "./tophat_out".
-r/--mate-inner-dist <int> 	This is the expected (mean) inner distance between mate pairs. For, example, for paired end runs with fragments selected at 300bp, where each end is 50bp, you should set -r to be 200. There is no default, and this parameter is required for paired end runs.
--mate-std-dev <int> 	The standard deviation for the distribution on inner distances between mate pairs. The default is 20bp.
-a/--min-anchor-length <int> 	The "anchor length". TopHat will report junctions spanned by reads with at least this many bases on each side of the junction. Note that individual spliced alignments may span a junction with fewer than this many bases on one side. However, every junction involved in spliced alignments is supported by at least one read with this many bases on each side. This must be at least 3 and the default is 8.
-m/--splice-mismatches <int> 	The maximum number of mismatches that may appear in the "anchor" region of a spliced alignment. The default is 0.
-i/--min-intron-length <int> 	The minimum intron length. TopHat will ignore donor/acceptor pairs closer than this many bases apart. The default is 70.
-I/--max-intron-length <int> 	The maximum intron length. When searching for junctions ab initio, TopHat will ignore donor/acceptor pairs farther than this many bases apart, except when such a pair is supported by a split segment alignment of a long read. The default is 500000.
--allow-indels 	Allow indel search. Indel search is disabled by default.
--max-insertion-length <int> 	The maximum insertion length. The default is 3.
--max-deletion-length <int> 	The maximum deletion length. The default is 3.
--solexa-quals 	Use the Solexa scale for quality values in FASTQ files.
--solexa1.3-quals 	As of the Illumina GA pipeline version 1.3, quality scores are encoded in Phred-scaled base-64. Use this option for FASTQ files from pipeline 1.3 or later.
-Q/--quals 	Separate quality value files - colorspace read files (CSFASTA) come with separate qual files.
--integer-quals 	Quality values are space-delimited integer values, this becomes default when you specify -C/--color.
-C/--color 	Colorspace reads, note that it uses a colorspace bowtie index and requires Bowtie 0.12.6 or higher.
Common usage: tophat --color --quals [other options]* <colorspace_index_base> <reads1_1[,...,readsN_1]> [reads1_2,...readsN_2] <quals1_1[,...,qualsN_1]> [quals1_2,...qualsN_2]
-F/--min-isoform-fraction <0.0-1.0> 	TopHat filters out junctions supported by too few alignments. Suppose a junction spanning two exons, is supported by S reads. Let the average depth of coverage of exon A be D, and assume that it is higher than B. If S / D is less than the minimum isoform fraction, the junction is not reported. A value of zero disables the filter. The default is 0.15.
-p/--num-threads <int> 	Use this many threads to align reads. The default is 1.
-g/--max-multihits <int> 	Instructs TopHat to allow up to this many alignments to the reference for a given read, and suppresses all alignments for reads with more than this many alignments. The default is 20 for read mapping (and it uses two time larger number (40) for segment mapping).
--no-closure-search 	Disables the mate pair closure-based search for junctions. Currently, has no effect - closure search is off by default.
--closure-search 	Enables the mate pair closure-based search for junctions. Closure-based search should only be used when the expected inner distance between mates is small (<= 50bp)
--no-coverage-search 	Disables the coverage based search for junctions.
--coverage-search 	Enables the coverage based search for junctions. Use when coverage search is disabled by default (such as for reads 75bp or longer), for maximum sensitivity.
--microexon-search 	With this option, the pipeline will attempt to find alignments incident to microexons. Works only for reads 50bp or longer.
--butterfly-search 	TopHat will use a slower but potentially more sensitive algorithm to find junctions in addition to its standard search. Consider using this if you expect that your experiment produced a lot of reads from pre-mRNA, that fall within the introns of your transcripts.
--library-type 	TopHat will treat the reads as strand specific. Every read alignment will have an XS attribute tag. Consider supplying library type options below to select the correct RNA-seq protocol.
Library Type	Examples	Description
fr-unstranded	Standard Illumina	Reads from the left-most end of the fragment (in transcript coordinates) map to the transcript strand, and the right-most end maps to the opposite strand.
fr-firststrand	dUTP, NSR, NNSR	Same as above except we enforce the rule that the right-most end of the fragment (in transcript coordinates) is the first sequenced (or only sequenced for single-end reads). Equivalently, it is assumed that only the strand generated during first strand synthesis is sequenced.
fr-secondstrand	Ligation, Standard SOLiD	Same as above except we enforce the rule that the left-most end of the fragment (in transcript coordinates) is the first sequenced (or only sequenced for single-end reads). Equivalently, it is assumed that only the strand generated during second strand synthesis is sequenced.

# Advanced Options: 	
--segment-mismatches 	Read segments are mapped independently, allowing up to this many mismatches in each segment alignment. The default is 2.
--bowtie-n 	TopHat uses "-v" in Bowtie for mapping (the default), but with this option, "-n" is used instead.
--segment-length 	Each read is cut up into segments, each at least this long. These segments are mapped independently. The default is 25.
--min-closure-exon 	During closure search for paired end reads, exonic hops in the potential splice graph must be at least this long. The default is 50.
--min-closure-intron 	The minimum intron length that may be found during closure search. The default is 50.
--max-closure-intron 	The maximum intron length that may be found during closure search. The default is 5000.
--min-coverage-intron 	The minimum intron length that may be found during coverage search. The default is 50.
--max-coverage-intron 	The maximum intron length that may be found during coverage search. The default is 20000.
--min-segment-intron 	The minimum intron length that may be found during split-segment search. The default is 50.
--max-segment-intron 	The maximum intron length that may be found during split-segment search. The default is 500000.
--keep-tmp 	Causes TopHat to preserve its intermediate files produced during the run. By default, they are deleted upon exit.
--no-sort-bam 	Output BAM is not coordinate-sorted.
--no-convert-bam 	Do not convert to bam format. Output is <output_dir>/accepted_hit.sam. Implies --no-sort-bam.
-z/--zpacker 	manually specify the program used for compression of temporary files; default is gzip; use -z0 to disable compression altogether. Any program that is option-compatible with gzip can be used (e.g. bzip2, pigz, pbzip2). 





# Setting up parameters
common parameters
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
bowtie.command <- paste(command.start, common.parameters, mode.parameters, output.parameters, command.end)
#stop(paste('CHIPSTER-NOTE: ', bowtie.command))
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
system("mv alignment.sorted.bam bowtie.bam")
system("mv alignment.sorted.bam.bai bowtie.bam.bai")
