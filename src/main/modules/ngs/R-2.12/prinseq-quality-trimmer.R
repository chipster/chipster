# TOOL prinseq-quality-trimmer.R: "Trim reads by quality" (Trims reads based on the quality values. This tool is based on the PRINSEQ package.)
# INPUT fastqfile: "Input reads set" TYPE GENERIC
# OUTPUT trimmed.fastq
# OUTPUT OPTIONAL trim.log
# PARAMETER OPTIONAL trim.qual.right: "Trim 3-prime-end by quality" TYPE INTEGER (Trim sequence by quality score from the 3-prime-end with the given threshold score.)		
# PARAMETER OPTIONAL trim.qual.left: "Trim 5-prime-end by quality" TYPE INTEGER (Trim sequence by quality score from the 5-prime-end with the given threshold score.)
# PARAMETER OPTIONAL trim.qual.type: "The quality score calculation method" TYPE [min: "minimum quality value", mean: "mean quality", max: "maximum quality value", sum: "total sum of quality values"] DEFAULT min (Type of quality score calculation to use.)
# PARAMETER OPTIONAL trim.qual.rule: "Quality score comparison condition" TYPE [ lt: "Less than", gt: "greater than", et: "equal to"] DEFAULT lt (Rule to use to compare quality score threshold to calculated value.)
# PARAMETER OPTIONAL trim.qual.window: "Window size for quality calculation" TYPE INTEGER DEFAULT 1 ( The sliding window size used to calculate quality score by type. To stop at the first base that fails the rule defined, use a window size of 1.)
# PARAMETER OPTIONAL trim.qual.step: "Step size used to move the quality window" TYPE INTEGER DEFAULT 1 (Step size used to move the sliding window. To move the window over all quality scores without missing any, the step size should be less or equal to the window size.)
# PARAMETER OPTIONAL phred64: "Quality data is in Phred+64 format" TYPE [ n: "no", y: "yes"] DEFAULT n ( You should select \"yes" option if the quality data in FASTQ file is in Phred+64 format. For Illumina 1.8+, Sanger, Roche/454, Ion Torrent, PacBio data, you should use the default value: no)
# PARAMETER OPTIONAL log.file: "Write a log file" TYPE [ n: "no", y: "yes"] DEFAULT n (Write a log file)



# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("fastqfile")

# binary
binary.prinseq <- c(file.path(chipster.tools.path, "prinseq", "prinseq-lite.pl" ))


trim.params <- paste(" ")

if (!is.na(trim.qual.right)) {
	trim.params <- paste(trim.params, "-trim_qual_right",  trim.qual.right, "-trim_qual_type", trim.qual.type, "-trim_qual_rule", trim.qual.rule, "-trim_qual_window", trim.qual.window, "-trim_qual_step",  trim.qual.step)
}

if (!is.na(trim.qual.left)) {
	trim.params <- paste(trim.params, "-trim_qual_left",  trim.qual.left, "-trim_qual_type", trim.qual.type, "-trim_qual_rule", trim.qual.rule, "-trim_qual_window", trim.qual.window, "-trim_qual_step",  trim.qual.step)
}

if ( phred64 == "y") {
	trim.params <- paste(trim.params, "-phred64")
}


trim.command <- paste(binary.prinseq, trim.params, "-fastq fastqfile -out_good trimmed")


if (log.file == "y") {
	system("echo Running PRINSEQ filtering with command: > trim.log")
	echo.command <- paste("echo '", trim.command, "'>> trim.log")
	system(echo.command)
	trim.command <- paste(trim.command, "-verbose 2>> trim.log")
}


system(trim.command)

#Make sure something is in the output
system("if [ ! -e  trimmed.fastq ] ; then echo 'Trimming produced an empty trimmed.fastq sequence set' >> trim.log ; echo '' > trimmed.fastq ; fi")

#stop

