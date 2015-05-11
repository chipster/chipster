# TOOL macs2.R: "Find peaks using MACS2" (Detects statistically significantly enriched genomic regions in ChIP-seq data, using a control sample if available. If you have several samples,you need to merge them first to one ChIP file and one control file. BAM files can be merged with the Utilities tool \"Merge BAM\".)
# INPUT treatment.bam: "Treatment data file" TYPE GENERIC 
# INPUT OPTIONAL control.bam: "Control data file" TYPE GENERIC 
# OUTPUT macs2-log.txt
# OUTPUT OPTIONAL macs2-peaks.tsv
# OUTPUT OPTIONAL macs2-summits.bed
# OUTPUT OPTIONAL macs2_narrowpeak.bed
# OUTPUT OPTIONAL macs2_broad_peaks.bed
# OUTPUT OPTIONAL macs2_model.pdf
# PARAMETER file.format: "Input file format" TYPE [ELAND, BAM, BED] DEFAULT BAM (The format of the input files.)
# PARAMETER precalculated.size: "Mappable genome size" TYPE [2.7e9: "human hg18 (2.7e9\)", 2.72e9: "human hg19 (2.72e9\)", 1.87e9: "mouse mm9 (1.87e9\)", 1.89e9: "mouse mm10 (1.89e9\)", 2.32e9: "rat rn5 (2.32e9\)", user_specified: "User specified"] DEFAULT 2.72e9 (Mappable genome size. You can use one of the precalculated ones or choose User specified and provide the size in the field below.)
# PARAMETER OPTIONAL userspecified.size: "User specified mappable genome size" TYPE STRING (You can also use scientific notation, e.g. 1.23e9 . Remember to select User specified as Mappable genome size.)
# PARAMETER OPTIONAL q.value.threshold: "q-value cutoff" TYPE DECIMAL FROM 0 TO 0.99 DEFAULT 0.01 (The minimum FDR for peak detection.)
# PARAMETER OPTIONAL read.length: "Read length" TYPE INTEGER FROM 0 TO 200 DEFAULT 0 (The read length in nucleotides. Read length is autodetected if you set this to 0.)
# PARAMETER OPTIONAL keep.dup: "Keep duplicate reads" TYPE [auto, all, 1] DEFAULT auto (Procedure used to handle duplicate reads. If auto, MACS computes the maximum reads at the exact same location based on binomal distribution using 1e-5 p-value cutoff. The option All option keeps all reads, while 1 keeps only one read per site.)
# PARAMETER OPTIONAL build.model: "Build peak model" TYPE [yes, no] DEFAULT yes (If enabled, a peak model is built from the data. Disabling model building means the extension size has to be guessed and set with the parameter.)
# PARAMETER OPTIONAL bandwidth: "Bandwidth" TYPE INTEGER FROM 1 TO 1000 DEFAULT 300 (The window size for picking regions to compute fragment size when building the shifting model.)
# PARAMETER OPTIONAL ext.size: "Extension size" TYPE INTEGER FROM 1 TO 1000 DEFAULT 200 (When model building has been switched off or when it fails, MACS will extend the reads to this length.)
# PARAMETER OPTIONAL m.fold.upper: "Upper M-fold cutoff" TYPE INTEGER FROM 1 TO 100 DEFAULT 30 (Sets the cutoff used to determine peak regions for model building. A too high value may result in not enough peaks being identified for building the model. Note that if the peak model is disabled this parameter has no effect.)
# PARAMETER OPTIONAL m.fold.lower: "Lower M-fold cutoff" TYPE INTEGER FROM 1 TO 100 DEFAULT 10 (Sets the cutoff used to determine peak regions for model building. A too low value may result in the inclusion of many false peaks being used for building the model. Note that if the peak model is disabled this parameter has no effect.)
# PARAMETER OPTIONAL broad: "Call broad peaks" TYPE [yes, no] DEFAULT no (Call broad peaks by linking nearby highly enriched region.)

# 26.05.2010 MG, Created
# 01.12.2012 MG, Modified to take BAM files as input. Modified to run version 1.4 of MACS.
# 08.03.2011 MG, Modified to disable wiggle output.
# 05.04.2014 MK, Polished. Added MACS2
# 10.07.2014 AMS, Updated genome sizes, added parameter userspecified.size
# 09.09.2014 EK, Made a separate script for MACS2 in order to cope with new parameters, fixed the bug in disabled model building, added the broad option and outputs, polished the script and output      
# 07.10.2014 AMS, Simplified script structure
# 13.10.2014 EK, Modified to use MACS2.1.0

# MACS binary
macs.binary <- file.path(chipster.tools.path, "macs", "macs2")

# Options

# treatment data file
options <- paste ("callpeak -t treatment.bam")

# control data file (optional)
if (file.exists("control.bam")){
	options <- paste(options, "-c control.bam")
}

# output file prefix
options <- paste(options, "-n macs2")

# input file format
options <- paste(options, "-f", file.format)

# mappable genome size
if (precalculated.size == "user_specified") {
	if (nchar(userspecified.size) < 1){
		stop(paste('CHIPSTER-NOTE: ', "You need to provide a value for mappable genome size or select one of the precalculated values."))
	}
	genome.size <- userspecified.size
}else{
	genome.size <- precalculated.size
}
options <- paste(options, "-g", genome.size)

# q-value cutoff
options <- paste(options, "-q", q.value.threshold)

# read length
if (read.length > 0) {
	options <- paste(options, "-s", read.length)
}

# keep duplicate reads
options <- paste(options, "--keep-dup", keep.dup)

# build peak model
if (build.model == "no") {
	options <- paste(options, "--nomodel")
} else{
	options <- paste(options, "--fix-bimodal")
}

# bandwidth (only applicable to model building)
if (build.model == "yes") {
	options <- paste(options, "--bw", bandwidth)
} 

# extension size
options <- paste(options, "--extsize", ext.size)

# Set up the m-fold limits
options <- paste(options, "-m", m.fold.lower, m.fold.upper)

# call broad peaks
if (broad == "yes") {
	options <- paste(options, "--broad")
} 

# common options
options <- paste(options, "--verbose=2")

# Run macs
macs.command <- paste(macs.binary, options, "2> macs2-log.txt")
# stop(paste('CHIPSTER-NOTE: ', macs.command))
system(macs.command)

# Read in and parse the results (rename and the p- and q-value columns, sort)
output <- try(read.table(file="macs2_peaks.xls", skip=0, header=TRUE, stringsAsFactors=FALSE))
if (class(output) != "try-error") {
	colnames(output)[7] <- "neglog10pvalue"
	colnames(output)[9] <- "neglog10qvalue"
	# output <- output[ order(output[,9], decreasing=TRUE), ]
	output <- output[order(output$chr, output$start),]
	write.table(output, file="macs2-peaks.tsv", sep="\t", quote=FALSE, row.names=FALSE)
}

# Sort the summit BED
source(file.path(chipster.common.path, "bed-utils.R"))
bed <- try(read.table(file="macs2_summits.bed", skip=0, sep="\t"))
if (class(bed) != "try-error") {
	colnames(bed)[1:2] <- c("chr", "start")
	bed <- sort.bed(bed)
	write.table(bed, file="macs2-summits.bed", sep="\t", row.names=F, col.names=F, quote=F)
}

# Add BED extension to the narrow peak format file
if (file.exists("macs2_peaks.narrowPeak") && file.info("macs2_peaks.narrowPeak")$size > 0){
	system ("mv macs2_peaks.narrowPeak macs2_narrowpeak.bed")
}

# Source the R code for plotting the MACS model
if (build.model == "yes") {
	try(source("macs2_model.r"), silent=TRUE)
}

