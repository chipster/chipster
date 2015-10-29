# TOOL unique-alignments-from-bam.R: "Retrieve unique alignments from BAM" (Retrieve unique alignments from BAM. Input BAM needs to be sorted. This tool is based on the SAMtools package.)
# INPUT alignment.bam: "BAM alignment file" TYPE GENERIC
# OUTPUT OPTIONAL unique_alignments.bam
# PARAMETER paired: "Does the alignment file contain paired-end data" TYPE [yes, no] DEFAULT no (Does the alignment data contain paired end or single end reads?)
# PARAMETER OPTIONAL nhtag: "Add NH tag" TYPE [yes, no] DEFAULT no (Add NH tag.)

# 2014.06.02 AMS 

# convert bam to sam, sort bam if the data is paired-end
samtools.binary <- file.path(chipster.tools.path, "samtools", "samtools")
picard.binary <- file.path(chipster.tools.path, "picard-tools", "picard.jar")

if(paired == "yes") {
	#List of IDs having mapping quality >= 4 and having proper pairing (0x0002). Uniq command make sure that both pairs must fulfill filtering criteria
	keep_id.command <- paste(samtools.binary, "view -q4 -F4 -f2 alignment.bam | cut -f1 | sort | uniq -d > alignment_id_keep.txt")
	system(keep_id.command)
	filtersamreads.command <- paste("java -Xmx4096m -jar",picard.binary, "FilterSamReads INPUT=alignment.bam OUTPUT=alignment_keep.bam FILTER=includeReadList READ_LIST_FILE=alignment_id_keep.txt  VALIDATION_STRINGENCY=LENIENT")
	system(filtersamreads.command)
	if (nhtag == "yes"){
		samtools.view <- paste(samtools.binary, "view -h alignment_keep.bam | perl -p -e 's/\n$/\tNH:i:1\n/g'")
	}else{
		samtools.view <- paste(samtools.binary, "view -h alignment_keep.bam")
	}
} else if(paired == "no") {
if (nhtag == "yes"){
	samtools.view <- paste(samtools.binary, "view -h -q4 -F4 alignment.bam | perl -p -e 's/\n$/\tNH:i:1\n/g'")
} else{
	samtools.view <- paste(samtools.binary, "view -h -q4 -F4 alignment.bam")
}
}
# run
command <- paste(samtools.view, "|", samtools.binary, "view -bS - > unique_alignments.bam")
system(command)

# Handle output names
source(file.path(chipster.common.path, "tool-utils.R"))

# read input names
inputnames <- read_input_definitions()

# Make a matrix of output names
outputnames <- matrix(NA, nrow=1, ncol=2)
outputnames[1,] <- c("unique_alignments.bam", paste(inputnames$alignment.bam))

# Write output definitions file
write_output_definitions(outputnames)
