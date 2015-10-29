# TOOL picard-markduplicates.R: "Mark duplicates in BAM" (Examines aligned records in the supplied BAM file to locate duplicate molecules. All records are then written to the output file with the duplicate records flagged. This tool is based on the Picard Tools package.)
# INPUT alignment.bam: "BAM alignment file" TYPE GENERIC
# OUTPUT OPTIONAL marked.bam
# OUTPUT OPTIONAL marked.bam.bai
# OUTPUT OPTIONAL duplicateMetrics.tsv 

# 2015.09.09 AMS 

picard.binary <- file.path(chipster.tools.path, "picard-tools", "picard.jar")
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# run
picard.command <- paste("java -Xmx2g -jar", picard.binary, "MarkDuplicates INPUT=alignment.bam OUTPUT=marked.bam METRICS_FILE=duplicateMetrics.txt VALIDATION_STRINGENCY=LENIENT")
system(picard.command)

system(paste(samtools.binary, "index marked.bam > marked.bam.bai"))

# duplicateMetrics 
system("grep -A2 LIBRARY duplicateMetrics.txt > duplicateMetrics.tsv")

# Handle output names
source(file.path(chipster.common.path, "tool-utils.R"))

# read input names
inputnames <- read_input_definitions()

basename  <- strip_name(inputnames$alignment.bam)

# Make a matrix of output names
outputnames <- matrix(NA, nrow=2, ncol=2)
outputnames[1,] <- c("marked.bam", paste(basename, "_dedup.bam", sep =""))
outputnames[2,] <- c("marked.bam.bai", paste(basename, "_dedup.bam.bai", sep =""))

# Write output definitions file
write_output_definitions(outputnames)
