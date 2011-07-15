# TOOL bedtools-test.R: "Bedtools test tool" (Heippa hei!)
# INPUT reads.fastq TYPE GENERIC
# OUTPUT result.txt 
# PARAMETER first: "First base to keep" TYPE INTEGER FROM 1 TO 100 DEFAULT 1 (First base to keep.)
# PARAMETER last: "Last base to keep" TYPE INTEGER FROM 1 TO 100 DEFAULT 50 (Last base to keep.)




# EK 17.6.2011

# binary
binary <- "ls > result.txt"

# run
system(binary)