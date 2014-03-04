# TOOL skipredundant.R: "Skip redundant sequences" (Remove redundant sequences from an input set)
# INPUT sequence: "Input sequence" TYPE GENERIC 
# OUTPUT OPTIONAL non_redundant.fasta
# OUTPUT OPTIONAL removed_redundant.fasta
# OUTPUT OPTIONAL skipredundant.log 
# PARAMETER mode: "Select number" TYPE [1: "Single threshold percentage sequence similarity", 2: "Outside a range of acceptable threshold percentage similarities"] FROM 1 TO 1 DEFAULT 1 (This option specifies whether to remove redundancy at a single threshold percentage sequence similarity or remove redundancy outside a range of acceptable threshold percentage similarity. All permutations of pair-wise sequence alignments are calculated for each set of input sequences in turn using the EMBOSS implementation of the Needleman and Wunsch global alignment algorithm. Redundant sequences are removed in one of two modes as follows: (i\) If a pair of proteins achieve greater than a threshold percentage sequence similarity (specified by the user\) the shortest sequence is discarded. (ii\) If a pair of proteins have a percentage sequence similarity that lies outside an acceptable range (specified by the user\) the shortest sequence is discarded.)
# PARAMETER threshold: "The percentage sequence identity redundancy threshold." TYPE DECIMAL DEFAULT 95.0 (This option specifies the percentage sequence identity redundancy threshold. The percentage sequence identity redundancy threshold  determines the redundancy calculation. If a pair of proteins achieve greater than this threshold the shortest sequence is discarded.)
# PARAMETER OPTIONAL minthreshold: "The % sequence identity redundancy threshold (lower limit\)." TYPE DECIMAL DEFAULT 30.0 (This option specifies the percentage sequence identity redundancy threshold (lower limit\). The percentage sequence identity redundancy threshold determines the redundancy calculation. If a pair of proteins have a percentage sequence similarity that lies outside an acceptable range the shortest sequence is discarded.)
# PARAMETER OPTIONAL maxthreshold: "The percentage sequence identity redundancy threshold (upper limit\)." TYPE DECIMAL DEFAULT 90.0 (This option specifies the percentage sequence identity redundancy threshold (upper limit\). The percentage sequence identity redundancy threshold determines the redundancy calculation. If a pair of proteins have a percentage sequence similarity that lies outside an acceptable range the shortest sequence is discarded.)
# PARAMETER gapopen: "Gap opening penalty" TYPE DECIMAL FROM 0.0 TO 100.0 DEFAULT 10.0 (The gap open penalty is the score taken away when a gap is created. The best value depends on the choice of comparison matrix. The default value assumes you are using the EBLOSUM62 matrix for protein sequences, and the EDNAFULL matrix for nucleotide sequences.)
# PARAMETER gapextend: "Gap extension penalty" TYPE DECIMAL FROM 0.0 TO 10.0 DEFAULT 0.5 (The gap extension, penalty is added to the standard gap penalty for each base or residue in the gap. This is how long gaps are penalized. Usually you will expect a few long gaps rather than many short gaps, so the gap extension penalty should be lower than the gap penalty. An exception is where one or both sequences are single reads with possible sequencing errors in which case you would expect many single base gaps. You can get this result by setting the gap open penalty to zero (or very low\) and using the gap extension penalty to control gap scoring.)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)


# KM 23.1. 2014
options(scipen=999)
emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("sequence")


#check sequece file type
sfcheck.binary <- file.path(chipster.module.path ,"/shell/sfcheck.sh")
sfcheck.command <- paste(sfcheck.binary, emboss.path, "sequence" )
str.filetype <- system(sfcheck.command, intern = TRUE )

if ( str.filetype == "Not an EMBOSS compatible sequence file"){
	stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
}

#count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount -filter sequence")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)
#round(num.queryseq)

if (num.queryseq > 50000) {
	stop(paste('CHIPSTER-NOTE: Too many query sequences. Maximun is 50000 but your file contains ', num.queryseq ))
}


emboss.binary <- file.path(emboss.path, "skipredundant")
emboss.parameters <- paste('sequence -auto -outseq non_redundant.fasta -redundantoutseq removed_redundant.fasta ' )
emboss.parameters <- paste(emboss.parameters, "-mode", mode)
emboss.parameters <- paste(emboss.parameters, "-threshold", threshold)
emboss.parameters <- paste(emboss.parameters, "-minthreshold", minthreshold)
emboss.parameters <- paste(emboss.parameters, "-maxthreshold", maxthreshold)
emboss.parameters <- paste(emboss.parameters, "-gapopen", gapopen)
emboss.parameters <- paste(emboss.parameters, "-gapextend", gapextend)



command.full <- paste(emboss.binary, emboss.parameters, ' >> skipredundant.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> skipredundant.log' )
system(echo.command)

system(command.full)



if ( save_log == "no") {
	system ("rm -f skipredundant.log")
}
