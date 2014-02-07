# TOOL extractseq.R: "Extract regions from a sequence" (Extract regions from a sequence)
# INPUT sequence: sequence TYPE GENERIC 
# OUTPUT OPTIONAL extractseq.fasta: extractseq.fasta 
# OUTPUT OPTIONAL extractseq.log
# PARAMETER regions: "Regions to extract" TYPE STRING (Regions to extract.  A set of regions is specified by a set of pairs of positions. The positions are integers.  They are separated by any non-digit, non-alpha character. Examples of region specifications are: 24-45, 56-78  1:45, 67=99;765..888 1,5,8,10,23,45,57,99)
# PARAMETER OPTIONAL separate: "Write regions to separate sequences" TYPE [ Y: Yes, N: No] DEFAULT N (If this is set true then each specified region is written out as a separate sequence. The name of the sequence is created from the name of the original sequence with the start and end positions of the range appended with underscore characters between them, eg: XYZ region 2 to 34 is written as: XYZ_2_34)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

#check sequece file type
inputfile.to.check <- ("sequence")
sfcheck.binary <- file.path(chipster.module.path ,"/shell/sfcheck.sh")
sfcheck.command <- paste(sfcheck.binary, emboss.path, inputfile.to.check )
str.filetype <- system(sfcheck.command, intern = TRUE )

if ( str.filetype == "Not an EMBOSS compatible sequence file"){
	stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
}

#count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount sequence -filter")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)

#round(num.queryseq)

if (num.queryseq > 100000){
	stop(paste("CHIPSTER-NOTE: Too many query sequences. Maximun is 100000 but your file contains ", num.queryseq ))
}

emboss.binary <- file.path(emboss.path, "extractseq")
emboss.parameters <- paste("-sequence sequence")
emboss.parameters <- paste(emboss.parameters, "-outseq extractseq.fasta")
emboss.parameters <- paste(emboss.parameters, "-regions", regions)
emboss.parameters <- paste(emboss.parameters, "-separate", separate)


command.full <- paste(emboss.binary, emboss.parameters, ' >> extractseq.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> extractseq.log' )
system(echo.command)

system(command.full)
system ("ls -l >> extractseq.log")

if ( save_log == "no") {
	system ("rm -f extractseq.log")
}
