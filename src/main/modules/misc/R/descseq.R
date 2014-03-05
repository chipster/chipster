# TOOL descseq.R: "Modify sequence name or description" (Alter the name or description of a sequence with EMBOSS tool descseq.)
# INPUT sequence: "Input sequence" TYPE GENERIC 
# OUTPUT OPTIONAL renamed.fasta
# OUTPUT OPTIONAL descseq.log 
# PARAMETER OPTIONAL name: "Name of the sequence" TYPE STRING (Name of the sequence)
# PARAMETER OPTIONAL description: "Description of the sequence" TYPE STRING (Description of the sequence)
# PARAMETER OPTIONAL append: "Append to the existing description" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (This allows you to append the name or description you have given on to the end of the existing name or description of the sequence.)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("sequence")

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

if (num.queryseq > 1){
	stop(paste("CHIPSTER-NOTE:Too many query sequences. Maximun is 1 but your file contains ", num.queryseq ))
}

emboss.binary <- file.path(emboss.path, "descseq")
emboss.parameters <- paste("-auto")
emboss.parameters <- paste(emboss.parameters, "-sequence sequence")
emboss.parameters <- paste(emboss.parameters, "-outseq renamed.fasta")
if (nchar(name) > 0 ) {
emboss.parameters <- paste(emboss.parameters, "-name", name)
}
if (nchar(description) > 0 ) {
emboss.parameters <- paste(emboss.parameters, "-description", description)
}
emboss.parameters <- paste(emboss.parameters, "-append", append)

command.full <- paste(emboss.binary, emboss.parameters, ' >> descseq.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> descseq.log' )
system(echo.command)

system(command.full)

if ( save_log == "no") {
	system ("rm -f descseq.log")
}
