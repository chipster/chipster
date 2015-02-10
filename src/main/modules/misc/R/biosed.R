# TOOL biosed.R: "Replace or delete sequence sections" (Replace or delete sequence sections)
# INPUT sequence: sequence TYPE GENERIC 
# OUTPUT OPTIONAL outseq.fasta
# OUTPUT OPTIONAL biosed.log
# PARAMETER targetregion: "Sequence section to match" TYPE STRING DEFAULT N (Sequence section to match)
# PARAMETER OPTIONAL delete: "Delete the target sequence sections" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Delete the target sequence sections)
# PARAMETER replace: "Replacement sequence section" TYPE STRING DEFAULT A (Replacement sequence section)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)

source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("sequence")

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
	stop(paste('CHIPSTER-NOTE: Too many query sequences. Maximun is 100000 but your file contains ', num.queryseq ))
}

emboss.binary <- file.path(emboss.path, "biosed")
emboss.parameters <- paste("-auto")
emboss.parameters <- paste(emboss.parameters, "-sequence sequence")
emboss.parameters <- paste(emboss.parameters, "-outseq outseq.fasta")
emboss.parameters <- paste(emboss.parameters, "-targetregion", targetregion)
emboss.parameters <- paste(emboss.parameters, "-delete", delete)
emboss.parameters <- paste(emboss.parameters, "-replace", replace)

command.full <- paste(emboss.binary, emboss.parameters, ' >> biosed.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> biosed.log' )
system(echo.command)

system(command.full)

if ( save_log == "no") {
	system ("rm -f biosed.log")
}

