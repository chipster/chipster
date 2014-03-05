# TOOL plotcon.R: "Alignment conservation plot" (Plot conservation of a sequence alignment)
# INPUT OPTIONAL sequences: sequences TYPE GENERIC 
# OUTPUT plotcon{...}.png
# OUTPUT OPTIONAL plotcon.log
# PARAMETER winsize: "Window size" TYPE INTEGER DEFAULT 4 (Number of columns to average alignment quality over. The larger this value is, the smoother the plot will be.)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)

# KM 8.11. 2013
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("sequences")

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

#check sequece file type
inputfile.to.check <- ("sequences")
sfcheck.binary <- file.path(chipster.module.path ,"/shell/sfcheck.sh")
sfcheck.command <- paste(sfcheck.binary, emboss.path, inputfile.to.check )
str.filetype <- system(sfcheck.command, intern = TRUE )

if ( str.filetype == "Not an EMBOSS compatible sequence file"){
	stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
}

#count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount sequences -filter")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)

#round(num.queryseq)

if (num.queryseq > 10000){
	stop(paste('CHIPSTER-NOTE: Too many query sequences. Maximun is 10000 but your file contains ', num.queryseq ))
}

emboss.binary <- file.path(emboss.path, "plotcon")
emboss.parameters <- paste("-auto")
emboss.parameters <- paste(emboss.parameters, "-sequences sequences")
emboss.parameters <- paste(emboss.parameters, "-graph png")
emboss.parameters <- paste(emboss.parameters, "-winsize", winsize)

command.full <- paste(emboss.binary, emboss.parameters, ' >> plotcon.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> plotcon.log' )
system(echo.command)

system(command.full)

if ( save_log == "no") {
	system ("rm -f plotcon.log")
}
