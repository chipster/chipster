# TOOL dotmatcher.R: "Threshold dotplot of two sequences" (Dotmatcher: Draw a threshold dotplot of two sequences.)
# INPUT asequence: asequence TYPE GENERIC 
# INPUT bsequence: bsequence TYPE GENERIC 
# OUTPUT OPTIONAL dotmatcher.png 
# OUTPUT OPTIONAL dotmatcher.log
# PARAMETER OPTIONAL windowsize: "Window size over which to test threshold" TYPE INTEGER FROM 3 DEFAULT 10 (Window size over which to test threshold)
# PARAMETER OPTIONAL threshold: Threshold TYPE INTEGER FROM 0 DEFAULT 23 (Threshold)
# PARAMETER OPTIONAL stretch: "Stretch plot" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Display a non-proportional graph)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

#check sequece file type
inputfile.to.check <- ("asequence")
sfcheck.binary <- file.path(chipster.module.path ,"/shell/sfcheck.sh")
sfcheck.command <- paste(sfcheck.binary, emboss.path, inputfile.to.check )
str.filetype <- system(sfcheck.command, intern = TRUE )

if ( str.filetype == "Not an EMBOSS compatible sequence file"){
	stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
}

#count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount asequence -filter")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)

#round(num.queryseq)

if (num.queryseq > 1){
	stop(paste("CHIPSTER-NOTE:Too many query sequences. Maximun is 1 but your file contains ", num.queryseq ))
}

#check sequece file type
inputfile.to.check <- ("bsequence")
sfcheck.binary <- file.path(chipster.module.path ,"/shell/sfcheck.sh")
sfcheck.command <- paste(sfcheck.binary, emboss.path, inputfile.to.check )
str.filetype <- system(sfcheck.command, intern = TRUE )

if ( str.filetype == "Not an EMBOSS compatible sequence file"){
	stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
}

#count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount bsequence -filter")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)

#round(num.queryseq)

if (num.queryseq > 1){
	stop(paste("CHIPSTER-NOTE:Too many query sequences. Maximun is 1 but your file contains ", num.queryseq ))
}

emboss.binary <- file.path(emboss.path, "dotmatcher")
emboss.parameters <- paste("-asequence asequence")
emboss.parameters <- paste(emboss.parameters, "-bsequence bsequence")
emboss.parameters <- paste(emboss.parameters, "-graph png")
emboss.parameters <- paste(emboss.parameters, "-xygraph png")
emboss.parameters <- paste(emboss.parameters, "-windowsize", windowsize)
emboss.parameters <- paste(emboss.parameters, "-threshold", threshold)
emboss.parameters <- paste(emboss.parameters, "-stretch", stretch)

command.full <- paste(emboss.binary, emboss.parameters, ' >> dotmatcher.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> dotmatcher.log' )
system(echo.command)

system(command.full)
system("mv dotmatcher.1.png dotmatcher.png")
system("ls -l >> dotmatcher.log")
if ( save_log == "no") {
	system ("rm -f dotmatcher.log")
}
