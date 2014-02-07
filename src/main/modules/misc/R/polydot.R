# TOOL polydot.R: "Dotplots for a sequence set" (Plydor. Draw dotplots for all-against-all comparison of a sequence set)
# INPUT OPTIONAL sequences: sequences TYPE GENERIC 
# OUTPUT OPTIONAL polydot.log
# OUTPUT OPTIONAL polydot.png
# PARAMETER wordsize: "Word size" TYPE INTEGER FROM 2 DEFAULT 6 (Word size)
# PARAMETER OPTIONAL gap: "Gap between dotplots" TYPE INTEGER FROM 0 DEFAULT 10 (This specifies the size of the gap that is used to separate the individual dotplots in the display. The size is measured in residues, as displayed in the output.)
# PARAMETER OPTIONAL boxit: "Draw a box around each dotplot" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Draw a box around each dotplot)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)

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

if (num.queryseq > 100){
	stop(paste('CHIPSTER-NOTE: Too many query sequences. Maximun is 100 but your file contains ', num.queryseq ))
}

emboss.binary <- file.path(emboss.path, "polydot")
emboss.parameters <- paste("-auto")
emboss.parameters <- paste(emboss.parameters, "-sequences sequences")
emboss.parameters <- paste(emboss.parameters, "-graph png")
emboss.parameters <- paste(emboss.parameters, "-wordsize", wordsize)
emboss.parameters <- paste(emboss.parameters, "-gap", gap)
emboss.parameters <- paste(emboss.parameters, "-boxit", boxit)


command.full <- paste(emboss.binary, emboss.parameters, ' >> polydot.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> polydot.log' )
system(echo.command)

system(command.full)
system("mv polydot.1.png polydot.png")

if ( save_log == "no") {
	system ("rm -f polydot.log")
}
