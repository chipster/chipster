# TOOL plotorf.R: "Plot open reading frames" (Plot potential open reading frames in a nucleotide sequence)
# INPUT sequence: "Input sequence" TYPE GENERIC 
# OUTPUT OPTIONAL plotorf.png
# OUTPUT OPTIONAL plotorf.svg 
# OUTPUT OPTIONAL plotorf.pdf 
# OUTPUT OPTIONAL plotorf.log 
# PARAMETER OPTIONAL start: "Start codons" TYPE STRING DEFAULT ATG (Start codons)
# PARAMETER OPTIONAL stop: "Stop codons" TYPE STRING DEFAULT "TAA,TAG,TGA" (Stop codons)
# PARAMETER OPTIONAL graph: "Output image file format" TYPE [png: "PNG", svg: "SVG", cps: "PDF"] DEFAULT png (Output image file format)
# PARAMETER OPTIONAL save_log: "Output a log file" TYPE [yes: yes, no: no] DEFAULT no (Collect a log file.)

# KM 8.11. 2013
options(scipen=999)
emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

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

if (num.queryseq > 1) {
	stop(paste('CHIPSTER-NOTE: Too many query sequences. Maximun is 1 but your file contains ', num.queryseq ))
}


emboss.binary <- file.path(emboss.path, "plotorf")
emboss.parameters <- paste('sequence -auto -goutfile plotorf -graph', graph)
emboss.parameters <- paste(emboss.parameters, '-start "', start, '"')
emboss.parameters <- paste(emboss.parameters, '-stop "', stop, '"')




command.full <- paste(emboss.binary, emboss.parameters, ' >> plotorf.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> plotorf.log' )
system(echo.command)

system(command.full)
if ( graph == "png") {
	system ("mv plotorf.1.png plotorf.png")
}

if ( graph == "cps") {
	system ("ps2pdf plotorf.ps >> plotorf.log ") 
}


if ( save_log == "no") {
	system ("rm -f plotorf.log")
}
