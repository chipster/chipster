# TOOL fuzzpro.R: "Protein pattern search" (Search for patterns in protein sequences)
# INPUT OPTIONAL sequence: sequence TYPE GENERIC 
# OUTPUT OPTIONAL pattern_matches.txt
# OUTPUT OPTIONAL pattern_matches.tsv
# OUTPUT OPTIONAL fuzznuc.log
# PARAMETER pattern TYPE STRING (pattern)
# PARAMETER OPTIONAL rformat: "Output format type" TYPE [excel: "Table", table: "Text formatted report", gff: "GFF3 formatted file", tagseq: "Tagseq format", listfile: "EMBOSS list file"] DEFAULT excel (Output format type)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")
options(scipen=999)

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

emboss.binary <- file.path(emboss.path, "fuzzpro")
emboss.parameters <- paste("-auto")
emboss.parameters <- paste(emboss.parameters, "-sequence sequence")
emboss.parameters <- paste(emboss.parameters, "-pattern", pattern)
emboss.parameters <- paste(emboss.parameters, "-rformat", rformat)
emboss.parameters <- paste(emboss.parameters, "-outfile pattern_matches.txt")

command.full <- paste(emboss.binary, emboss.parameters, ' >> fuzznuc.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> fuzznuc.log' )
system(echo.command)

system(command.full)


if (rformat == "excel") {
	system ("mv pattern_matches.txt pattern_matches.tsv")	
}

system("ls -l >> fuzznuc.log")

if ( save_log == "no") {
	system ("rm -f fuzznuc.log")
}
