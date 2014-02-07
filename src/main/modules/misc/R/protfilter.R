# TOOL protfilter.R: "Filter protein sequence sets" (Select sequences from protein sequence set using different criteria)
# INPUT sequence: sequence TYPE GENERIC
# OUTPUT OPTIONAL filtered_prot.fasta
# OUTPUT OPTIONAL protfilter.log
# PARAMETER OPTIONAL desctext: "Enter a pattern to search sequence descriotion field" TYPE STRING DEFAULT "" (The sequcences, that have matching string in sequence descriotion will be selected.)
# PARAMETER OPTIONAL casesensitive: "Do a case-sensitive search" TYPE [Y: Yes, N: No] DEFAULT N (Do a case-sensitive decriotion text search)
# PARAMETER OPTIONAL pattern: "Search based on sequence patterns" TYPE STRING DEFAULT "" (The sequcences, that match the given sequence pattern are selected. The pattern search is performed with EMBOSS command fuzzpro. Pattern example\:\[FY\]-\[LIV\]-G-\[DE\]-E-A-Q-x-\[RKQ\]\(2\)-G )
# PARAMETER OPTIONAL maxlength: "Maximum length of the protein" TYPE INTEGER DEFAULT 100000 (Sequences that are longer than the given value are filterd out.)
# PARAMETER OPTIONAL minlength: "Minimum length of the protein" TYPE INTEGER DEFAULT 1 (Sequences that are shorter than the given value are filterd out.)
# PARAMETER OPTIONAL maxmass: "Maximum mass of the protein" TYPE INTEGER DEFAULT 100000  (Sequences that are heavier than the given value are filterd out.)
# PARAMETER OPTIONAL minmass: "Minimum mass of the protein" TYPE INTEGER DEFAULT 1 (Sequences that are lighter than the given value are filterd out.)
# PARAMETER OPTIONAL save_log: "Output a log file" TYPE [yes: yes, no: no] DEFAULT no (Collect a log file.)

# KM 23.1. 2014
#settings
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
seqcount.exe <- file.path(emboss.path, "seqcount -filter sequence")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)
#round(num.queryseq)

if (num.queryseq > 50000){
	stop(paste('CHIPSTER-NOTE: Too many query sequences. Maximun is 50000 but your file contains ', num.queryseq ))
}

command.path <- file.path(chipster.module.path ,"/shell/protfilter.sh")
protfilter.command <- paste(command.path, "-emboss_path ", emboss.path, " -sequence sequence -outfile filtered_prot.fasta")
   emboss.parameters <- paste('-casesensitive' ,casesensitive)
if (nchar(desctext) > 0 ) {
	emboss.parameters <- paste(emboss.parameters, '-desctext "', desctext, '"')
}

if (nchar(pattern) > 0 ) {
	emboss.parameters <- paste(emboss.parameters, '-pattern "', pattern,'"')
}

if ( maxlength < 100000 ) {
	emboss.parameters <- paste(emboss.parameters, '-maxlength', maxlength)
}

if (minlength > 1 ) {
	emboss.parameters <- paste(emboss.parameters, '-minlength', minlength)
}

if (maxmass < 100000 ) { 
	emboss.parameters <- paste(emboss.parameters, '-maxmass ', maxmass)
}

if (minmass > 1  ) {
	emboss.parameters <- paste(emboss.parameters, '-minmass ', minmass)
}

command.full <- paste(protfilter.command, emboss.parameters, ' >> protfilter.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> protfilter.log' )
system(echo.command)
system(command.full)

if ( save_log == "no") {
	system ("rm -f protfilter.log")
}