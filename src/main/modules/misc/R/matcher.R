# TOOL matcher.R: "Local pairwise sequence alignment" (Smith-Waterman type local pairwise sequence alignment of two sequences with matcher EMBOSS command)
# INPUT asequence.fa TYPE GENERIC 
# INPUT bsequence.fa TYPE GENERIC 
# OUTPUT OPTIONAL alignment.pair.html
# OUTPUT OPTIONAL alignment.pair.txt
# OUTPUT OPTIONAL alignment.fasta
# OUTPUT OPTIONAL alignment.clustal.txt
# OUTPUT OPTIONAL alignment.phylip.txt
# OUTPUT OPTIONAL matcher.log
# PARAMETER OPTIONAL alternatives: "How many alignments to show" TYPE INTEGER FROM 1 DEFAULT 1 (This sets the number of alternative matches output. By default only the highest scoring alignment is shown. A value of 2 gives you other reasonable alignments. In some cases, for example multidomain proteins of cDNA and genomic DNA comparisons, there may be other interesting and significant alignments.)
# PARAMETER OPTIONAL gapopen: "Gap penalty" TYPE INTEGER FROM 0 DEFAULT 12 (Gap penalty)
# PARAMETER OPTIONAL gapextend: "Gap length penalty" TYPE INTEGER FROM 0 DEFAULT 2 (Gap length penalty)
# PARAMETER OPTIONAL datafile: "Scoring matrix" TYPE [ def: "Default (BLOSUM62 for protein, EDNAFULL for DNA\)", EBLOSUM30: BLOSUM30, EBLOSUM50: BLOSUM50, EBLOSUM80: BLOSUM80] DEFAULT def (Choose the scoring matrix file used when comparing sequences. By default BLOSUM62 (for proteins\) or the file EDNAFULL (for nucleic sequences\) are used)
# PARAMETER OPTIONAL aformat: "Output type for outfile" TYPE [pair: "Pairvise", fasta: FASTA, clustal: ClustalW, phylip: Phylip] DEFAULT pair (Choose format for output file)
# PARAMETER OPTIONAL awidth: "Row length in the alignmnet file" TYPE INTEGER FROM 3 DEFAULT 100 (This parameter defines the row length used in the alinment file)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)

options(scipen=999)
# K.M 28.10.2013
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("asequence.fa")
unzipIfGZipFile("bsequence.fa")

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

#check sequece file type
inputfile.to.check <- ("asequence.fa")
sfcheck.binary <- file.path(chipster.module.path ,"/shell/sfcheck.sh")
sfcheck.command <- paste(sfcheck.binary, emboss.path, inputfile.to.check )
str.filetype <- system(sfcheck.command, intern = TRUE )

if ( str.filetype == "Not an EMBOSS compatible sequence file"){
	stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
}

#count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount asequence.fa -filter")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)

#round(num.queryseq)

if (num.queryseq > 1){
	stop(paste("CHIPSTER-NOTE:Too many query sequences. Maximun is 1 but your file contains ", num.queryseq ))
}

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

#check sequece file type
inputfile.to.check <- ("bsequence.fa")
sfcheck.binary <- file.path(chipster.module.path ,"/shell/sfcheck.sh")
sfcheck.command <- paste(sfcheck.binary, emboss.path, inputfile.to.check )
str.filetype <- system(sfcheck.command, intern = TRUE )

if ( str.filetype == "Not an EMBOSS compatible sequence file"){
	stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
}

#count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount bsequence.fa -filter")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)

#round(num.queryseq)

if (num.queryseq > 1){
	stop(paste("CHIPSTER-NOTE:Too many query sequences. Maximun is 1 but your file contains ", num.queryseq ))
}





outfile <- paste("alignment",aformat ,"txt", sep="." )

if ( aformat == "fasta"){
	outfile <- paste("alignment",aformat, sep="." )	
}

#emboss settings
ecommand <- ("matcher")
emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")
emboss.binary <- file.path(emboss.path,ecommand)
emboss.options <- paste("-asequence asequence.fa -bsequence bsequence.fa -gapopen", gapopen, "-gapextend",  gapextend,"-awidth", awidth, "-aformat", aformat)
emboss.options <- paste(emboss.options, "-alternatives", alternatives )
emboss.options <- paste(emboss.options, "-outfile", outfile, "-auto")
if ( datafile != "def"){
	emboss.options <- paste(emboss.options, "-datafile", datafile )	
}
	

command.full <- paste(emboss.binary, emboss.options, ' >> matcher.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> matcher.log' )
system(echo.command)

system(command.full)
system ("ls -l >>  matcher.log")

if ( aformat == "pair"){
	system("echo '<html>' > alignment.pair.html")
	system("echo '<body>' >> alignment.pair.html")
	system("echo '<pre>' >> alignment.pair.html")
	system("cat alignment.pair.txt >> alignment.pair.html ")
	system("echo '</pre>' >> alignment.pair.html")
	system("echo '</body>' >> alignment.pair.html")
	system("echo '</html>' >> alignment.pair.html")
	system("rm -f alignment.pair.txt")
}



if ( save_log == "no") {
	system ("rm -f matcher.log")
}
