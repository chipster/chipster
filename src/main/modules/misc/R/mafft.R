# TOOL mafft.R: "Multiple sequence alignment with MAFFT" (Multiple sequence alignmnet with MAFFT software.)
# INPUT sequence: "Query sequences" TYPE GENERIC 
# OUTPUT OPTIONAL alignment.aln.txt
# OUTPUT OPTIONAL alignment.phylip 
# OUTPUT OPTIONAL alignment.fasta 
# OUTPUT OPTIONAL alignment.log
# PARAMETER OPTIONAL strategy: "Alignment strategy" TYPE [ mafft-fftns1: "FFT-NS-1 (Very fast. Recommended for more than 2000 sequences. progressive method\)", mafft-fftns: "FFT-NS-2 (Fast. Progressive method\) ", mafft-fftnsi: "FFT-NS-i (Slow. Iterative refinement method\) ", mafft-einsi: "E-INS-i (Very slow. recommended for less than 200 sequences with multiple conserved domains and long gaps\)", mafft-linsi: "L-INS-i (Very slow. Recommended for less than 200 sequences with one conserved domain and long gaps\)", mafft-ginsi: "G-INS-i (Very slow. Recommended for less than 200 sequences with global homology\)", "auto": "Automatic (FFT-NS-1, FFT-NS-2, FFT-NS-i or L-INS-i. Depends on data size\)"] DEFAULT "auto" (Select the alignment strategy. )
# PARAMETER OPTIONAL outformat: "Output format" TYPE [ aln: "Clustal", fasta: "Fasta" ] DEFAULT aln (Multiple sequence alingment file format)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: yes, no: no] DEFAULT no (Collect a log file.)


source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("sequence")

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")
mafft.path <- file.path(chipster.tools.path, "mafft" ,"bin")
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

if (num.queryseq > 50000){
	stop(paste('CHIPSTER-NOTE: Too many query sequences. Maximun is 50000 but your file contains ', num.queryseq ))
}


if ( str.filetype != "fasta"){
	seqret.binary <- file.path(emboss.path,"seqret")
	seqret.command <- paste(seqret.binary, "sequence sequence.fasta -auto")
	system(seqret.command)
	system("rm -f sequence")
	system("mv sequence.fasta sequence")
}


mafft.binary <- file.path(mafft.path, "mafft")

if ( strategy == "mafft-fftns1" ){
	mafft.options <- paste("--retree 1")
}

if ( strategy == "mafft-fftns" ){
	mafft.options <- paste("--retree 2")
}

if ( strategy == "mafft-fftnsi" ){
	mafft.options <- paste("--maxiterate 2")
}

if ( strategy == "mafft-einsi" ){
	mafft.options <- paste("--genafpair --maxiterate 1000")
}

if ( strategy == "mafft-linsi" ){
	mafft.options <- paste("--localpair --maxiterate 1000")
}

if ( strategy == "mafft-ginsi" ){
	mafft.options <- paste("--globalpair --maxiterate 1000")
}

if ( strategy == "auto" ){
	mafft.options <- paste("--auto")
}

outfile <- ("alignment.fasta")

if (outformat == "aln" ){
	mafft.options <- paste(mafft.options, "--clustalout")
	outfile <- ("alignment.aln.txt")
}

mafft.options <- paste(mafft.options, "--thread", chipster.threads.max)

command.full <- paste(mafft.binary, mafft.options, "sequence >", outfile, "2>alignment.log" )
system(command.full)

if ( save_log == "no") {
	system ("rm -f alignment.log")
}
