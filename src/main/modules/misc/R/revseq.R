# TOOL revseq.R: "Reverse sequence" (Reverse and complement a nucleotide sequence.)
# INPUT sequence: "Query sequences" TYPE GENERIC
# OUTPUT reverse.fasta
# PARAMETER OPTIONAL reverse: "Reverse sequence" TYPE [ Y: Yes, N: No] DEFAULT Y (Set this to be false if you do not wish to reverse the output sequence)
# PARAMETER OPTIONAL complement: "Complement sequence" TYPE [ Y: Yes, N: No] DEFAULT Y (Set this to be false if you do not wish to complement the output sequence)

# KM 8.11. 2013
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("sequence")

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

if (num.queryseq > 50000){
	stop(paste('CHIPSTER-NOTE: Too many query sequences. Maximun is 50000 but your file contains ', num.queryseq ))
}

emboss.binary <- file.path(emboss.path, "revseq")
command.full <- paste(emboss.binary, 'sequence -outseq reverse.fasta -reverse', reverse ,'-complement', complement )
system(command.full)

