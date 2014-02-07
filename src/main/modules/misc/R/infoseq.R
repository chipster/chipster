# TOOL infoseq.R: "Display basic information about sequences" (Tool to calculate basic properties of a sequence file.)
# INPUT sequence: "Query sequences" TYPE GENERIC
# OUTPUT outfile.infoseq.tsv


# KM 8.11. 2013

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


infoseq.binary <- file.path(emboss.path, "infoseq")
command.full <- paste(infoseq.binary, 'sequence -outfile outfile.infoseq.tsv -nodatabase -nousa -nocolumns -delimiter "\t"' )
system(command.full)

