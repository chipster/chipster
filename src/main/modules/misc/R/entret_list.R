# TOOL entret_list.R: "Retrieve sequences from sequence file" (Tool to retrieve a set of sequences from a given sequence set based on a list of sequence IDs or names)
# INPUT idlist.txt: "Names or IDs of the database entries" TYPE GENERIC
# INPUT datafile.txt: "Sequence file from which the sequeces will be picked" TYPE GENERIC
# OUTPUT OPTIONAL data.txt
# OUTPUT OPTIONAL data.fasta
# OUTPUT OPTIONAL data.html
# OUTPUT OPTIONAL log.txt


# KM 8.11. 2013

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

#check sequece file type
inputfile.to.check <- ("datafile.txt")
sfcheck.binary <- file.path(chipster.module.path ,"/shell/sfcheck.sh")
sfcheck.command <- paste(sfcheck.binary, emboss.path, inputfile.to.check )
str.filetype <- system(sfcheck.command, intern = TRUE )

if ( str.filetype == "Not an EMBOSS compatible sequence file"){
	stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
}

#count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount datafile.txt -filter")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)

#round(num.queryseq)

if (num.queryseq > 100000){
	stop(paste("CHIPSTER-NOTE: Too many query sequences. Maximun is 100000 but your file contains ", num.queryseq ))
}


list_command <- paste("awk '{ print \"datafile.txt:\"$1 }' idlist.txt > idusa.txt" ,sep="" )
#stop("CHIPSTER-NOTE:",list_command)
system(list_command)

emboss.binary <- file.path(emboss.path, "entret")
command.full <- paste(emboss.binary, ' @idusa.txt -outfile data.txt > log.txt 2>&1' )
system(command.full)

log.row.num <- system("cat log.txt | wc -l ", intern = TRUE )

if ( log.row.num < 2 ){
	system("rm -f log.txt")
}

