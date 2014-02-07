# TOOL seqor.R: "Merge two sequence sets using logical operators" (This tool joins two sequence sets based selected logical operator. Note that this tool joins the two sequence sets based only on the sequence data. Names or ID:s are not used in comparison.)
# INPUT sequences1.txt: "First sequence set" TYPE GENERIC
# INPUT sequences2.txt: "Second sequence set" TYPE GENERIC
# OUTPUT OPTIONAL results.fasta
# OUTPUT OPTIONAL log.txt
# PARAMETER OPTIONAL operator: "Logical operator for merging" TYPE [OR: OR, AND: AND, XOR: XOR, NOT: NOT] DEFAULT OR (The following logical operators combine the sequences in the following ways. OR - gives all that occur in one set or the other. AND - gives only those which occur in both sets.  XOR - gives those which only occur in one set or the other, but not in both. NOT - gives those which occur in the first set except for those that also occur in the second )
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)



# KM 8.11. 2013

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

#check sequece file type
inputfile.to.check <- ("sequences1.txt")
sfcheck.binary <- file.path(chipster.module.path ,"/shell/sfcheck.sh")
sfcheck.command <- paste(sfcheck.binary, emboss.path, inputfile.to.check )
str.filetype <- system(sfcheck.command, intern = TRUE )

if ( str.filetype == "Not an EMBOSS compatible sequence file"){
	stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
}

#check sequece file type
inputfile.to.check <- ("sequences2.txt")
sfcheck.binary <- file.path(chipster.module.path ,"/shell/sfcheck.sh")
sfcheck.command <- paste(sfcheck.binary, emboss.path, inputfile.to.check )
str.filetype <- system(sfcheck.command, intern = TRUE )

if ( str.filetype == "Not an EMBOSS compatible sequence file"){
	stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
}

#count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount sequences1.txt -filter")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)

#round(num.queryseq)

if (num.queryseq > 100000){
	stop(paste("CHIPSTER-NOTE: Too many query sequences. Maximun is 100000 but your file contains ", num.queryseq ))
}

#count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount sequences2.txt -filter")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)

#round(num.queryseq)

if (num.queryseq > 100000){
	stop(paste("CHIPSTER-NOTE: Too many query sequences. Maximun is 100000 but your file contains ", num.queryseq ))
}


emboss.binary <- file.path(emboss.path, "listor")
command.full <- paste(emboss.binary, ' -firstsequences sequences1.txt -secondsequences sequences2.txt -outfile resultlist.txt -operator', operator,'  > log.txt 2>&1' )
system(command.full)

system("ls -l >> log.txt")

str.queryseq <- system("cat resultlist.txt | wc -l" )
num.queryseq <- as.integer(str.queryseq)

if (num.queryseq < 1){
	stop(paste("CHIPSTER-NOTE: Joining the two sequce set with opeator:", operator, "produced an empty sequence set." ))
}

emboss.binary <- file.path(emboss.path, "seqret")
command.full <- paste(emboss.binary, ' @resultlist.txt results.fasta >> log.txt 2>&1' )
system(command.full)

if ( save_log == "no" ){
	system("rm -f log.txt")
}

