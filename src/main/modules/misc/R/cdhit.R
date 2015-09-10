# TOOL cdhit.R: "remove redundant sequences" (Cluster the given sequences and create a non-redundant sequence set)
# INPUT sequence: "Input sequence" TYPE GENERIC 
# OUTPUT OPTIONAL non_redundant.fasta
# OUTPUT OPTIONAL sequence_clusters.txt
# OUTPUT OPTIONAL remove_redundant.log 
# PARAMETER thresholdp: "The percentage sequence identity redundancy threshold." TYPE DECIMAL DEFAULT 95.0 (This option specifies the percentage sequence identity redundancy threshold. The percentage sequence identity redundancy threshold  determines the redundancy calculation. If a pair of proteins achieve greater than this threshold the shortest sequence is discarded.)
# PARAMETER OPTIONAL inptype: "Input sequence type" TYPE [prot: Protein, nuc: Nucloetide] DEFAULT prot (Define wether the input sequences are proteins or nucleotides.)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)


# KM 25.8. 2015
options(scipen=999)
emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")
cdhit.path <- file.path(chipster.tools.path, "cdhit")

source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("sequence")


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

if (num.queryseq > 100000) {
	stop(paste('CHIPSTER-NOTE: Too many query sequences. Maximun is 100000 but your file contains ', num.queryseq ))
}


threshold <- (thresholdp / 100)

if ( inptype == "prot"){
  cdhit.binary <- file.path(cdhit.path, "cd-hit")
  if ( threshold > 0.4 ){
	wlen <- (2)
  }

  if ( threshold > 0.5 ){
	wlen <- (3)	
  }

  if ( threshold > 0.6 ){
	wlen <- (4)
  }

  if ( threshold > 0.7 ){
	wlen <- (5)
  }
}

if ( inptype == "nuc"){
	cdhit.binary <- file.path(cdhit.path, "cd-hit-est")
	if ( threshold >= 0.74){
		wlen <- (4)
	}
	
	if ( threshold >= 0.8 ){
		wlen <- (5)	
	}
	
	if ( threshold >= 0.85 ){
		wlen <- (6)
	}
	
	if ( threshold >= 0.88 ){
		wlen <- (7)
	}
	if ( threshold >= 0.9 ){
		wlen <- (9)
	}	
}


cdhit.parameters <- paste('-i sequence -o non_redundant.fasta -M 2000 -d 0  -T ', chipster.threads.max )
cdhit.parameters <- paste(cdhit.parameters, "-n", wlen)
cdhit.parameters <- paste(cdhit.parameters, "-c", threshold)

command.full <- paste(cdhit.binary, cdhit.parameters, ' >> remove_redundant.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> remove_redundant.log' )
system(echo.command)

system(command.full)
system('mv non_redundant.fasta.clstr sequence_clusters.txt')
system('ls -l >> remove_redundant.log')

if ( save_log == "no") {
	system ("rm -f remove_redundant.log")
}
