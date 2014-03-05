# TOOL cons.R: "Create a consensus sequence" (Create a consensus sequence from a multiple sequence alignment)
# INPUT OPTIONAL sequence: sequence TYPE GENERIC 
# OUTPUT OPTIONAL consensus.fasta
# OUTPUT OPTIONAL cons.log
# PARAMETER OPTIONAL ctype: "Consensus sequence type" TYPE [cons: "Normal", consambig:"Ambiguous consensus sequence"] DEFAULT cons (Type of the consensus sequence to calculate. In the case of Ambiguous consensus sequence, EMBOSS command consambig is used in stead of cons. In the case of Ambiguous consensus sequence the other paramters of this tool are ignored )
# PARAMETER OPTIONAL plurality: "Plurality check value" TYPE DECIMAL (Set a cut-off for the number of positive matches below which there is no consensus. The default plurality is taken as half the total weight of all the sequences in the alignment.)
# PARAMETER OPTIONAL identity: "Required number of identities at a position" TYPE INTEGER FROM 0 DEFAULT 0 (Provides the facility of setting the required number of identities at a site for it to give a consensus at that position. Therefore, if this is set to the number of sequences in the alignment only columns of identities contribute to the consensus.)
# PARAMETER OPTIONAL name: "Name of the consensus sequence" TYPE STRING DEFAULT Consensus (Name of the consensus sequence)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("sequence")

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

if (num.queryseq > 10000){
	stop(paste("CHIPSTER-NOTE: Too many query sequences. Maximun is 10000 but your file contains ", num.queryseq ))
}

if (ctype == "cons"){

  emboss.binary <- file.path(emboss.path, "cons")
  emboss.parameters <- paste("-auto")
  emboss.parameters <- paste(emboss.parameters, "-sequence sequence")
  emboss.parameters <- paste(emboss.parameters, "-outseq consensus.fasta")
  if (is.numeric(plurality)){
  emboss.parameters <- paste(emboss.parameters, "-plurality", plurality)
  }
  emboss.parameters <- paste(emboss.parameters, "-identity", identity)
  if (nchar(name) > 0 ) {
     emboss.parameters <- paste(emboss.parameters, "-name", name)
  }
  
}

if (ctype == "consambig"){
	
	emboss.binary <- file.path(emboss.path, "consambig")
	emboss.parameters <- paste("-auto")
	emboss.parameters <- paste(emboss.parameters, "-sequence sequence")
	emboss.parameters <- paste(emboss.parameters, "-outseq consensus.fasta")
	if (nchar(name) > 0 ) {
		emboss.parameters <- paste(emboss.parameters, "-name", name)
	}
	
}	

command.full <- paste(emboss.binary, emboss.parameters, ' >> cons.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> cons.log' )	

system(echo.command)

system(command.full)


if ( save_log == "no") {
	system ("rm -f cons.log")
}
