# TOOL getorf.R: "Find open reading farmes" (Finds and extracts open reading frames (ORFs\))
# INPUT sequence: "Input sequence" TYPE GENERIC 
# OUTPUT OPTIONAL getorf.fasta
# OUTPUT OPTIONAL getorf.log 
# PARAMETER OPTIONAL table: "Code to use" TYPE [0: Standard, 1: "Standard (with alternative initiation codons\)", 2: "Vertebrate Mitochondrial", 3: "Yeast Mitochondrial", 4: "Mold, Protozoan, Coelenterate Mitochondrial and Mycoplasma/Spiroplasma", 5: "Invertebrate Mitochondrial", 6: "Ciliate  Macronuclear and Dasycladacean", 9: "Echinoderm Mitochondrial", 10: "Euplotid Nuclear", 11: Bacterial, 12: "Alternative Yeast Nuclear", 13: "Ascidian Mitochondrial", 14: "Flatworm Mitochondrial", 15: "Blepharisma Macronuclear", 16: "Chlorophycean Mitochondrial", 21: "Trematode Mitochondrial", 22: "Scenedesmus obliquus", 23: "Thraustochytrium Mitochondrial"] FROM 1 TO 1 DEFAULT 0 (Code to use)
# PARAMETER OPTIONAL minsize: "Minimum nucleotide size of ORF to report" TYPE INTEGER DEFAULT 30 (Minimum nucleotide size of ORF to report)
# PARAMETER OPTIONAL maxsize: "Maximum nucleotide size of ORF to report" TYPE INTEGER DEFAULT 1000000 (Maximum nucleotide size of ORF to report)
# PARAMETER OPTIONAL find: "Type of output" TYPE [0: "Translation of regions between STOP codons", 1: "Translation of regions between START and STOP codons", 2: "Nucleic sequences between STOP codons", 3: "Nucleic sequences between START and STOP codons", 4: "Nucleotides flanking START codons", 5: "Nucleotides flanking initial STOP codons", 6: "Nucleotides flanking ending STOP codons"] FROM 1 TO 1 DEFAULT 0 (This is a small menu of possible output options. The first four options are to select either the protein translation or the original nucleic acid sequence of the open reading frame. There are two possible definitions of an open reading frame: it can either be a region that is free of STOP codons or a region that begins with a START codon and ends with a STOP codon. The last three options are probably only of interest to people who wish to investigate the statistical properties of the regions around potential START or STOP codons. The last option assumes that ORF lengths are calculated between two STOP codons.)
# PARAMETER OPTIONAL methionine: "Change initial START codons to Methionine" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (START codons at the beginning of protein products will usually code for Methionine, despite what the codon will code for when it is internal to a protein. This qualifier sets all such START codons to code for Methionine by default.)
# PARAMETER OPTIONAL circular: "Is the sequence circular" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Is the sequence circular)
# PARAMETER OPTIONAL reverse: "Find ORFs in the reverse sequence" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Set this to be false if you do not wish to find ORFs in the reverse complement of the sequence.)
# PARAMETER OPTIONAL flanking: "Number of flanking nucleotides to report" TYPE INTEGER DEFAULT 100 (If you have chosen one of the options of the type of sequence to find that gives the flanking sequence around a STOP or START codon, this allows you to set the number of nucleotides either side of that codon to output. If the region of flanking nucleotides crosses the start or end of the sequence, no output is given for this codon.)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)


# KM 8.11. 2013
options(scipen=999)
emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

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

if (num.queryseq > 50000){
	stop(paste('CHIPSTER-NOTE: Too many query sequences. Maximun is 50000 but your file contains ', num.queryseq ))
}

emboss.binary <- file.path(emboss.path, "getorf")
emboss.parameters <- paste('sequence -auto -outseq getorf.fasta')

emboss.parameters <- paste(emboss.parameters, "-table", table)
emboss.parameters <- paste(emboss.parameters, "-minsize", minsize)
emboss.parameters <- paste(emboss.parameters, "-maxsize", maxsize)
emboss.parameters <- paste(emboss.parameters, "-find", find)
emboss.parameters <- paste(emboss.parameters, "-methionine", methionine)
emboss.parameters <- paste(emboss.parameters, "-circular", circular)
emboss.parameters <- paste(emboss.parameters, "-reverse", reverse)
emboss.parameters <- paste(emboss.parameters, "-flanking", flanking)


command.full <- paste(emboss.binary, emboss.parameters, ' >> getorf.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> getorf.log' )
system(echo.command)

system(command.full)

if ( save_log == "no") {
	system ("rm -f getorf.log")
}
