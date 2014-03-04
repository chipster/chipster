#TOOL transeq.R: Translate  (Translate nucleic acid sequences)
#INPUT sequence TYPE GENERIC 
#OUTPUT OPTIONAL transeq.fasta
#OUTPUT OPTIONAL log.txt
#PARAMETER OPTIONAL frame: "Frame(s\) to translate" TYPE [1: 1, 2: 2, 3: 3, F: "Forward three frames", -1: -1, -2: -2, -3: -3, R: "Reverse three frames", 6: "All six frames"] FROM 1 TO 6 DEFAULT 1 (Frame(s\) to translate)
#PARAMETER OPTIONAL table: "Code to use" TYPE [0: Standard, 1: "Standard (with alternative initiation odons\)", 2: "Vertebrate Mitochondrial", 3: "Yeast Mitochondrial", 4: "Mold, Protozoan, Coelenterate Mitochondrial and Mycoplasma/Spiroplasma", 5: "Invertebrate Mitochondrial", 6: "Ciliate Macronuclear and Dasycladacean", 9: "Echinoderm Mitochondrial", 10: "Euplotid Nuclear", 11: Bacterial, 12: "Alternative Yeast Nuclear", 13: "Ascidian Mitochondrial", 14: "Flatworm Mitochondrial", 15: "Blepharisma Macronuclear", 16: "Chlorophycean Mitochondrial", 21: "Trematode Mitochondrial", 22: "Scenedesmus obliquus", 23: "Thraustochytrium Mitochondrial"] DEFAULT 0 (Code to use)
#PARAMETER OPTIONAL regions: "Regions to translate (eg: 4-57,78-94\)" TYPE STRING DEFAULT "full length" (Regions to translate. \ A set of regions is specified by a set of pairs of positions. \ The positions are integers. \ They are separated by any non-digit, non-alpha character. \ Examples of region specifications are: \ 24-45, 56-78 \ 1:45, 67=99;765..888 \ 1,5,8,10,23,45,57,99 \ Note: you should not try to use this option with any other frame than the default, -frame=1)
#PARAMETER OPTIONAL trim: "Trim trailing X's and *'s" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (This removes all 'X' and '*' characters from the right end of the translation. The trimming process starts at the end and continues until the next character is not a 'X' or a '*')
#PARAMETER OPTIONAL clean: "Change all *'s to X's" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (This changes all STOP codon positions from the '*' character to 'X' (an unknown residue\). This is useful because some programs will not accept protein sequences with '*' characters in them.)
#PARAMETER OPTIONAL alternative: "Define frame '-1' as starting in the last codon" TYPE [Y: Yes, N: No] DEFAULT N (The default definition of frame '-1' is the reverse-complement of the set of codons used in frame 1. (Frame -2 is the set of codons used by frame 2, similarly frames -3 and 3\). This is a common standard, used by the Staden package and other programs. If you prefer to define frame '-1' as using the set of codons starting with the last codon of the sequence, then set this to be true.)



# KM 8.11. 2013

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

emboss.parameters <- paste('sequence -outseq transeq.fasta -frame', frame, '-table', table, "-trim", trim, "-clean", clean, "-alternative", alternative)

emboss.binary <- file.path(emboss.path, "transeq")
if ( regions !=  "full length"){
	emboss.parameters <- paste(emboss.parameters, "-regions", regions )
}

command.full <- paste(emboss.binary, emboss.parameters,"-auto > log.txt 2>&1")
system(command.full)

log.row.num <- system("cat log.txt | wc -l ", intern = TRUE )

if ( log.row.num < 2 ){
	system("rm -f log.txt")
}
