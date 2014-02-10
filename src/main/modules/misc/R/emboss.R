# TOOL emboss.R: "Other sequence utilities" (Selection of simple EMBOSS sequence analysis tools.)
# INPUT input.txt: "Query sequences" TYPE GENERIC
# OUTPUT output.txt
# PARAMETER OPTIONAL ecommand: "Select task" TYPE [chips: "(chips\) Calculate Nc codon usage statistic", degapseq: "(degapseq\) Remove non-alphabetic e.g. gap characters from sequences",  garnier: "(garnier\) Predict protein secondary structure using GOR method", geecee: "(geecee\) Calculate fractional GC content of nucleic acid sequences", helixturnhelix: "(helixturnhelix\) Identify nucleic acid-binding motifs in protein sequences", maskambignuc: "(maskambignuc\) Mask all ambiguity characters in nucleotide sequences with N", maskambigprot: "(maskambigprot\) Mask all ambiguity characters in protein sequences with X", prettyseq: "(prettyseq\) Write a nucleotide sequence and its translation to file", showseq: "(showseq\) Display nucleic sequences with features in pretty format", showpep: "(showpep\) Display protein sequences with features", showfeat: "(showfeat\) Display features of a sequence in pretty format", seqcount: "(seqcount\) Count sequences", sizeseq:  "(sizeseq\) Sort sequences by size"] DEFAULT seqcount (Choose EMBOSS commmand to execute) 


# K.M 28.10.2013


# pb settings

if ( ecommand == "showseq"){
	ecommand <- paste(ecommand,"-width 100000"  )
}

if ( ecommand == "showpep"){
	ecommand <- paste(ecommand,"-width 100000"  )
}

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

emboss.binary <- file.path(emboss.path,ecommand)
command.full <- paste(emboss.binary, " input.txt -filter -auto > output.txt 2>&1" )
system(command.full)
