# TOOL emboss.R: "Other sequence utilities" (Selection of simple EMBOSS sequence analysis tools.)
# INPUT input.txt: "Query sequences" TYPE GENERIC
# OUTPUT OPTIONAL output.txt
# OUTPUT OPTIONAL output.html
# OUTPUT OPTIONAL output.fasta
# OUTPUT OPTIONAL emboss.log
# PARAMETER OPTIONAL ecommand: "Select task" TYPE [chips: "(chips\) Calculate Nc codon usage statistic", degapseq: "(degapseq\) Remove non-alphabetic e.g. gap characters from sequences",  garnier: "(garnier\) Predict protein secondary structure using GOR method", geecee: "(geecee\) Calculate fractional GC content of nucleic acid sequences", helixturnhelix: "(helixturnhelix\) Identify nucleic acid-binding motifs in protein sequences", maskambignuc: "(maskambignuc\) Mask all ambiguity characters in nucleotide sequences with N", maskambigprot: "(maskambigprot\) Mask all ambiguity characters in protein sequences with X", prettyseq: "(prettyseq\) Write a nucleotide sequence and its translation to file", showseq: "(showseq\) Display nucleic sequences with features in pretty format", showpep: "(showpep\) Display protein sequences with features", showfeat: "(showfeat\) Display features of a sequence in pretty format", seqcount: "(seqcount\) Count sequences", sizeseq:  "(sizeseq\) Sort sequences by size"] DEFAULT seqcount (Choose EMBOSS commmand to execute) 
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)


# K.M 28.10.2013


# pb settings
use.html=(0)

if ( ecommand == "showseq" || ecommand == "showpep" || ecommand == "prettyseq" ){
	ecommand <- paste(ecommand,"-width 100000"  )
	use.html=(1)
}

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

emboss.binary <- file.path(emboss.path,ecommand)
command.full <- paste(emboss.binary, " input.txt -filter -auto > output.txt 2>&1" )

echo.command <- paste('echo "',command.full, ' "> emboss.log' )
system(echo.command)


system(command.full)

system("ls -l >> emboss.log")

#test.command <- paste('echo "',ecommand, ' ">> emboss.log')
#system(test.command)

if ( use.html == 1){
system("echo showseq >> emboss.log")
system("echo '<html>' > output.html")
system("echo '<body>' >> output.html")
system("echo '<pre>' >> output.html")
system("cat output.txt >> output.html")
system("echo '</pre>' >> output.html")
system("echo '</body>' >> output.html")
system("echo '</html>' >> output.html")
system("rm -f output.txt")
}


if ( ecommand == "degapseq" || ecommand == "maskambignuc" || ecommand == "maskambigprot"){
	system("mv output.txt output.fasta")
}


if ( save_log == "no") {
	system ("rm -f emboss.log")
}