# TOOL restrict.R: "Find restriction enzyme cleavage sites" (Report restriction enzyme cleavage sites in a nucleotide sequence)
# INPUT sequence: "Input sequence" TYPE GENERIC 
# OUTPUT OPTIONAL restrict.tsv 
# OUTPUT OPTIONAL restrict.txt 
# OUTPUT OPTIONAL restrict.log 
# PARAMETER sitelen: "Minimum recognition site length" TYPE INTEGER FROM 2 TO 20 DEFAULT 4 (This sets the minimum length of the restriction enzyme recognition site. Any enzymes with sites shorter than this will be ignored.)
# PARAMETER enzymes: "Comma separated enzyme list" TYPE STRING DEFAULT all (the name 'all' reads in all enzyme names from the rebase database. you can specify enzymes by giving their names with commas between then, such as: 'hincii,hinfi,ppii,hindiii'. \ the case of the names is not important. you can specify a file of enzyme names to read in by giving the name of the file holding the enzyme names with a '@' character in front of it, for example, '@enz.list'. \ blank lines and lines starting with a hash character or '!' are ignored and all other lines are concatenated together with a comma character ',' and then treated as the list of enzymes to search for. \ an example of a file of enzyme names is: \ ! my enzymes \ hincii, ppiii \ ! other enzymes \ hindiii \ hinfi \ ppii)
# PARAMETER OPTIONAL rformat: "Output format type" TYPE [excel: "Table", table: "Text formatted report", gff: "GFF3 formatted file", tagseq: "Tagseq format", listfile: "EMBOSS list file"] DEFAULT excel (Output format type)
# PARAMETER OPTIONAL mincuts: "Minimum cuts per RE" TYPE INTEGER FROM 1 TO 1000 DEFAULT 1 (This sets the minimum number of cuts for any restriction enzyme that will be considered. Any enzymes that cut fewer times than this will be ignored.)
# PARAMETER OPTIONAL maxcuts: "Maximum cuts per RE" TYPE INTEGER FROM 1 TO 2000000 DEFAULT 1000000 (This sets the maximum number of cuts for any restriction enzyme that will be considered. Any enzymes that cut more times than this will be ignored.)
# PARAMETER OPTIONAL solofragment: "List individual enzymes with their fragments" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (This gives the fragment lengths of the forward sense strand produced by complete restriction by each restriction enzyme on its own. Results are added to the tail section of the report.)
# PARAMETER OPTIONAL single: "Force single site only cuts" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (If this is set then this forces the values of the mincuts and maxcuts qualifiers to both be 1. Any other value you may have set them to will be ignored.)
# PARAMETER OPTIONAL blunt: "Allow blunt end cutters" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (This allows those enzymes which cut at the same position on the forward and reverse strands to be considered.)
# PARAMETER OPTIONAL sticky: "Allow sticky end cutters" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (This allows those enzymes which cut at different positions on the forward and reverse strands, leaving an overhang, to be considered.)
# PARAMETER OPTIONAL ambiguity: "Allow ambiguous matches" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (This allows those enzymes which have one or more 'N' ambiguity codes in their pattern to be considered)
# PARAMETER OPTIONAL plasmid: "Allow circular DNA" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (If this is set then this allows searches for restriction enzyme recognition site and cut positions that span the end of the sequence to be considered.)
# PARAMETER OPTIONAL methylation: "Use methylation data" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (If this is set then RE recognition sites will not match methylated bases.)
# PARAMETER OPTIONAL commercial: "Only enzymes with suppliers" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (If this is set, then only those enzymes with a commercial supplier will be searched for. This qualifier is ignored if you have specified an explicit list of enzymes to search for, rather than searching through 'all' the enzymes in the REBASE database. It is assumed that, if you are asking for an explicit enzyme, then you probably know where to get it from and so all enzymes names that you have asked to be searched for, and which cut, will be reported whether or not they have a commercial supplier.)
# PARAMETER OPTIONAL limit: "Limits reports to one isoschizomer" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (This limits the reporting of enzymes to just one enzyme from each group of isoschizomers. The enzyme chosen to represent an isoschizomer group is the prototype indicated in the data file 'embossre.equ', which is created by the program 'rebaseextract'. If you prefer different prototypes to be used, make a copy of embossre.equ in your home directory and edit it. If this value is set to be false then all of the input enzymes will be reported. You might like to set this to false if you are supplying an explicit set of enzymes rather than searching 'all' of them.)
# PARAMETER OPTIONAL alphabetic: "Sort output alphabetically" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Sort output alphabetically)
# PARAMETER OPTIONAL fragments: "Show fragment lengths" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (This gives the fragment lengths of the forward sense strand produced by complete restriction using all of the input enzymes together. Results are added to the tail section of the report.)
# PARAMETER OPTIONAL name: "Show sequence name" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Show sequence name)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)


# KM 8.11. 2013
# KM 8.11. 2013
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("sequence")

options(scipen=999)
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


emboss.binary <- file.path(emboss.path, "restrict")
emboss.parameters <- paste('sequence -auto -outfile restrict.txt -sitelen', sitelen, '-enzymes', enzymes)
emboss.parameters <- paste(emboss.parameters, '-rformat', rformat)
emboss.parameters <- paste(emboss.parameters, '-min', mincuts) 
emboss.parameters <- paste(emboss.parameters, '-max', maxcuts)
emboss.parameters <- paste(emboss.parameters, '-solofragment', solofragment )
emboss.parameters <- paste(emboss.parameters, '-single', single )
emboss.parameters <- paste(emboss.parameters, '-blunt', blunt )
emboss.parameters <- paste(emboss.parameters, '-sticky', sticky )
emboss.parameters <- paste(emboss.parameters, '-ambiguity', ambiguity )
emboss.parameters <- paste(emboss.parameters, '-plasmid', plasmid )
emboss.parameters <- paste(emboss.parameters, '-methylation', methylation )
emboss.parameters <- paste(emboss.parameters, '-commercial', commercial )
emboss.parameters <- paste(emboss.parameters, '-limit', limit )
emboss.parameters <- paste(emboss.parameters, '-alphabetic', alphabetic )
emboss.parameters <- paste(emboss.parameters, '-fragments', fragments )
emboss.parameters <- paste(emboss.parameters, '-name', name)

command.full <- paste(emboss.binary, emboss.parameters, ' >> restrict.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> restrict.log' )
system(echo.command)

system(command.full)

if (rformat == "excel") {
	system ("mv restrict.txt restrict.tsv")	
}
if ( save_log == "no") {
	system ("rm -f restrict.log")
}
