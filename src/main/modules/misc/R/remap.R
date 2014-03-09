# TOOL remap.R: "Display restriction enzyme binding sites" (Display restriction enzyme binding sites in a nucleotide sequence)
# INPUT sequence: "Input sequence" TYPE GENERIC 
# OUTPUT OPTIONAL remap.html 
# OUTPUT OPTIONAL remap.log 
# PARAMETER enzymes: "Comma separated enzyme list" TYPE STRING DEFAULT all (Rhe name 'all' reads in all enzyme names from the rebase database. You can specify enzymes by giving their names with commas between then, such as: 'hincii,hinfi,ppii,hindiii'.  The case of the names is not important.)
# PARAMETER sitelen: "Minimum recognition site length" TYPE INTEGER FROM 2 TO 20 DEFAULT 4 (This sets the minimum length of the restriction enzyme recognition site. Any enzymes with sites shorter than this will be ignored.)
# PARAMETER OPTIONAL mincuts: "Minimum cuts per RE" TYPE INTEGER FROM 1 TO 1000 DEFAULT 1 (This sets the minimum number of cuts for any restriction enzyme that will be considered. Any enzymes that cut fewer times than this will be ignored.)
# PARAMETER OPTIONAL maxcuts: "Maximum cuts per RE" TYPE INTEGER FROM 1 DEFAULT 2000000000 (This sets the maximum number of cuts for any restriction enzyme that will be considered. Any enzymes that cut more times than this will be ignored.)
# PARAMETER OPTIONAL single: "Force single site only cuts" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (If this is set then this forces the values of the mincuts and maxcuts qualifiers to both be 1. Any other value you may have set them to will be ignored.)
# PARAMETER OPTIONAL blunt: "Allow blunt end cutters" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (This allows those enzymes which cut at the same position on the forward and reverse strands to be considered.)
# PARAMETER OPTIONAL sticky: "Allow sticky end cutters" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (This allows those enzymes which cut at different positions on the forward and reverse strands, leaving an overhang, to be considered.)
# PARAMETER OPTIONAL ambiguity: "Allow ambiguous matches" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (This allows those enzymes which have one or more 'N' ambiguity codes in their pattern to be considered)
# PARAMETER OPTIONAL plasmid: "Allow circular DNA" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (If this is set then this allows searches for restriction enzyme recognition site and cut positions that span the end of the sequence to be considered.)
# PARAMETER OPTIONAL methylation: "Use methylation data" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (If this is set then RE recognition sites will not match methylated bases.)
# PARAMETER OPTIONAL commercial: "Only enzymes with suppliers" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (If this is set, then only those enzymes with a commercial supplier will be searched for. This qualifier is ignored if you have specified an explicit list of enzymes to search for, rather than searching through 'all' the enzymes in the REBASE database. It is assumed that, if you are asking for an explicit enzyme, then you probably know where to get it from and so all enzymes names that you have asked to be searched for, and which cut, will be reported whether or not they have a commercial supplier.)
# PARAMETER OPTIONAL table: "Genetic code to use" TYPE [0: Standard, 1: "Standard (with alternative initiation codons\)", 2: "Vertebrate Mitochondrial", 3: "Yeast Mitochondrial", 4: "Mold, Protozoan, Coelenterate Mitochondrial and Mycoplasma/Spiroplasma", 5: "Invertebrate Mitochondrial", 6: "Ciliate Macronuclear and Dasycladacean", 9: "Echinoderm Mitochondrial", 10: "Euplotid Nuclear", 11: Bacterial, 12: "Alternative Yeast Nuclear", 13: "Ascidian Mitochondrial", 14: "Flatworm Mitochondrial", 15: "Blepharisma Macronuclear", 16: "Chlorophycean Mitochondrial", 21: "Trematode Mitochondrial", 22: "Scenedesmus obliquus", 23: "Thraustochytrium Mitochondrial"] FROM 1 TO 1 DEFAULT 0 (Genetic code to use)
# PARAMETER OPTIONAL frame: "Frame(s\) to translate" TYPE [1: 1, 2: 2, 3: 3, F: "Forward three frames", -1: -1, -2: -2, -3: -3, R: "Reverse three frames", 6: "All six frames"] FROM 1 TO 6 DEFAULT 6 (This allows you to specify the frames that are translated. If you are not displaying cut sites on the reverse sense, then the reverse sense translations will not be displayed even if you have requested frames 4, 5 or 6. By default, all six frames will be displayed.)
# PARAMETER OPTIONAL cutlist: "List the enzymes that cut" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (This produces lists in the output of the enzymes that cut, those that cut but are excluded because that cut fewer times than mincut or more times than maxcut and those enzymes that do not cut.)
# PARAMETER OPTIONAL flatreformat: "Display RE sites in flat format" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (This changes the output format to one where the recognition site is indicated by a row of '===' characters and the cut site is pointed to by a '>' character in the forward sense, or a '<' in the reverse sense strand.)
# PARAMETER OPTIONAL limit: "Limits reports to one isoschizomer" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (This limits the reporting of enzymes to just one enzyme from each group of isoschizomers. The enzyme chosen to represent an isoschizomer group is the prototype indicated in the data file 'embossre.equ', which is created by the program 'rebaseextract'. If you prefer different prototypes to be used, make a copy of embossre.equ in your home directory and edit it. If this value is set to be false then all of the input enzymes will be reported. You might like to set this to false if you are supplying an explicit set of enzymes rather than searching 'all' of them.)
# PARAMETER OPTIONAL translation: "Display translation" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (This displays the 6-frame translations of the sequence in the output.)
# PARAMETER OPTIONAL reverse: "Display cut sites and translation of reverse sense" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (This displays the cut sites and translation of the reverse sense.)
# PARAMETER OPTIONAL orfminsize: "Minimum size of ORFs" TYPE INTEGER FROM 0 DEFAULT 0 (This sets the minimum size of Open Reading Frames (ORFs\) to display in the translations. All other translation regions are masked by changing the amino acids to '-' characters.)
# PARAMETER OPTIONAL uppercase: "Regions to put in uppercase (eg: 4-57,78-94\)" TYPE STRING (Regions to put in uppercase. \ If this is left blank, then the sequence case is left alone. \ A set of regions is specified by a set of pairs of positions. \ The positions are integers. \ They are separated by any non-digit, non-alpha character. \ Examples of region specifications are: \ 24-45, 56-78 \ 1:45, 67=99;765..888 \ 1,5,8,10,23,45,57,99)
# PARAMETER OPTIONAL highlight: "Regions to colour in HTML (eg: 4-57 red 78-94 green\)" TYPE STRING (regions to colour if formatting for html. \ if this is left blank, then the sequence is left alone. \ a set of regions is specified by a set of pairs of positions. \ the positions are integers. \ they are followed by any valid html font colour. \ examples of region specifications are: \ 24-45 blue 56-78 orange \ 1-100 green 120-156 red.)
# PARAMETER OPTIONAL threeletter: "Display protein sequences in three-letter code" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Display protein sequences in three-letter code)
# PARAMETER OPTIONAL number: "Number the sequences" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Number the sequences)
# PARAMETER OPTIONAL width: "Width of sequence to display" TYPE INTEGER FROM 1 DEFAULT 100000 (Width of sequence to display)
# PARAMETER OPTIONAL length: "Line length of page (0 for indefinite\)" TYPE INTEGER FROM 0 DEFAULT 0 (Line length of page (0 for indefinite\))
# PARAMETER OPTIONAL margin: "Margin around sequence for numbering" TYPE INTEGER FROM 0 DEFAULT 10 (Margin around sequence for numbering)
# PARAMETER OPTIONAL name: "Display sequence ID" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Set this to be false if you do not wish to display the ID name of the sequence)
# PARAMETER OPTIONAL description: "Display description" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Set this to be false if you do not wish to display the description of the sequence)
# PARAMETER OPTIONAL offset: "Offset to start numbering the sequence from" TYPE INTEGER DEFAULT 1 (Offset to start numbering the sequence from)
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

if (num.queryseq > 100 ){
	stop(paste('CHIPSTER-NOTE: Too many query sequences. Maximun is 100 but your file contains ', num.queryseq ))
}

seqlen.exe <- file.path(emboss.path, "infoseq -nohead -only -length  sequence -filter | sort -n |tail -1")
str.seqlen <- system(seqlen.exe, intern = TRUE )
num.seqlen <- as.integer(str.seqlen)

if ( num.seqlen <  width ){
	width <- num.seqlen 
}

emboss.binary <- file.path(emboss.path, "remap")
emboss.parameters <- paste('sequence -auto -outfile remap.html -html Y')
emboss.parameters <- paste(emboss.parameters, "-enzymes", enzymes)
emboss.parameters <- paste(emboss.parameters, "-sitelen", sitelen)
emboss.parameters <- paste(emboss.parameters, "-mincuts", mincuts)
emboss.parameters <- paste(emboss.parameters, "-maxcuts", maxcuts)
emboss.parameters <- paste(emboss.parameters, "-single", single)
emboss.parameters <- paste(emboss.parameters, "-blunt", blunt)
emboss.parameters <- paste(emboss.parameters, "-sticky", sticky)
emboss.parameters <- paste(emboss.parameters, "-ambiguity", ambiguity)
emboss.parameters <- paste(emboss.parameters, "-plasmid", plasmid)
emboss.parameters <- paste(emboss.parameters, "-methylation", methylation)
emboss.parameters <- paste(emboss.parameters, "-commercial", commercial)
emboss.parameters <- paste(emboss.parameters, "-table", table)
emboss.parameters <- paste(emboss.parameters, "-frame", frame)
emboss.parameters <- paste(emboss.parameters, "-cutlist", cutlist)
emboss.parameters <- paste(emboss.parameters, "-flatreformat", flatreformat)
emboss.parameters <- paste(emboss.parameters, "-limit", limit)
emboss.parameters <- paste(emboss.parameters, "-translation", translation)
emboss.parameters <- paste(emboss.parameters, "-reverse", reverse)
emboss.parameters <- paste(emboss.parameters, "-orfminsize", orfminsize)

if (nchar(uppercase) > 0 ) {
emboss.parameters <- paste(emboss.parameters, '-uppercase "', uppercase, '"')
}

if (nchar(highlight) > 0 ) {
emboss.parameters <- paste(emboss.parameters, '-highlight "', highlight,'"')
}

emboss.parameters <- paste(emboss.parameters, "-threeletter", threeletter)
emboss.parameters <- paste(emboss.parameters, "-number", number)
emboss.parameters <- paste(emboss.parameters, "-width", width)
emboss.parameters <- paste(emboss.parameters, "-length", length)
emboss.parameters <- paste(emboss.parameters, "-margin", margin)
emboss.parameters <- paste(emboss.parameters, "-name", name)
emboss.parameters <- paste(emboss.parameters, "-description", description)
emboss.parameters <- paste(emboss.parameters, "-offset", offset)



command.full <- paste(emboss.binary, emboss.parameters, ' >> remap.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> remap.log' )
system(echo.command)

system(command.full)

if ( save_log == "no") {
	system ("rm -f remap.log")
}
