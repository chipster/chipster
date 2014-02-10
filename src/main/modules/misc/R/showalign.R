# TOOL showalign.R: "Display a multiple sequence alignment as formatted text" (Display a multiple sequence alignment in pretty format)
# INPUT OPTIONAL sequence: sequence TYPE GENERIC 
# OUTPUT OPTIONAL alignment.html
# OUTPUT OPTIONAL showalign.log
# PARAMETER OPTIONAL refseq: "The number or the name of the reference sequence" TYPE STRING DEFAULT 0 (If you give the number in the alignment or the name of a sequence, it will be taken to be the reference sequence. The reference sequence is always shown in full and is the one against which all the other sequences are compared. If this is set to 0 then the consensus sequence will be used as the reference sequence. By default the consensus sequence is used as the reference sequence.)
# PARAMETER OPTIONAL bottom: "Display the reference sequence at the bottom" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (If this is true then the reference sequence is displayed at the bottom of the alignment instead of the top.)
# PARAMETER OPTIONAL show: "What to show" TYPE [A: "All of the sequences", I: "Identities between the sequences", N: "Non-identities between the sequences", S: "Similarities between the sequences", D: "Dissimilarities between the sequences"] FROM 1 TO 1 DEFAULT N (What to show)
# PARAMETER OPTIONAL order: "Output order of the sequences" TYPE [I: "Input order - no change", A: "Alphabetical order of the names", S: "Similarity to the reference sequence"] FROM 1 TO 1 DEFAULT I (Output order of the sequences)
# PARAMETER OPTIONAL similarcase: "Show similar residues in lower-case" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (If this is set True, then when -show is set to 'Similarities' or 'Non-identities' and a residue is similar but not identical to the reference sequence residue, it will be changed to lower-case. If -show is set to 'All' then non-identical, non-similar residues will be changed to lower-case. If this is False then no change to the case of the residues is made on the basis of their similarity to the reference sequence.)
# PARAMETER OPTIONAL consensus: "Display the consensus line" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (If this is true then the consensus line is displayed.)
# PARAMETER OPTIONAL uppercase: "Regions to put in uppercase (eg: 4-57,78-94\)" TYPE STRING (Regions to put in uppercase. \ If this is left blank, then the sequence case is left alone. \ A set of regions is specified by a set of pairs of positions. \ The positions are integers. \ They are separated by any non-digit, non-alpha character. \ Examples of region specifications are: \ 24-45, 56-78 \ 1:45, 67=99;765..888 \ 1,5,8,10,23,45,57,99)
# PARAMETER OPTIONAL number: "Number the sequences" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (If this option is true then a line giving the positions in the alignment is displayed every 10 characters above the alignment.)
# PARAMETER OPTIONAL ruler: "Display ruler" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (If this option is true then a ruler line marking every 5th and 10th character in the alignment is displayed.)
# PARAMETER OPTIONAL width: "Width of sequence to display" TYPE INTEGER FROM 1 DEFAULT 100000 (Width of sequence to display)
# PARAMETER OPTIONAL margin: "Length of margin for sequence names" TYPE INTEGER FROM -1 DEFAULT -1 (This sets the length of the left-hand margin for sequence names. If the margin is set at 0 then no margin and no names are displayed. If the margin is set to a value that is less than the length of a sequence name then the sequence name is displayed truncated to the length of the margin. If the margin is set to -1 then the minimum margin width that will allow all the sequence names to be displayed in full plus a space at the end of the name will automatically be selected.)
# PARAMETER OPTIONAL html: "Use HTML formatting" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Use HTML formatting)
# PARAMETER OPTIONAL highlight: "Regions to colour in HTML (eg: 4-57 red 78-94 green\)" TYPE STRING (regions to colour if formatting for html. \ if this is left blank, then the sequence is left alone. \ a set of regions is specified by a set of pairs of positions. \ the positions are integers. \ they are followed by any valid html font colour. \ examples of region specifications are: \ 24-45 blue 56-78 orange \ 1-100 green 120-156 red \ a file of ranges to colour (one range per line\) can be specified as '@filename'.)
# PARAMETER OPTIONAL plurality: "Plurality check % for consensus" TYPE DECIMAL FROM 0.0 TO 100.0 DEFAULT 50.0 (Set a cut-off for the % of positive scoring matches below which there is no consensus. The default plurality is taken as 50% of the total weight of all the sequences in the alignment.)
# PARAMETER OPTIONAL identity: "Required % of identities at a position for consensus" TYPE DECIMAL FROM 0.0 TO 100.0 DEFAULT 0.0 (Provides the facility of setting the required number of identities at a position for it to give a consensus. Therefore, if this is set to 100% only columns of identities contribute to the consensus.)
# PARAMETER OPTIONAL gaps: "Use gap characters in consensus" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (If this option is true then gap characters can appear in the consensus. The alternative is 'N' for nucleotide, or 'X' for protein)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")
options(scipen=999)

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

if (num.queryseq > 1000){
	stop(paste('CHIPSTER-NOTE: Too many query sequences. Maximun is 1000 but your file contains ', num.queryseq ))
}


emboss.binary <- file.path(emboss.path, "showalign")
emboss.parameters <- paste("-auto")
emboss.parameters <- paste(emboss.parameters, "-sequence sequence")
emboss.parameters <- paste(emboss.parameters, "-outfile alignment.html")
emboss.parameters <- paste(emboss.parameters, "-refseq", refseq)
emboss.parameters <- paste(emboss.parameters, "-bottom", bottom)
emboss.parameters <- paste(emboss.parameters, "-show", show)
emboss.parameters <- paste(emboss.parameters, "-order", order)
emboss.parameters <- paste(emboss.parameters, "-similarcase", similarcase)
emboss.parameters <- paste(emboss.parameters, "-consensus", consensus)
if (nchar(uppercase) > 0 ) {
	emboss.parameters <- paste(emboss.parameters, '-uppercase "', uppercase, '"')
}
if (nchar(highlight) > 0 ) {
	emboss.parameters <- paste(emboss.parameters, '-highlight "', highlight,'"')
}
emboss.parameters <- paste(emboss.parameters, "-number", number)
emboss.parameters <- paste(emboss.parameters, "-ruler", ruler)
emboss.parameters <- paste(emboss.parameters, "-width", width)
emboss.parameters <- paste(emboss.parameters, "-margin", margin)
emboss.parameters <- paste(emboss.parameters, "-html", html)
emboss.parameters <- paste(emboss.parameters, "-plurality", plurality)
emboss.parameters <- paste(emboss.parameters, "-identity", identity)
emboss.parameters <- paste(emboss.parameters, "-gaps", gaps)

command.full <- paste(emboss.binary, emboss.parameters, ' >> showalign.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> showalign.log' )
system(echo.command)

system(command.full)

if ( save_log == "no") {
	system ("rm -f showalign.log")
}
