# TOOL prettyplot.R: "Plot out sequence alignment" (Draw a sequence alignment with pretty formatting)
# INPUT sequences: sequences TYPE GENERIC 
# OUTPUT OPTIONAL prettyplot{...}.png 
# OUTPUT OPTIONAL prettyplot.log
# PARAMETER OPTIONAL residuesperline: "Number of residues to be displayed on each line" TYPE INTEGER DEFAULT 50 (The number of residues to be displayed on each line)
# PARAMETER OPTIONAL blocksperline: "Blocks of residues on each line" TYPE INTEGER FROM 1 TO 50 DEFAULT 1 (Blocks of residues on each line)
# PARAMETER OPTIONAL ccolours: "Colour residues by their consensus value." TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Colour residues by their consensus value.)
# PARAMETER OPTIONAL cidentity: "Colour to display identical residues (RED\)" TYPE STRING DEFAULT RED (Colour to display identical residues (RED\))
# PARAMETER OPTIONAL csimilarity: "Colour to display similar residues (GREEN\)" TYPE STRING DEFAULT GREEN (Colour to display similar residues (GREEN\))
# PARAMETER OPTIONAL cother: "Colour to display other residues (BLACK\)" TYPE STRING DEFAULT BLACK (Colour to display other residues (BLACK\))
# PARAMETER OPTIONAL docolour: "Colour residues by table oily, amide etc." TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Colour residues by table oily, amide etc.)
# PARAMETER OPTIONAL shade: Shading TYPE STRING DEFAULT BPLW (Set to BPLW for normal shading \ (black, pale, light, white\) \ so for pair = 1.5,1.0,0.5 and shade = BPLW \ Residues score Colour \ 1.5 or over... BLACK (B\) \ 1.0 to 1.5 ... BROWN (P\) \ 0.5 to 1.0 ... WHEAT (L\) \ under 0.5 .... WHITE (W\) \ The only four letters allowed are BPLW, in any order.)
# PARAMETER OPTIONAL pair: "Values to represent identical similar related" TYPE STRING FROM 0.0 DEFAULT "1.5,1.0,0.5" (Values to represent identical similar related)
# PARAMETER OPTIONAL identity: "Only match those which are identical in all sequences." TYPE INTEGER FROM 0 DEFAULT 0 (Only match those which are identical in all sequences.)
# PARAMETER OPTIONAL doboxes: "Display prettyboxes" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Display prettyboxes)
# PARAMETER OPTIONAL boxcol: "Colour the background in the boxes" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Colour the background in the boxes)
# PARAMETER OPTIONAL boxuse: "Colour to be used for background. (GREY\)" TYPE STRING DEFAULT GREY (Colour to be used for background. (GREY\))
# PARAMETER OPTIONAL name: "Display the sequence names" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Display the sequence names)
# PARAMETER OPTIONAL maxnamelen: "Margin size for the sequence name." TYPE INTEGER DEFAULT 10 (Margin size for the sequence name.)
# PARAMETER OPTIONAL number: "Display the residue number" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Display the residue number)
# PARAMETER OPTIONAL listoptions: "Display the date and options used" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Display the date and options used)
# PARAMETER OPTIONAL ratio: "Plurality ratio for a consensus match" TYPE DECIMAL FROM 0.0 TO 1.0 DEFAULT 0.5 (Plurality ratio for a consensus match)
# PARAMETER OPTIONAL consensus: "Display the consensus" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Display the consensus)
# PARAMETER OPTIONAL collision: "Allow collisions in calculating consensus" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Allow collisions in calculating consensus)
# PARAMETER OPTIONAL alternative: "Use alternative collisions routine" TYPE [0: "Normal collision check. (default\)", 1: "Compares identical scores with the max score found. So if any other residue matches the identical score then a collision has  occurred.", 2: "If another residue has a greater than or equal to  matching score and these do not match then a collision has   occurred.", 3: "Checks all those not in the current consensus. If any of these give a top score for matching or identical scores then a collision has occured."] FROM 1 TO 1 DEFAULT 0 (Values are 0:Normal collision check. (default\) \ 1:Compares identical scores with the max score found. So if any other residue matches the identical score then a collision has occurred. \ 2:If another residue has a greater than or equal to matching score and these do not match then a collision has occurred. \ 3:Checks all those not in the current consensus.If any of these give a top score for matching or identical scores then a collision has occured.)
# PARAMETER OPTIONAL showscore: "Print residue scores" TYPE INTEGER DEFAULT -1 (Print residue scores)
# PARAMETER OPTIONAL portrait: "Set page to Portrait" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Set page to Portrait)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

# KM 8.11. 2013
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("sequences")

#check sequece file type
inputfile.to.check <- ("sequences")
sfcheck.binary <- file.path(chipster.module.path ,"/shell/sfcheck.sh")
sfcheck.command <- paste(sfcheck.binary, emboss.path, inputfile.to.check )
str.filetype <- system(sfcheck.command, intern = TRUE )

if ( str.filetype == "Not an EMBOSS compatible sequence file"){
	stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
}

#count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount sequences -filter")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)

#round(num.queryseq)

if (num.queryseq > 1000){
	stop(paste('CHIPSTER-NOTE: Too many sequences in the alignmnet. Maximun is 1000 but your file contains ', num.queryseq ))
}

emboss.binary <- file.path(emboss.path, "prettyplot")
emboss.parameters <- paste("-auto")
emboss.parameters <- paste(emboss.parameters, "-sequences sequences")
emboss.parameters <- paste(emboss.parameters, "-graph png")
emboss.parameters <- paste(emboss.parameters, "-residuesperline", residuesperline)
emboss.parameters <- paste(emboss.parameters, "-blocksperline", blocksperline)
emboss.parameters <- paste(emboss.parameters, "-ccolours", ccolours)
emboss.parameters <- paste(emboss.parameters, "-cidentity", cidentity)
emboss.parameters <- paste(emboss.parameters, "-csimilarity", csimilarity)
emboss.parameters <- paste(emboss.parameters, "-cother", cother)
emboss.parameters <- paste(emboss.parameters, "-docolour", docolour)
emboss.parameters <- paste(emboss.parameters, "-shade", shade)
emboss.parameters <- paste(emboss.parameters, "-pair", pair)
emboss.parameters <- paste(emboss.parameters, "-identity", identity)
emboss.parameters <- paste(emboss.parameters, "-doboxes", doboxes)
emboss.parameters <- paste(emboss.parameters, "-boxcol", boxcol)
emboss.parameters <- paste(emboss.parameters, "-boxuse", boxuse)
emboss.parameters <- paste(emboss.parameters, "-name", name)
emboss.parameters <- paste(emboss.parameters, "-maxnamelen", maxnamelen)
emboss.parameters <- paste(emboss.parameters, "-number", number)
emboss.parameters <- paste(emboss.parameters, "-listoptions", listoptions)
emboss.parameters <- paste(emboss.parameters, "-ratio", ratio)
emboss.parameters <- paste(emboss.parameters, "-consensus", consensus)
emboss.parameters <- paste(emboss.parameters, "-collision", collision)
emboss.parameters <- paste(emboss.parameters, "-alternative", alternative)
emboss.parameters <- paste(emboss.parameters, "-showscore", showscore)
emboss.parameters <- paste(emboss.parameters, "-portrait", portrait)

command.full <- paste(emboss.binary, emboss.parameters, ' >> prettyplot.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> prettyplot.log' )
system(echo.command)

system(command.full)
system("ls -l >> prettyplot.log")


if ( save_log == "no") {
	system ("rm -f prettyplot.log")
}
