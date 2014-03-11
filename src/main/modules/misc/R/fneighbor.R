# TOOL fneighbor.R: "Phylogenies from distance matrix by N-J or UPGMA method" (Phylogenies from distance matrix by N-J or UPGMA method)
# INPUT OPTIONAL datafile: datafile TYPE GENERIC 
# OUTPUT OPTIONAL neighbor.txt 
# OUTPUT OPTIONAL treefile.txt 
# OUTPUT OPTIONAL fneighbor.log
# PARAMETER OPTIONAL matrixtype: "Type of data matrix" TYPE [s: Square, u: "Upper triangular", l: "Lower triangular"] FROM 1 TO 1 DEFAULT s (Type of data matrix)
# PARAMETER OPTIONAL treetype: "Neighbor-joining or UPGMA tree" TYPE [n: Neighbor-joining, u: UPGMA] FROM 1 TO 1 DEFAULT n (Neighbor-joining or UPGMA tree)
# PARAMETER OPTIONAL outgrno: "Species number to use as outgroup" TYPE INTEGER FROM 0 DEFAULT 0 (Species number to use as outgroup)
# PARAMETER OPTIONAL jumble: "Randomise input order of species" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Randomise input order of species)
# PARAMETER OPTIONAL seed: "Random number seed between 1 and 32767 (must be odd\)" TYPE INTEGER FROM 1 TO 32767 DEFAULT 1 (Random number seed between 1 and 32767 (must be odd\))
# PARAMETER OPTIONAL replicates: Subreplicates TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Subreplicates)
# PARAMETER OPTIONAL trout: "Write out trees to tree file" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Write out trees to tree file)
# PARAMETER OPTIONAL printdata: "Print data at start of run" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Print data at start of run)
# PARAMETER OPTIONAL treeprint: "Print out tree" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Print out tree)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("datafile")


emboss.binary <- file.path(emboss.path, "fneighbor")
emboss.parameters <- paste("-auto")
emboss.parameters <- paste(emboss.parameters, "-datafile datafile" )
emboss.parameters <- paste(emboss.parameters, "-outfile neighbor.txt")
emboss.parameters <- paste(emboss.parameters, "-outtreefile treefile.txt")
emboss.parameters <- paste(emboss.parameters, "-matrixtype", matrixtype)
emboss.parameters <- paste(emboss.parameters, "-treetype", treetype)
emboss.parameters <- paste(emboss.parameters, "-outgrno", outgrno)
emboss.parameters <- paste(emboss.parameters, "-jumble", jumble)
emboss.parameters <- paste(emboss.parameters, "-seed", seed)
emboss.parameters <- paste(emboss.parameters, "-replicates", replicates)
emboss.parameters <- paste(emboss.parameters, "-trout", trout)
emboss.parameters <- paste(emboss.parameters, "-printdata", printdata)
emboss.parameters <- paste(emboss.parameters, "-treeprint", treeprint)

command.full <- paste(emboss.binary, emboss.parameters, ' >> fneighbor.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> fneighbor.log' )
system(echo.command)

system(command.full)

if ( save_log == "no") {
	system ("rm -f fneighbor.log")
}
