# TOOL fdnaml.R: "Estimate nucleotide phylogeny by maximum likelihood" (Estimate nucleotide phylogeny by maximum likelihood)
# INPUT OPTIONAL sequence: sequence TYPE GENERIC  
# OUTPUT OPTIONAL dnaml.txt
# OUTPUT OPTIONAL outtree.txt
# OUTPUT OPTIONAL fdnaml.log
# PARAMETER OPTIONAL ncategories: "Number of substitution rate categories" TYPE INTEGER FROM 1 TO 9 DEFAULT 1 (Number of substitution rate categories)
# PARAMETER OPTIONAL rate: "Rate for each category" TYPE STRING DEFAULT 1.0 (Rate for each category)
# PARAMETER OPTIONAL lengths: "Use branch lengths from user trees" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Use branch lengths from user trees)
# PARAMETER OPTIONAL ttratio: "Transition/transversion ratio" TYPE DECIMAL FROM 0.001 DEFAULT 2.0 (Transition/transversion ratio)
# PARAMETER OPTIONAL freqsfrom: "Use empirical base frequencies from seqeunce input" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Use empirical base frequencies from seqeunce input)
# PARAMETER OPTIONAL basefreq: "Base frequencies for A C G T/U (use blanks to separate\)" TYPE STRING FROM 0.0 TO 1.0 DEFAULT "0.25 0.25 0.25 0.25" (Base frequencies for A C G T/U (use blanks to separate\))
# PARAMETER OPTIONAL gammatype: "Rate variation among sites" TYPE [g: "Gamma distributed rates", i: "Gamma+invariant sites", h: "User defined HMM of rates", n: "Constant rate"] FROM 1 TO 1 DEFAULT n (Rate variation among sites)
# PARAMETER OPTIONAL gammacoefficient: "Coefficient of variation of substitution rate among sites" TYPE DECIMAL FROM 0.001 DEFAULT 1 (Coefficient of variation of substitution rate among sites)
# PARAMETER OPTIONAL ngammacat: "Number of categories (1-9\)" TYPE INTEGER FROM 1 TO 9 DEFAULT 1 (Number of categories (1-9\))
# PARAMETER OPTIONAL invarcoefficient: "Coefficient of variation of substitution rate among sites" TYPE DECIMAL FROM 0.001 DEFAULT 1 (Coefficient of variation of substitution rate among sites)
# PARAMETER OPTIONAL ninvarcat: "Number of categories (1-9\) including one for invariant sites" TYPE INTEGER FROM 1 TO 9 DEFAULT 1 (Number of categories (1-9\) including one for invariant sites)
# PARAMETER OPTIONAL invarfrac: "Fraction of invariant sites" TYPE DECIMAL FROM 0.0 TO 0.9999 DEFAULT 0.0 (Fraction of invariant sites)
# PARAMETER OPTIONAL nhmmcategories: "Number of HMM rate categories" TYPE INTEGER FROM 1 TO 9 DEFAULT 1 (Number of HMM rate categories)
# PARAMETER OPTIONAL hmmrates: "HMM category rates" TYPE STRING FROM 0.0 DEFAULT 1.0 (HMM category rates)
# PARAMETER OPTIONAL hmmprobabilities: "Probability for each HMM category" TYPE STRING FROM 0.0 TO 1.0 DEFAULT 1.0 (Probability for each HMM category)
# PARAMETER OPTIONAL adjsite: "Rates at adjacent sites correlated" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Rates at adjacent sites correlated)
# PARAMETER OPTIONAL lambda: "Mean block length of sites having the same rate" TYPE DECIMAL FROM 1.0 DEFAULT 1.0 (Mean block length of sites having the same rate)
# PARAMETER OPTIONAL njumble: "Number of times to randomise" TYPE INTEGER FROM 0 DEFAULT 0 (Number of times to randomise)
# PARAMETER OPTIONAL seed: "Random number seed between 1 and 32767 (must be odd\)" TYPE INTEGER FROM 1 TO 32767 DEFAULT 1 (Random number seed between 1 and 32767 (must be odd\))
# PARAMETER OPTIONAL global: "Global rearrangements" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Global rearrangements)
# PARAMETER OPTIONAL outgrno: "Species number to use as outgroup" TYPE INTEGER FROM 0 DEFAULT 0 (Species number to use as outgroup)
# PARAMETER OPTIONAL rough: "Speedier but rougher analysis" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Speedier but rougher analysis)
# PARAMETER OPTIONAL trout: "Write out trees to tree file" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Write out trees to tree file)
# PARAMETER OPTIONAL printdata: "Print data at start of run" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Print data at start of run)
# PARAMETER OPTIONAL treeprint: "Print out tree" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Print out tree)
# PARAMETER OPTIONAL hypstate: "Reconstruct hypothetical sequence" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Reconstruct hypothetical sequence)
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

if (num.queryseq > 100000){
	stop(paste('Too many query sequences. Maximun is 100000 but your file contains ', num.queryseq ))
}

emboss.binary <- file.path(emboss.path, "fdnaml")
emboss.parameters <- paste("-auto")
emboss.parameters <- paste(emboss.parameters, "-sequence", sequence)
emboss.parameters <- paste(emboss.parameters, "-outfile dnaml.txt")
emboss.parameters <- paste(emboss.parameters, "-outtree treefile.txt")
emboss.parameters <- paste(emboss.parameters, "-ncategories", ncategories)
emboss.parameters <- paste(emboss.parameters, '-rate "', rate, '"')
emboss.parameters <- paste(emboss.parameters, "-lengths", lengths)
emboss.parameters <- paste(emboss.parameters, "-ttratio", ttratio)
emboss.parameters <- paste(emboss.parameters, "-freqsfrom", freqsfrom)
emboss.parameters <- paste(emboss.parameters, '-basefreq "', basefreq,'"')
emboss.parameters <- paste(emboss.parameters, "-gammatype", gammatype)
emboss.parameters <- paste(emboss.parameters, "-gammacoefficient", gammacoefficient)
emboss.parameters <- paste(emboss.parameters, "-ngammacat", ngammacat)
emboss.parameters <- paste(emboss.parameters, "-invarcoefficient", invarcoefficient)
emboss.parameters <- paste(emboss.parameters, "-ninvarcat", ninvarcat)
emboss.parameters <- paste(emboss.parameters, "-invarfrac", invarfrac)
emboss.parameters <- paste(emboss.parameters, "-nhmmcategories", nhmmcategories)
emboss.parameters <- paste(emboss.parameters, "-hmmrates", hmmrates)
emboss.parameters <- paste(emboss.parameters, "-hmmprobabilities", hmmprobabilities)
emboss.parameters <- paste(emboss.parameters, "-adjsite", adjsite)
emboss.parameters <- paste(emboss.parameters, "-lambda", lambda)
emboss.parameters <- paste(emboss.parameters, "-njumble", njumble)
emboss.parameters <- paste(emboss.parameters, "-seed", seed)
emboss.parameters <- paste(emboss.parameters, "-global", global)
emboss.parameters <- paste(emboss.parameters, "-outgrno", outgrno)
emboss.parameters <- paste(emboss.parameters, "-rough", rough)
emboss.parameters <- paste(emboss.parameters, "-trout", trout)
emboss.parameters <- paste(emboss.parameters, "-printdata", printdata)
emboss.parameters <- paste(emboss.parameters, "-treeprint", treeprint)
emboss.parameters <- paste(emboss.parameters, "-hypstate", hypstate)

command.full <- paste(emboss.binary, emboss.parameters, ' >> fdnaml.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> fdnaml.log' )
system(echo.command)

system(command.full)

if ( save_log == "no") {
	system ("rm -f fdnaml.log")
}
