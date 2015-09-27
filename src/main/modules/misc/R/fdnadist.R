# TOOL fdnadist.R: "Nucleic acid distance algorithm" (Calculate nucleic acid sequence matrix with PHYLIP program protdist)
# INPUT OPTIONAL sequence: sequence TYPE GENERIC 
# OUTPUT OPTIONAL dnadist.txt
# OUTPUT OPTIONAL fdnadist.log
# PARAMETER method: "Choose the method to use" TYPE [f: "F84 distance model", k: "Kimura 2-parameter distance", j: "Jukes-Cantor distance", l: "LogDet distance", s: "Similarity table"] FROM 1 TO 1 DEFAULT f (Choose the method to use)
# PARAMETER OPTIONAL gammatype: "Gamma distribution" TYPE [g: "Gamma distributed rates", i: "Gamma+invariant sites", n: "No distribution parameters used"] FROM 1 TO 1 DEFAULT n (Gamma distribution)
# PARAMETER OPTIONAL ncategories: "Number of substitution rate categories" TYPE INTEGER FROM 1 TO 9 DEFAULT 1 (Number of substitution rate categories)
# PARAMETER OPTIONAL rate: "Category rates" TYPE STRING FROM 0.0 DEFAULT 1.0 (Category rates)
# PARAMETER OPTIONAL gammacoefficient: "Coefficient of variation of substitution rate among sites" TYPE DECIMAL FROM 0.001 DEFAULT 1 (Coefficient of variation of substitution rate among sites)
# PARAMETER OPTIONAL invarfrac: "Fraction of invariant sites" TYPE DECIMAL FROM 0.0 TO 0.9999 DEFAULT 0.0 (Fraction of invariant sites)
# PARAMETER OPTIONAL ttratio: "Transition/transversion ratio" TYPE DECIMAL FROM 0.001 DEFAULT 2.0 (Transition/transversion ratio)
# PARAMETER OPTIONAL freqsfrom: "Use empirical base frequencies from seqeunce input" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Use empirical base frequencies from seqeunce input)
# PARAMETER OPTIONAL basefreq: "Base frequencies for A C G T/U (use blanks to separate\)" TYPE STRING FROM 0.0 TO 1.0 DEFAULT "0.25 0.25 0.25 0.25" (Base frequencies for A C G T/U (use blanks to separate\))
# PARAMETER OPTIONAL lower: "Output as a lower triangular distance matrix" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Output as a lower triangular distance matrix)
# PARAMETER OPTIONAL humanreadable: "Output as a human-readable distance matrix" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT <undefined> (Output as a human-readable distance matrix)
# PARAMETER OPTIONAL printdata: "Print data at start of run" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Print data at start of run)
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

emboss.binary <- file.path(emboss.path, "fdnadist")
emboss.parameters <- paste("-auto")
emboss.parameters <- paste(emboss.parameters, "-sequence", sequence)
emboss.parameters <- paste(emboss.parameters, "-outfile dnadist.txt")
emboss.parameters <- paste(emboss.parameters, "-method", method)
emboss.parameters <- paste(emboss.parameters, "-gammatype", gammatype)
emboss.parameters <- paste(emboss.parameters, "-ncategories", ncategories)
emboss.parameters <- paste(emboss.parameters, "-rate", rate)
emboss.parameters <- paste(emboss.parameters, "-gammacoefficient", gammacoefficient)
emboss.parameters <- paste(emboss.parameters, "-invarfrac", invarfrac)
emboss.parameters <- paste(emboss.parameters, "-ttratio", ttratio)
emboss.parameters <- paste(emboss.parameters, "-freqsfrom", freqsfrom)
emboss.parameters <- paste(emboss.parameters, "-basefreq", basefreq)
emboss.parameters <- paste(emboss.parameters, "-lower", lower)
emboss.parameters <- paste(emboss.parameters, "-humanreadable", humanreadable)


command.full <- paste(emboss.binary, emboss.parameters, ' >> fdnadist.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> fdnadist.log' )
system(echo.command)

system(command.full)

if ( save_log == "no") {
	system ("rm -f fdnadist.log")
}

