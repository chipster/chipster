# TOOL fprotdist.R: "Protein distance algorithm" (Calculate protein distance matrix with PHYLIP program protdist)
# INPUT OPTIONAL sequence: sequence TYPE GENERIC 
# OUTPUT OPTIONAL protdist.txt
# OUTPUT OPTIONAL fprotdist.log
# PARAMETER OPTIONAL ncategories: "Number of substitution rate categories" TYPE INTEGER FROM 1 TO 9 DEFAULT 1 (Number of substitution rate categories)
# PARAMETER OPTIONAL rate: "Rate for each category" TYPE STRING DEFAULT 1.0 (Rate for each category)
# PARAMETER OPTIONAL method: "Choose the method to use" TYPE [j: "Jones-Taylor-Thornton matrix", h: "Henikoff/Tiller PMB matrix", d: "Dayhoff PAM matrix", k: "Kimura formula", s: "Similarity table", c: "Categories model"] FROM 1 TO 1 DEFAULT j (Choose the method to use)
# PARAMETER OPTIONAL gammatype: "Rate variation among sites" TYPE [g: "Gamma distributed rates", i: "Gamma+invariant sites", c: "Constant rate"] FROM 1 TO 1 DEFAULT c (Rate variation among sites)
# PARAMETER OPTIONAL gammacoefficient: "Coefficient of variation of substitution rate among sites" TYPE DECIMAL FROM 0.001 DEFAULT 1 (Coefficient of variation of substitution rate among sites)
# PARAMETER OPTIONAL invarcoefficient: "Coefficient of variation of substitution rate among sites" TYPE DECIMAL FROM 0.001 DEFAULT 1 (Coefficient of variation of substitution rate among sites)
# PARAMETER OPTIONAL aacateg: "Choose the category to use" TYPE [G: "George/Hunt/Barker (Cys\), (Met Val Leu Ileu\), (Gly Ala Ser Thr Pro\)", C: "Chemical (Cys Met\), (Val Leu Ileu Gly Ala Ser Thr\), (Pro\)", H: "Hall (Cys\), (Met Val Leu Ileu\), (Gly Ala Ser Thr\),(Pro\)"] FROM 1 TO 1 DEFAULT G (Choose the category to use)
# PARAMETER OPTIONAL whichcode: "Which genetic code" TYPE [u: Universal, c: Ciliate, m: "Universal mitochondrial", v: "Vertebrate mitochondrial", f: "Fly mitochondrial", y: "Yeast mitochondrial"] FROM 1 TO 1 DEFAULT u (Which genetic code)
# PARAMETER OPTIONAL ease: "Prob change category (1.0=easy\)" TYPE DECIMAL FROM 0.0 TO 1.0 DEFAULT 0.457 (Prob change category (1.0=easy\))
# PARAMETER OPTIONAL ttratio: "Transition/transversion ratio" TYPE DECIMAL FROM 0.0 DEFAULT 2.0 (Transition/transversion ratio)
# PARAMETER OPTIONAL basefreq: "Base frequencies for A C G T/U (use blanks to separate\)" TYPE STRING FROM 0.0 TO 1.0 DEFAULT "0.25 0.25 0.25 0.25" (Base frequencies for A C G T/U (use blanks to separate\))
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

emboss.binary <- file.path(emboss.path, "fprotdist")
emboss.parameters <- paste("-auto")
emboss.parameters <- paste(emboss.parameters, "-sequence sequence")
emboss.parameters <- paste(emboss.parameters, "-outfile protdist.txt")
emboss.parameters <- paste(emboss.parameters, "-ncategories", ncategories)
emboss.parameters <- paste(emboss.parameters, "-rate", rate)
emboss.parameters <- paste(emboss.parameters, "-method", method)
emboss.parameters <- paste(emboss.parameters, "-gammatype", gammatype)
emboss.parameters <- paste(emboss.parameters, "-gammacoefficient", gammacoefficient)
emboss.parameters <- paste(emboss.parameters, "-invarcoefficient", invarcoefficient)
emboss.parameters <- paste(emboss.parameters, "-aacateg", aacateg)
emboss.parameters <- paste(emboss.parameters, "-whichcode", whichcode)
emboss.parameters <- paste(emboss.parameters, "-ease", ease)
emboss.parameters <- paste(emboss.parameters, "-ttratio", ttratio)
emboss.parameters <- paste(emboss.parameters, '-basefreq "', basefreq,'"')
emboss.parameters <- paste(emboss.parameters, "-printdata", printdata)

command.full <- paste(emboss.binary, emboss.parameters, ' >> fprotdist.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> fprotdist.log' )
system(echo.command)

system(command.full)

if ( save_log == "no") {
	system ("rm -f fprotdist.log")
}
