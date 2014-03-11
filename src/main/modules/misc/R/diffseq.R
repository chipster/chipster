# TOOL diffseq.R: "Compare two similar sequences" (Compare and report features of two similar sequences using EMBOSS tool diffseq)
# INPUT OPTIONAL asequence: asequence TYPE GENERIC 
# INPUT OPTIONAL bsequence: bsequence TYPE GENERIC 
# OUTPUT OPTIONAL diffseq.txt
# OUTPUT OPTIONAL diffseq_feat1.txt 
# OUTPUT OPTIONAL diffseq_feat2.txt 
# PARAMETER wordsize: "Word size" TYPE INTEGER FROM 2 DEFAULT 10 (The similar regions between the two sequences are found by creating a hash table of 'wordsize'd subsequences. 10 is a reasonable default. Making this value larger (20?\) may speed up the program slightly, but will mean that any two differences within 'wordsize' of each other will be grouped as a single region of difference. This value may be made smaller (4?\) to improve the resolution of nearby differences, but the program will go much slower.)
# PARAMETER OPTIONAL globaldifferences: "Force reporting of differences at the start and end" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Normally this program will find regions of identity that are the length of the specified word-size or greater and will then report the regions of difference between these matching regions. This works well and is what most people want if they are working with long overlapping nucleic acid sequences. You are usually not interested in the non-overlapping ends of these sequences. If you have protein sequences or short RNA sequences however, you will be interested in differences at the very ends . It this option is set to be true then the differences at the ends will also be reported.)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.) 

#check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("asequence")
unzipIfGZipFile("bsequence")

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

#check sequece file type
inputfile.to.check <- ("asequence")
sfcheck.binary <- file.path(chipster.module.path ,"/shell/sfcheck.sh")
sfcheck.command <- paste(sfcheck.binary, emboss.path, inputfile.to.check )
str.filetype <- system(sfcheck.command, intern = TRUE )

if ( str.filetype == "Not an EMBOSS compatible sequence file"){
	stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
}

#check sequece file type
inputfile.to.check <- ("bsequence")
sfcheck.binary <- file.path(chipster.module.path ,"/shell/sfcheck.sh")
sfcheck.command <- paste(sfcheck.binary, emboss.path, inputfile.to.check )
str.filetype <- system(sfcheck.command, intern = TRUE )

if ( str.filetype == "Not an EMBOSS compatible sequence file"){
	stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
}

#count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount asequence -filter")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)

#round(num.queryseq)
if (num.queryseq > 1){
	stop(paste('Too many query sequences. Maximun is 1 but your file contains ', num.queryseq ))
}

#count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount bsequence -filter")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)

#round(num.queryseq)
if (num.queryseq > 1){
	stop(paste('Too many query sequences. Maximun is 1 but your file contains ', num.queryseq ))
}

emboss.binary <- file.path(emboss.path, "diffseq")
emboss.parameters <- paste("-auto")
emboss.parameters <- paste(emboss.parameters, "-asequence asequence")
emboss.parameters <- paste(emboss.parameters, "-bsequence bsequence")
emboss.parameters <- paste(emboss.parameters, "-outfile diffseq.txt")
emboss.parameters <- paste(emboss.parameters, "-aoutfeat diffseq_feat1.txt")
emboss.parameters <- paste(emboss.parameters, "-boutfeat diffseq_feat2.txt")
emboss.parameters <- paste(emboss.parameters, "-wordsize", wordsize)
emboss.parameters <- paste(emboss.parameters, "-globaldifferences", globaldifferences)


command.full <- paste(emboss.binary, emboss.parameters, ' >> diffseq.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> diffseq.log' )
system(echo.command)

system(command.full)

if ( save_log == "no") {
	system ("rm -f diffseq.log")
}
