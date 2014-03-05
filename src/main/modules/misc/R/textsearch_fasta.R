# TOOL textsearch_fasta.R: "Select sequences based on descriptions" (Select sequences based on the sequence descriptions)
# INPUT sequence: sequence TYPE GENERIC
# OUTPUT OPTIONAL matching.fasta
# OUTPUT OPTIONAL textsearch_log.txt
# PARAMETER pattern: "Enter a pattern to search for" TYPE STRING (The search pattern is a regular expression.)
# PARAMETER OPTIONAL casesensitive: "Do a case-sensitive search" TYPE [Y: Yes, N: No] DEFAULT N (Do a case-sensitive search)
# PARAMETER OPTIONAL save_log: "Collect a log file about the text search run" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the run.)

# K.M 29.10.2013


#settings
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




command.path <- file.path(chipster.module.path ,"/shell/textsearch_fasta.sh")
command.full <- paste(command.path, "-emboss_path ", emboss.path, " -sequence sequence -outfile matching.fasta -casesensitive" ,casesensitive, " -pattern ", pattern , " >> textsearch_log.txt  2>&1")
echo.command <- paste("echo", command.full )
system(echo.command)
system(command.full)

if ( save_log == "no") {
	system ("rm -f textsearch_log.txt")
}