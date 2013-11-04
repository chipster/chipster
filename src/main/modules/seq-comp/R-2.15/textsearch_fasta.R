# TOOL textsearch_fasta.R: "Search sequence descriptions" (Select sequences based on the sequence descriptions)
# INPUT sequence: sequence TYPE GENERIC
# OUTPUT matching.fasta
# OUTPUT OPTIONAL textsearch_log.txt
# PARAMETER pattern: "Enter a pattern to search for" TYPE STRING (The search pattern is a regular expression. Use a | to indicate OR. \ For example: \ human|mouse \ will find text with either 'human' OR 'mouse' in the text)
# PARAMETER OPTIONAL casesensitive: "Do a case-sensitive search" TYPE [Y: Yes, N: No] DEFAULT N (Do a case-sensitive search)
# PARAMETER OPTIONAL save_log: "Collect a log file about the BLAST run" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the BLAST run.)

# K.M 29.10.2013


#settings
emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")
command.path <- file.path(emboss.path, "textsearch_fasta.bash")
command.full <- paste(command.path, "-emboss_path ", emboss.path, " -sequence sequence -outfile matching.fasta -casesensitive" ,casesensitive, " -pattern ", pattern , " >> textsearch_log.txt  2>&1")
echo.command <- paste("echo", command.full )
system(echo.command)
system(command.full)

if ( save_log == "no") {
	system ("rm -f textsearch_log.txt")
}