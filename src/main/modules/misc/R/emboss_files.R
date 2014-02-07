# TOOL emboss_files.R: "Simple file operations" (Tool to remove certain characters from text files.)
# INPUT input.txt: "Input file" TYPE GENERIC
# OUTPUT output.txt
# PARAMETER OPTIONAL ecommand: "Select task" TYPE [nohtml:  "(nohtml\) Remove mark-up e.g. HTML tags from an ASCII text file", noreturn: "(noreturn\) Remove carriage return from ASCII file", nospace: "(nospace\) Remove whitespace from an ASCII text file", notab: "(notab\) Replace tabs with spaces in an ASCII text file", trimspace: "(trimspace\) Remove extra whitespace from an ASCII text file"] DEFAULT nohtml (Choose operateion to execute) 


# K.M 28.10.2013


# pb settings

if ( ecommand == "other"){
	ecommand <- ocommand
}


emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")
emboss.binary <- file.path(emboss.path,ecommand)
command.full <- paste(emboss.binary, " input.txt -filter -auto > output.txt 2>&1" )
system(command.full)
