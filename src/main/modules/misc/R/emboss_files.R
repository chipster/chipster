# TOOL emboss_files.R: "Simple file operations" (Tool to remove certain characters from text files.)
# INPUT input.txt: "Input file" TYPE GENERIC
# OUTPUT output.txt
# PARAMETER OPTIONAL ecommand: "Select task" TYPE [nohtml: "(nohtml\) Remove mark-up e.g. HTML tags from a text file", noreturn: "(noreturn\) Remove carriage return from a text file", nospace: "(nospace\) Remove whitespace from a text file", notab: "(notab\) Replace tabs with spaces in a text file", trimspace: "(trimspace\) Remove extra whitespace from a text file", semicolon: "Remove semicolons (;\) from a text file", dot: "Remove dots (.\) from a text file", slash: "Remove slash characters (\/\) from a text file", at: "Remove at characters (\@\) from a text file", hash: "Remove hash characters (#\) from a text file", sqbrc: "Remove square bracket characters (\[\]\) from a text file", curbrc: "Remove curly bracket characters (\{\}\) from a text file", percent: "Remove percent characters (\%\) from a text file", dollar: "Remove dollar characters (\$\) from a text file", dq: "Remove double quotation characters (\"\) from a text file", sq: "Remove singe quotation characters (\'\) from a text file", vb: "Remove vertical bar (\|\) from a text file", et: "Remove et characters (\&\) from a text file" , gt: "Remove greater-than characters (\>\) from a text file", lt: "Remove less-than characters (\<\) from a text file"] DEFAULT nohtml (Choose operateion to execute) 


source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("input.txt")

if ( ecommand == "nohtml" || ecommand == "noreturn" ||  ecommand == "nospace" || ecommand == "notab" || ecommand == "trimspace"){
	
  emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")
  emboss.binary <- file.path(emboss.path,ecommand)
  command.full <- paste(emboss.binary, " input.txt -filter -auto > output.txt 2>&1" )
  system(command.full)

}
if ( ecommand == "semicolom"){
  system("tr -d ';' <  input.txt  > output.txt ")
}

if ( ecommand == "dot"){
	system("tr -d '.' <  input.txt  > output.txt ")
}

if ( ecommand == "slash"){
	system("tr -d '/' <  input.txt  > output.txt ")
}

if ( ecommand == "at"){
	system("tr -d '@' <  input.txt  > output.txt ")
}

if ( ecommand == "hash"){
	system("tr -d '#' <  input.txt  > output.txt ")
}

if ( ecommand == "sqbrc"){
	system("tr -d '[]' <  input.txt  > output.txt ")
}

if ( ecommand == "curbrc"){
	system("tr -d '{}' <  input.txt  > output.txt ")
}

if ( ecommand == "percent"){
	system("tr -d '%' <  input.txt  > output.txt ")
}

if ( ecommand == "dollar"){
	system("tr -d '$' <  input.txt  > output.txt ")
}

if ( ecommand == "dq"){
	system('tr -d \" <  input.txt  > output.txt ')
}

if ( ecommand == "sq"){
	system("tr -d \' <  input.txt  > output.txt ")
}

if ( ecommand == "vb"){
	system("tr -d '|' <  input.txt  > output.txt ")
}

if ( ecommand == "et"){
	system("tr -d '&' <  input.txt  > output.txt ")
}

if ( ecommand == "gt"){
	system("tr -d '>' <  input.txt  > output.txt ")
}

if ( ecommand == "lt"){
	system("tr -d '<' <  input.txt  > output.txt ")
}


