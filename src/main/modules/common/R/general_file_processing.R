# TOOL general_file_processing.R: "Modify text" (This tool can be used to modify txt, tsv, BED and GTF files. For example, you can replace text, or extract rows containing a given text. Note that the search string is interpreted as a regular expression were .,* and + characters have special functions)
# INPUT input: "Input file" TYPE GENERIC
# OUTPUT OPTIONAL selected.tsv
# OUTPUT OPTIONAL selected.txt
# OUTPUT OPTIONAL selected.bed
# OUTPUT OPTIONAL selected.gtf
# OUTPUT OPTIONAL file_operation.log
# PARAMETER operation: "Operation" TYPE [select: "Select rows with a regular expression", exclude: "Exclude rows with a regular expression", replace: "Replace text", pick_rows: "Select a set of rows from the file" ] DEFAULT replace (Operation to be performed for the selected text or table file)
# PARAMETER OPTIONAL sstring: "Search string" TYPE STRING (Search expression)
# PARAMETER OPTIONAL rstring: "Replacement string" TYPE STRING (Replacement string)
# PARAMETER OPTIONAL startrow: "First row to select" TYPE INTEGER DEFAULT 1 (Number of the first row to be selected. Note that in table files, the header row is considered as the first row.)
# PARAMETER OPTIONAL stoprow: "Last row to select" TYPE INTEGER DEFAULT 10000000 (Number of the last row to be selected.)
# PARAMETER OPTIONAL fstyle: "Input file format" TYPE [txt: Text, tsv: Table, bed: BED, gtf: GTF] DEFAULT txt (Is the input file a text file, tab-delimited table, BED file or GTF file)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)

# KM 10.4.2015
# AMS 9.11.2015 Added support for compressed input files

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("input")

if ( nchar(sstring)>50 ){
stop(paste("CHIPSTER-NOTE:", "Too long search string"))
}

if ( nchar(rstring)>50 ){
	stop(paste("CHIPSTER-NOTE:", "Too long replacement string"))
}

if (operation=="select"){
	command <- paste("grep '",sstring, "' input > output.tmp", sep="")
}

if (operation=="exclude"){
	command <- paste("grep -v '",sstring, "' input > output.tmp", sep="")
}
if (operation=="replace"){
	command <- paste('sed -e s/"',sstring, '"/"', rstring, '"/g input > output.tmp', sep="")
}

if (operation=="delete"){
	command <- paste('tr -d s/"',sstring,'" < input > output.tmp', sep="")
}
if (operation=="pick_rows"){
	command <- paste("awk '{if ( NR >= ", startrow," ) if ( NR <= ",stoprow, ") print $0}' input > output.tmp", sep="")
}

echo.command <- paste('echo "', command, '" > file_operation.log')
system(echo.command)
command.full <- paste(command,' 2>>file_operation.log ' )
system(command.full)



if ( save_log == "no"){
	system ("rm -f file_operation.log")
}else{
	system ('echo "Number of rows in the input file:" >> file_operation.log  ' )
	system ("wc -l input >> file_operation.log" )
	system ('echo "Number of rows in the output file:" >> file_operation.log  ' )
	system ("wc -l output.tmp >> file_operation.log" )
}


if (fstyle=="txt"){
	system('mv output.tmp selected.txt')
}
if (fstyle=="tsv"){
	system('mv output.tmp selected.tsv')
}
if (fstyle=="bed"){
	system('mv output.tmp selected.bed')
}
if (fstyle=="gtf"){
	system('mv output.tmp selected.gtf')
}
