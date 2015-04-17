# TOOL select_columns.R: "Table converter" (Selects and renames table columns, and converts text files into tables.)
# INPUT input: "Query sequences" TYPE GENERIC
# OUTPUT selected.tsv
# PARAMETER cols: "Column list" TYPE STRING DEFAULT "1,2" (Select columns to be printed. Use comma separated lists of column numbers like 3,1,5 )
# PARAMETER OPTIONAL colnames: "New column names" TYPE STRING (Add a new header row for the selected columns. Use comma separated lists of column names like source,date,name. By default the first imported row is used as a header row. If you rename the columns with this tool, you probably want to start reading the data from row 2. )
# PARAMETER OPTIONAL sep: "Column separator in input file" TYPE [tab: "tabulator", space: "space or tabulator", semic: "semicolon (\;\)", doubp: "colon (\:\)", comma: "comma (\,\)", pipe: "pipe (\|\)"] DEFAULT tab (Select the column separator used to parse the input data. By default, Chipster uses tabulator.)  
# PARAMETER OPTIONAL startrow: "First row to read" TYPE INTEGER DEFAULT 1 (Skip the first few lines of the input dataset. Note that in table files, the header row is considered as the first row.)
# PARAMETER OPTIONAL skiprows: "Number of rows to remove from the end of the file" TYPE INTEGER DEFAULT 0 (Remove the given number number of lines from the end of the file or table.)
# PARAMETER OPTIONAL rstyle: "First column has rownames" TYPE [yes: Yes, no: No] DEFAULT no (Choose Yes if the first column does not have a title.)

# KM 8.11.2013
# KM 16.4.2015 Support for rownames.

# Add tabulator to the first row if this is an R-style table with rownames
if (rstyle=="yes"){
	system('printf "\t%s" "" > input.tmp')
	system('cat input >> input.tmp')
	system ('rm -f input')
	system ('mv input.tmp input')
}


# create new header if defined
if ( nchar(colnames) > 0 ){
	command.full <- paste('echo ')
	for ( i in unlist(strsplit(colnames, split=","))){
		command.full <- paste(command.full,i,'"\\t"', sep="")
	}
	command.full <- paste(command.full,' > selected.tsv1', sep="")
	#stop(paste("CHIPSTER-NOTE:", command.full))
	system(command.full)
}

# AWK command to create the datarows


# replace tabs with spaces if other separators are used
if ( sep != "tab"){
	system ('cat input | tr "\t" " " > input.tmp ')	
	system ('rm -f input')
	system ('mv input.tmp input')
}


if ( sep == "tab"){
	command.start <- paste('awk -F "\\t" \'{ if ( NR >=', startrow, ' )  print ')	
}
if ( sep == "space"){
	command.start <- paste('awk \'{if ( NR >=', startrow, ' )  print ') 
}
if ( sep == "semic"){
	command.start <- paste('awk -F ";" \'{if ( NR >=', startrow, ' )  print ') 
}
if ( sep == "doubp"){
	command.start <- paste('awk -F ":" \'{if ( NR >=', startrow, ' )  print ') 
}
if ( sep == "comma"){
	command.start <- paste('awk -F "," \'{if ( NR >=', startrow, ' )  print ') 
}
if ( sep == "pipe"){
	command.start <- paste('awk -F "|" \'{if ( NR >=', startrow, ' )  print ') 
}


command.full <- paste(command.start)

for ( i in unlist(strsplit(cols, split=","))){
	command.full <- paste(command.full,'$',i,'"\\t"', sep="")				
}

command.full <- paste(command.full,"}' input >> selected.tsv1  2>&1")
#stop(paste("CHIPSTER-NOTE:", command.full))
system(command.full)

system('sed -e s/"\t$"/""/g selected.tsv1 > selected.tsv' )

# Remove the added tabulator on the first rwo if this is a R-style table
if (rstyle=="yes"){
	system(" sed -i -e '0,/\t/s/\t//' selected.tsv ")
}


if ( skiprows > 0 ){
	num.rows.str <- system("cat selected.tsv | wc -l", intern = TRUE )
	num.rows <- as.integer(num.rows.str)
	num.rows.to.print = ( num.rows - skiprows )
	command.full <- paste('awk  \'{if ( NR <=', num.rows.to.print, ' )  print $0',"}' selected.tsv > selected.tsv2 " )
	system(command.full)
	system("rm -l selected.tsv")
	system("mv selected.tsv2 selected.tsv")
}