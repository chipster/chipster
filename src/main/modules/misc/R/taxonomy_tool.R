# TOOL taxonomy_tool.R: "Get taxonomy information" (This tool prints out taxonomy infomation based on taxonomy names or ID numbers.)
# OUTPUT OPTIONAL taxonomy.tsv
# OUTPUT OPTIONAL taxonomy.log
# OUTPUT OPTIONAL taxonomy.txt
# PARAMETER taxon: "Taxon" TYPE STRING DEFAULT 1234 (Species name or taxonomy ID number, e.g Mus musculus or 10090) 
# PARAMETER OPTIONAL tool: "Analysis type" TYPE [ taxget: "Get taxon", taxgetdown: "Get descendants of taxon", taxgetup: "Get parents of taxon", taxgetspecies: "Get all species under taxon"] DEFAULT taxget (Select the taxonomy information to be collected.)
# PARAMETER OPTIONAL oformat: "Output format type" TYPE [excel: "Table", ncbi: "NCBI formatted report", ebi: "EBI formatted report", tax: "TAX formatted report"] DEFAULT excel (Output format type)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)

# KM 8.11. 2013

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

#count the query sequeces
test.taxon <- paste("echo", taxon, ' | tr -d "1234567890" | wc -c ')
str.taxtest <- system(test.taxon, intern = TRUE )
num.tax <- as.integer(str.taxtest)

# check if taxon is a number or a string
etaxon <- paste('"taxon:',taxon ,'"', sep="" )

if ( num.tax > 1) {
 etaxon <- paste('"taxon-tax:',taxon ,'"', sep="" )
}

emboss.binary <- file.path(emboss.path, tool)
emboss.parameters <- paste("-auto")
emboss.parameters <- paste(emboss.parameters, '-taxon ', etaxon)
emboss.parameters <- paste(emboss.parameters, "-oformat", oformat)
emboss.parameters <- paste(emboss.parameters, "-outfile taxonomy.txt")
command.full <- paste(emboss.binary, emboss.parameters, ' >> taxonomy.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> taxonomy.log' )
system(echo.command)


system(command.full)
system ("ls -l >> taxonomy.log ")
if ( oformat == "excel"){
system ('curl ftp://ftp.ensemblgenomes.org/pub/current/species.txt | cut -f4 | grep -i -v "[a-z]" > ensembl_taxid')
system ('curl ftp://ftp.ensembl.org/pub/release-81/mysql/ensembl_stable_ids_81/species.txt.gz | gunzip | cut -f3 | grep -i -v "[a-z]" >> ensembl_taxid')
system ('printf "%s\t%s\t%s\t%s\t%s\t%s\n" taxid "parent taxid" level species name "available in ensembl" > taxonomy.tsv ')
system ('grep -f ensembl_taxid -w taxonomy.txt | awk \'{print $0"\tYes"}\' > taxonomy2.txt')
system ('grep -v -f ensembl_taxid -w taxonomy.txt | awk \'{print $0"\tNo"}\' >> taxonomy2.txt')
system ("cat taxonomy2.txt >> taxonomy.tsv")
system ("ls -l >> taxonomy.log ")
system ("rm -f taxonomy.txt ")
system ("rm -f taxonomy2.txt ")
}
if ( save_log == "no") {
	system ("rm -f taxonomy.log")
}