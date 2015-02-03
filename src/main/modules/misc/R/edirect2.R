# TOOL edirect2.R: "Retrieve sequences from NCBI" (Tool to retrieve sequences from NCBI databases using different search criteria.)
# OUTPUT OPTIONAL sequences.txt: sequences.txt (Retrieved sequence set)
# OUTPUT OPTIONAL sequences.fasta: sequeces.fasta (Result sequences in fasta format)
# OUTPUT OPTIONAL edirect.log
# PARAMETER db: "Sequence type" TYPE [protein: Protein, nucleotide: Nucleotide] DEFAULT protein
# PARAMETER q1field: "Serch filed for first query term" TYPE [ALL: "All fileds", TITL: Title, KYWD: Keywords, AUTH: Author, ORGN: Organism, ACCN: Accession, GENE: "Gene name", PROT: "Protein name"] DEFAULT ALL (Select the search field.)
# PARAMETER q1term: "Search term or word" TYPE STRING (Search term or word.)
# PARAMETER OPTIONAL log_op: "Logical operator" TYPE [AND: AND, OR: OR, NOT: NOT] DEFAULT AND (Logical operators used to build the serch terms)
# PARAMETER OPTIONAL q2field: "Search filed for second query term" TYPE [ALL: "All fileds", TITL: Title, KYWD: Keywords, AUTH: Author, ORGN: Organism, ACCN: Accession, GENE: "Gene name", PROT: "Protein name"] DEFAULT ALL (Select the search field.)
# PARAMETER OPTIONAL q2term: "Second search term or word" TYPE STRING (Search term or word.)
# PARAMETER OPTIONAL q3field: "Search filed for second query term" TYPE [ALL: "All fileds", TITL: Title, KYWD: Keywords, AUTH: Author, ORGN: Organism, ACCN: Accession, GENE: "Gene name", PROT: "Protein name"] DEFAULT ALL (Select the search field.)
# PARAMETER OPTIONAL q3term: "Second search term or word" TYPE STRING (Search term or word.)
# PARAMETER outformat: "Output format" TYPE [fasta: FASTA, gp: "Genbank proteins", gb: "Genbank"] DEFAULT fasta (Logical operators used to build the search terms)
# PARAMETER OPTIONAL save_log: "Output a log file" TYPE [yes: yes, no: no] DEFAULT no (Collect a log file.)

# edirect.path <- file.path(chipster.tools.path, "edirect" )
#To make edirect work:
# sudo aptitude install libhttp-parser-perl
#
chipster.tools.path = "/opt/chipster/tools"
seqret.path <- file.path(chipster.tools.path, "emboss", "bin", "seqret")
edirect.path <- file.path("/opt/chipster/tools/edirect")
esearch.path <- file.path(edirect.path, "esearch" )
efetch.path <- file.path(edirect.path, "efetch" )
xtract.path <- file.path(edirect.path, "xtract" )


outfile <- ("sequences.txt")
if ( outformat == "fasta"){
	outfile <- ("sequences.fasta")
}

query.parameters <-  paste("-db", db, '-query " ')

q1field <- paste("[",q1field,"]", sep="")
query.parameters <-  paste(query.parameters, q1term, q1field)

q2field <- paste("[",q2field,"]", sep="")
query.parameters <-  paste(query.parameters, log_op)
query.parameters <-  paste(query.parameters, q2term, q2field )

q3field <- paste("[",q3field,"]", sep="")
query.parameters <-  paste(query.parameters, log_op)
query.parameters <-  paste(query.parameters, q3term, q3field )

query.parameters <-  paste(query.parameters, '"')

#command.full <- paste(esearch.path, query.parameters, "|", efetch.path ,"-format", outformat, ">", outfile, "2>> edirect.log")

command.full <- paste(esearch.path, query.parameters, "> query.xml 2>> edirect.log" )
echo.command <- paste("echo '",command.full,"' > edirect.log" )
#echo.command <- paste("echo '",command.full,"' > command.bash" )
#system(echo.command)
system2(esearch.path, args = query.parameters, stdout = "query.xml", stderr = "edirect.log" )
system("echo step 2 >> edirect.log ")
system("ls -l >> edirect.log")

command.full <- paste("cat query.xml | ", xtract.path, " -pattern ENTREZ_DIRECT -element Count" )
echo.command <- paste("echo '",command.full,"' >> edirect.log" )
#system(echo.command)
#str.hits <- system(command.full, intern = TRUE )
#num.hits <- as.integer(str.hits)

#if (num.hits > 10000){
# stop(paste('CHIPSTER-NOTE: Too many hit sequences. Maximun 10000 sequeces can be retrieved, but your query would retrieve', num.hits ))
#}
#system("ls -l >> edirect.log")

command.full <- paste("cat query.xml | ", efetch.path, "-format", outformat, ' | grep "." >', outfile, "2>> edirect.log" )
echo.command <- paste("echo '",command.full,"' >> edirect.log" )
#system(echo.command )
system(command.full)
system("echo step 3 >> edirect.log ")
system("ls -l >> edirect.log")

if ( save_log == "no") {
	system ("rm -f edirect.log")
}



