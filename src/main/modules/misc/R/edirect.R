# TOOL edirect.R: "Retrieve sequences from NCBI" (Tool to retrieve sequences from NCBI databases using different search criteria.)
# OUTPUT OPTIONAL sequences.txt: sequences.txt (Retrieved sequence set)
# OUTPUT OPTIONAL sequences.fasta: sequeces.fasta (Result sequences in fasta format)
# OUTPUT OPTIONAL sequences.gff: sequeces.gff (Result sequences in gff format)
# OUTPUT OPTIONAL edirect.log: edirect.log (log file)
# PARAMETER db: "Sequence type" TYPE [protein: Protein, nucleotide: Nucleotide] DEFAULT protein
# PARAMETER q1field: "Serch filed for first query term" TYPE [ALL: "All fileds", TITL: Title, KYWD: Keywords, AUTH: Author, ORGN: Organism, ACCN: Accession, GENE: "Gene name", PROT: "Protein name"] DEFAULT ALL (Select the search field.)
# PARAMETER q1term: "Search term or word" TYPE STRING (Search term or word.)
# PARAMETER OPTIONAL log_op: "Logical operator" TYPE [AND: AND, OR: OR, NOT: NOT] DEFAULT AND (Logical operators used to build the serch terms)
# PARAMETER OPTIONAL q2field: "Serch filed for first query term" TYPE [ALL: "All fileds", TITL: Title, KYWD: Keywords, AUTH: Author, ORGN: Organism, ACCN: Accession, GENE: "Gene name", PROT: "Protein name"] DEFAULT ALL (Select the search field.)
# PARAMETER OPTIONAL q2term: "Search term or word" TYPE STRING (Search term or word.)
# PARAMETER OPTIONAL q3field: "Serch filed for first query term" TYPE [ALL: "All fileds", TITL: Title, KYWD: Keywords, AUTH: Author, ORGN: Organism, ACCN: Accession, GENE: "Gene name", PROT: "Protein name"] DEFAULT ALL (Select the search field.)
# PARAMETER OPTIONAL q3term: "Search term or word" TYPE STRING (Search term or word.)
# PARAMETER OPTIONAL outformat: "Output format" TYPE [fasta: FASTA, fagff: "FASTA and GFF", gp: "Genbank proteins", gb: "Genbank"] DEFAULT fasta (Logical operators used to build the search terms)
# PARAMETER OPTIONAL save_log: "Output a log file" TYPE [yes: yes, no: no] DEFAULT no (Collect a log file.)


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
edoutformat <- outformat
if ( outformat == "fasta"){
	outfile <- ("sequences.fasta")
}

# if ( outformat == "fagff"){
# 	outfile <- ("sequences.fasta")
#	 edoutformat <- ("gb")
# 	
# 	if ( db == "protein"){
#	     edoutformat <- ("gp")
#     }
# }

query.parameters <-  paste("-db", db, '-query " ')

q1field <- paste("[",q1field,"]", sep="")
q1rule <- paste(q1term, q1field, sep="")
query.parameters <- paste(query.parameters, q1rule)
q2field <- paste("[",q2field,"]", sep="")
q2rule <- paste(q2term, q2field, sep="")
query.parameters <- paste(query.parameters, q2rule)
q3field <- paste("[",q3field,"]", sep="")
q3rule <- paste(q3term, q3field, sep="")
query.parameters <- paste(query.parameters, q3rule)   
query.parameters <- paste(query.parameters, '"')

command.full <- paste(esearch.path, query.parameters, "> query.xml 2>> edirect.log" )
echo.command <- paste("echo '",command.full,"' > edirect.log" )
system(echo.command)
system(command.full)

#command.full <- paste("cat query.xml | ", xtract.path, " -pattern ENTREZ_DIRECT -element Count" )
command.full <- paste('grep "<Count>" query.xml | awk -F "[<,>]"', "'{print $3}'| head -1")
echo.command <- paste("echo '",command.full,"' >> edirect.log" )
system(echo.command)
str.hits <- system(command.full, intern = TRUE )
num.hits <- as.integer(str.hits)
echo.command <- paste("echo num.hits'", num.hits,"' >> edirect.log" )

if (num.hits > 10000){
stop(paste('CHIPSTER-NOTE: Too many hit sequences. Maximun 10000 sequeces can be retrieved, but your query would retrieve', num.hits ))
}
system("ls -l >> edirect.log")




if ( save_log == "no") {
	system ("rm -f edirect.log")
}
