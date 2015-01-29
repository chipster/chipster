# TOOL edirect_fetch.R: "Retrieve sequences from NCBI" (Retrieves sequences from NCBI databases using different search criteria.)
# OUTPUT OPTIONAL sequences.txt: sequences.txt (Retrieved sequence set)
# OUTPUT OPTIONAL sequences.fasta: sequeces.fasta (Result sequences in fasta format)
# OUTPUT OPTIONAL edirect.log
# PARAMETER db: "Sequence type" TYPE [protein: Protein, nucleotide: Nucleotide] DEFAULT protein
# PARAMETER q1field: "Search field for first query term" TYPE [ALL: "All fields", TITL: Title, KYWD: Keywords, AUTH: Author, ORGN: Organism, ACCN: Accession, GENE: "Gene name", PROT: "Protein name"] DEFAULT ALL (Select the search field.)
# PARAMETER q1term: "Query term" TYPE STRING (Search term or word.)
# PARAMETER OPTIONAL log_op: "Logical operator" TYPE [AND: AND, OR: OR, NOT: NOT] DEFAULT AND (Logical operators used to combine the search terms)
# PARAMETER OPTIONAL q2field: "Search field for second query term" TYPE [ALL: "All fields", TITL: Title, KYWD: Keywords, AUTH: Author, ORGN: Organism, ACCN: Accession, GENE: "Gene name", PROT: "Protein name", SLEN: "Sequence length"] DEFAULT ALL (Select the search field )
# PARAMETER OPTIONAL q2term: "Second search term" TYPE STRING (Search term or word.)
# PARAMETER OPTIONAL q3field: "Search field for third query term" TYPE [ALL: "All fields", TITL: Title, KYWD: Keywords, AUTH: Author, ORGN: Organism, ACCN: Accession, GENE: "Gene name", PROT: "Protein name"] DEFAULT ALL (Select the search field.)
# PARAMETER OPTIONAL q3term: "Third search term" TYPE STRING (Search term or word.)
# PARAMETER outformat: "Output format" TYPE [fasta: FASTA, gp: "Genbank proteins", gb: "Genbank"] DEFAULT fasta (Logical operators used to build the search terms)


#To make edirect work:
# sudo aptitude install libhttp-parser-perl
#

outdata <- system2("hostname",args = "", stdout = TRUE, stderr = TRUE)

chipster.tools.path = "/opt/chipster/tools"
edirect.path <- file.path("/opt/chipster/tools/edirect")
esearch.path <- file.path(edirect.path, "esearch" )
efetch.path <- file.path(edirect.path, "efetch" )
xtract.path <- file.path(edirect.path, "xtract" )

outfile <- ("sequences.txt")
edoutformat <- outformat
if ( outformat == "fasta"){
	outfile <- ("sequences.fasta")
}

query.parameters <-  paste("-db", db, '-query " ')

q1field <- paste("[",q1field,"]", sep="")
q1rule <- paste(q1term, q1field, sep="")
query.parameters <- paste(query.parameters, q1rule)
q2field <- paste("[",q2field,"]", sep="")
q2rule <- paste(q2term, q2field, sep="")
if ( q2term != ""){
	query.parameters <- paste(query.parameters, log_op, q2rule)
}
q3field <- paste("[",q3field,"]", sep="")
q3rule <- paste(q3term, q3field, sep="")
if ( q3term != ""){
	query.parameters <- paste(query.parameters, log_op, q3rule) 
}  
query.parameters <- paste(query.parameters, '"')


system2( esearch.path, args = query.parameters, stdout = "query.xml" , stderr = "1search.tlog")
str.hits <- system2( xtract.path, args = " -pattern ENTREZ_DIRECT -element Count", stdin = "query.xml",  stdout = TRUE, stderr = "2xtract.tlog" )
num.hits <- as.integer(str.hits)
if (num.hits < 50000){
	if (num.hits > 0){	
		fetch.parameters <- paste( "-format", outformat)
		system2(efetch.path, args = fetch.parameters, stdout = "hits.txt",  stderr = "3fetch.tlog", stdin = "query.xml" )
		system2("grep", args = "'.' hits.txt", stdout = outfile, stderr = "4grep.tlog")
	} else {
		message.str <- paste("No hits found for query: esearch", query.parameters)
		system2("echo", args = message.str, stdout = "edirect.log")
	}
} else {
	system2("echo", args = "Query produced more than 50000 hits.", stdout = "edirect.log")
}



	

