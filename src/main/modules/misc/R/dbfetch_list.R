# TOOL dbfetch_list.R: "Retrieve listed entries from public databases" (Tool to retrieve data from public databases based on a list of entry IDs or names)
# INPUT idlist.txt: "Names or IDs of the database entries" TYPE GENERIC
# OUTPUT OPTIONAL data.txt
# OUTPUT OPTIONAL data.fasta
# OUTPUT OPTIONAL data.html
# OUTPUT OPTIONAL log.txt
# PARAMETER db: "Database" TYPE [uniprotkb: "UniProt", embl: "EMBL", ensemblgene: "Ensembl Gene", ensemblgenomesgene: "Ensembl Genomes Gene", ensemblgenomesranscrip: "Ensembl Genomes Transcript", ensemblranscrip: "Ensembl Transcript", interpro: "InterPro", medline: "MEDLINE", pdb: "PDB", refseqn : "RefSeq nucleoide",refseqp : "RefSeq protein", taxonomy: "Taxonomy", tracearchive: "Trace Archive" ] DEFAULT uniprotkb (Database to be used.)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file.)


# KM 8.11. 2013

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

list_command <- paste("awk '{ print \"dbfetch:" ,db ,":\"$1 }' idlist.txt > idusa.txt" ,sep="" )
#stop("CHIPSTER-NOTE:",list_command)
system(list_command)
emboss.binary <- file.path(emboss.path, "textget")
command.full <- paste(emboss.binary, ' @idusa.txt -outfile data.txt > log.txt 2>&1' )
system(command.full)

system("cat idusa.txt>> log.txt")

file.type <- system("file data.txt", intern = TRUE )

if ( file.type == "data.txt: HTML document text"){
	system("mv data.txt data.html")
}
if ( save_log == "no") {
	system ("rm -f log.txt")
}	
