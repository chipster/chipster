# TOOL ncbi_blastn.R: "NCBI BLASTN" (Heuristic tool to search hits for a nucleotide sequence from the NCBI nucleotide sequence databases.)
# INPUT query.fa: "Query sequences" TYPE GENERIC
# OUTPUT OPTIONAL blast_results.txt
# OUTPUT OPTIONAL blast_results.xml
# OUTPUT OPTIONAL blast_results.tsv
# OUTPUT OPTIONAL blast_results.fasta
# OUTPUT OPTIONAL blast_results.csv
# OUTPUT OPTIONAL blast_results.asn1
# OUTPUT OPTIONAL blast_results.html
# OUTPUT OPTIONAL blast.log
# PARAMETER db: "Database" TYPE [nt: "NCBI non-redundant nucleoties: nt", refseq_rna: "Reference RNA sequences: refseq_rna",  refseq_genomic: "Reference genomic sequences: refseq_genomic", pdb: "PDB sequences"] DEFAULT nt (Choose a database)
# PARAMETER evalue: "Expectation value (E) threshold for saving hits" TYPE DECIMAL DEFAULT 1 (Expectation value)
# PARAMETER num_hits: "Maximun number of hits to collect per sequence" TYPE INTEGER DEFAULT 100 (Number of database sequences to show one-line descriptions for. Default 100)
# PARAMETER outfmt: "Output format type" TYPE [0: "Normal BLAST report with pairwise alignments", 1: "Query-anchored alignments showing identities", 2: "Query-anchored alignments with no identities", 3: "Flat query-anchored, show identities", 4: "flat query-anchored, no identities", 5: "XML Blast output", 6: "Tabular", 10: "Comma-separated values", 11: "BLAST archive format", 14: "Hit regions in fasta format"] DEFAULT 0 (Output format type)
# PARAMETER OPTIONAL task: "Task to execute" TYPE [blastn: "blastn", blastn-short: "blastn-short",  dc-megablast: "dc-megablast",  megablast: "megablast", rmblastn: "rmblastn"] DEFAULT megablast (Task to execute)
# PARAMETER OPTIONAL entrez_query: "Restrict search with the given Entrez query" TYPE STRING DEFAULT "None" (You can use Entrez query syntax to search a subset of the selected BLAST database. This can be helpful to limit searches to molecule types, sequence lengths or to exclude organisms.)
# PARAMETER OPTIONAL query_loc: "Location on the query sequence" TYPE STRING DEFAULT "full length" (Location of the search region on the query sequence. Format: start-stop, for example: 23-66.  Default: the whole query sequence) 
# PARAMETER OPTIONAL reward: "Reward for a nucleotide match" TYPE STRING DEFAULT "Default" (Reward for a nucleotide match)
# PARAMETER OPTIONAL penalty: "Penaltyfor a nucleotide mismatch" TYPE STRING DEFAULT "Default" (Penaly for a nucleotide mismatch)
# PARAMETER OPTIONAL gapopen: "Gap opening penalty" TYPE STRING DEFAULT "Default" (Cost to open a gap. Integer value from 6 to 25. The default value of this parameter depends on the selected scoring matrix. Note that if you assign this value, you must define also the gap extension penalty )
# PARAMETER OPTIONAL gapextend: "Gap extension penalty" TYPE STRING DEFAULT "Default" (Gap extension penalty  Integer value from 1 to 3.The default value of this parameter depends on the selected scoring matrix. Note that if you assign this value, you must define also the gap opening penalty )
# PARAMETER OPTIONAL word_size: "Word size for wordfinder algorithm" TYPE STRING DEFAULT "Default" (Word size for wordfinder algorithm)
# PARAMETER OPTIONAL dust: "Filter query sequence with DUST" TYPE [yes: Yes, no: No] DEFAULT yes (Use SEG filtering to ignore low cmoplexity regions in the query sequence) 
# PARAMETER OPTIONAL save_log: "Collect a log file about the BLAST run" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the BLAST run.)

# KM 31.10.2013

# check out if the file is compressed and if so unzip it
#source(file.path(chipster.common.path, "zip-utils.R"))
#unzipIfGZipFile("query.fa")
#unzipIfGZipFile("dbprot.fa")

# pb settings
pb.binary <- file.path(chipster.tools.path, "blast", "/ncbi-blast-2.2.28+", "bin", "pb_for_chipster")
command.start <- paste(pb.binary, "blastn")

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

#count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount -filter query.fa")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)

outfmt <- as.integer(outfmt)
#round(num.queryseq)

if (num.queryseq > 10){
	stop(paste('Too many query sequences. Maximun is 10 but your file contains ', num.queryseq ))
}

#Modify table format
outfmt.is.table <- paste("no")

if (outfmt == 6)  {
   outfmt <- paste('"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"')	
   outfmt.is.table <- paste("yes")
}
#These parameters have allways some value
general.parameters <- paste("-chipster_path /opt/chipster -remote -no_slurm -query query.fa -out blast_results -db", db )
general.parameters <- paste( general.parameters, "-task", task)
general.parameters <- paste( general.parameters, "-evalue ", evalue)
general.parameters <- paste( general.parameters, "-dust" , dust)
general.parameters <- paste( general.parameters, "-outfmt" , outfmt )

optional.parameters <- paste(" ")

if ( nchar(entrez_query) > 4 ) {
	optional.parameters <- paste(optional.parameters, '-entrez_query "', entrez_query, '"')
}


if (outfmt < 5) {
	optional.parameters <- paste(optional.parameters, " -html -num_descriptions ", num_hits, "-num_alignments", num_hits )
}

if (outfmt > 5) {
	optional.parameters <- paste(optional.parameters, " -max_target_seqs ", num_hits )
}


if ( gapopen != "Default" ) {
	if ( gapextend == "Default" ){
		stop(paste("If you define gap opening penalty you must also define gap extension pelalty"))
	}
}

if ( gapextend != "Default" ) {
	if ( gapopen == "Default" ){
		stop(paste("If you define gap opening penalty you must also define gap extension pelalty"))
	}
}

if ( gapopen !=  "Default"){
	if ( gapextend != "Default" ){
		optional.parameters <- paste(optional.parameters, " -gapopen ", gapopen, "-gapextend" , gapextend )	
	}
}

if ( word_size !=  "Default"){
	optional.parameters <- paste(optional.parameters, "-word_size", word_size )
}

if ( query_loc != "full length"){
	optional.parameters <- paste(optional.parameters, " -query_loc ", query_loc )
}

if ( reward != "Default"){
	optional.parameters <- paste(optional.parameters, " -reward ", reward )
}

if ( penalty != "Default"){
	optional.parameters <- paste(optional.parameters, " -penalty ", penalty )
}
command.end <- (paste(" >>", "blast.log 2>&1"))

#run the task
system("ls -l > blast.log")
blast.command <- paste(command.start, general.parameters, optional.parameters, command.end)
echo.command <- paste("echo '", blast.command , "' >> blast.log" )
system(echo.command)
system(blast.command)
system("echo xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
system("ls -l >> blast.log")
#rename the results
if (outfmt < 5) {
	test.command <- paste("echo renaming as html. outformat is ", outfmt, ">> blast.log")
	system (test.command)
	system ("mv blast_results blast_results.html")	
}

if (outfmt == 5) {
	system ("mv blast_results blast_results.xml")	
}

if (outfmt.is.table == "yes") {
	system ('printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "query id" "subject id" "identity percent" "alignment length" mismatches "gap opens" "query start" "query end" "subject start" "subject end" evalue "bit score" "subject description"> blast_results.tsv ')
	system ("cat blast_results >> blast_results.tsv")	
}

if (outfmt == 7) {
	system ("mv blast_results blast_results.tsv")	
}

if (outfmt == 8) {
	system ("mv blast_results blast_results.txt")	
}

if (outfmt == 10) {
	system ("mv blast_results blast_results.csv")	
}

if (outfmt == 11) {
	system ("mv blast_results blast_results.asn1")	
}

if (outfmt == 12) {
	system ("mv blast_results blast_results.txt")	
}

if (outfmt == 13) {
	system ("mv blast_results blast_results.fasta")	
}

if (outfmt == 14) {
	system ("mv blast_results blast_results.fasta")	
}

if ( save_log == "no") {
	system ("rm -f blast.log")
}
