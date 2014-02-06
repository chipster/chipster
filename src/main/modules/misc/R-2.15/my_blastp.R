# TOOL my_blastp.R: "Protein BLAST against users own protein sequence set" (Heuristic tool to search hits for a protein sequenceses from users own protein sequence set. Your query sequence set can contain 10 sequeces in maximum.)
# INPUT query.fa: "Query sequences" TYPE GENERIC
# INPUT dbprot.fa: "Database sequences" TYPE GENERIC
# OUTPUT OPTIONAL blast_results.txt
# OUTPUT OPTIONAL blast_results.xml
# OUTPUT OPTIONAL blast_results.tsv
# OUTPUT OPTIONAL blast_results.fasta
# OUTPUT OPTIONAL blast_results.csv
# OUTPUT OPTIONAL blast_results.asn1
# OUTPUT OPTIONAL blast_results.html
# OUTPUT OPTIONAL blast.log
# PARAMETER task: "BLAST program to use" TYPE [blastp: "blastp", blastp-short: "blastp short"] DEFAULT blastp (BLAST algorithm to use. Use blastp short for queries shorter than 30 residues.)
# PARAMETER OPTIONAL evalue: "Expectation threshold for saving hits" TYPE DECIMAL DEFAULT 10 (E-value specifies the statistical significance threshold for reporting matches against database sequences. The default value 10 means that 10 such matches are expected to be found merely by chance. Lower thresholds are more stringent, leading to fewer chance matches being reported.)
# PARAMETER OPTIONAL word_size: "Word size" TYPE INTEGER FROM 2 TO 10 DEFAULT 3 (The length of the seed that initiates an alignment. BLAST works by finding word-matches between the query and database sequences. One may think of this process as finding hot-spots that BLAST can then use to initiate extensions that might eventually lead to full-blown alignments. For BLASTP searches non-exact word matches are taken into account based upon the similarity between words. The amount of similarity can be varied so one normally uses just the word-sizes 2 and 3 for these searches.)
# PARAMETER OPTIONAL num_hits: "Maximun number of hits to collect per sequence" TYPE INTEGER DEFAULT 100 (Number of database sequences to show one-line descriptions for.)
# PARAMETER OPTIONAL outfmt: "Output format type" TYPE [0: "normal BLAST report with pairwise alignments", 1: "query-anchored alignments showing identities", 2: "query-anchored alignments with no identities", 3: "flat query-anchored, show identities", 4: "flat query-anchored, no identities", 5: "XML Blast output", 6: "tabular", 10: "comma-separated values", 11: "BLAST archive format", 13: "Hit sequences in fasta format", 14: "hit regions in fasta format"] DEFAULT 0 (Output format type)
# PARAMETER OPTIONAL seg: "Filter low complexity regions" TYPE [yes: Yes, no: No] DEFAULT yes (Use SEG program for filtering low complexity regions in the query sequence) 
# PARAMETER OPTIONAL entrez_query: "Entrez query to limit search" TYPE STRING DEFAULT "none" (You can use Entrez query syntax to search a subset of the selected BLAST database. This can be helpful to limit searches to molecule types, sequence lengths or to exclude organisms.)
# PARAMETER OPTIONAL query_loc: "Location on the query sequence" TYPE STRING DEFAULT "full length" (Location of the search region in the query sequence, for example: 23-66.) 
# PARAMETER OPTIONAL matrix: "Matrix" TYPE [BLOSUM45: "BLOSUM45", BLOSUM50: "BLOSUM50", BLOSUM62: "BLOSUM62", BLOSUM80: "BLOSUM80", BLOSUM90: "BLOSUM90"] DEFAULT BLOSUM62 (Weight matrix assigns a score for aligning pairs of residues, and determines overall alignment score. Experimentation has shown that the BLOSUM62 matrix is among the best for detecting most weak protein similarities. For particularly long and weak alignments, the BLOSUM45 matrix may prove superior. For proteins, shorter than 85 residues, the BLOSUM80 matrix may provide better hits"  )
# PARAMETER OPTIONAL gapopen: "Gap opening penalty" TYPE STRING DEFAULT "default" (Cost to open a gap. Integer value from 6 to 25. The default value of this parameter depends on the selected scoring matrix. Note that if you assign this value, you must define also the gap extension penalty )
# PARAMETER OPTIONAL gapextend: "Gap extension penalty" TYPE STRING DEFAULT "default" (Gap extension penalty  Integer value from 1 to 3.The default value of this parameter depends on the selected scoring matrix. Note that if you assign this value, you must define also the gap opening penalty )
# PARAMETER OPTIONAL save_log: "Output a log file" TYPE [yes: yes, no: no] DEFAULT no (Collect a log file for the BLAST run.)



# KM 24.10.2013

# check out if the file is compressed and if so unzip it
#source(file.path(chipster.common.path, "zip-utils.R"))
#unzipIfGZipFile("query.fa")
#unzipIfGZipFile("dbprot.fa")

# pb settings
pb.binary <- file.path(chipster.module.path, "/shell/pb_for_chipster.sh")
#pb.binary <- file.path(chipster.tools.path, "blast", "/ncbi-blast-2.2.29+", "bin", "pb_for_chipster")
command.start <- paste(pb.binary, "blastp")

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

#check sequece file type
inputfile.to.check <- ("query.fa")
sfcheck.binary <- file.path(chipster.module.path ,"/shell/sfcheck.sh")
sfcheck.command <- paste(sfcheck.binary, emboss.path, inputfile.to.check )
str.filetype <- system(sfcheck.command, intern = TRUE )

if ( str.filetype == "Not an EMBOSS compatible sequence file"){
	stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
}
#count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount -filter query.fa")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)

outfmt <- as.integer(outfmt)
#round(num.queryseq)

if (num.queryseq > 100000){
	stop(paste('Too many query sequences. Maximun is 100000 but your file contains ', num.queryseq ))
}


#Modify table format
outfmt.is.table <- paste("no")

if (outfmt == 6)  {
   outfmt <- paste('"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"')	
   outfmt.is.table <- paste("yes")
}
#These parameters have allways some value
general.parameters <- paste("-chipster_path /opt/chipster -no_slurm -query query.fa -out blast_results -dbprot dbprot.fa -task", task, "-evalue ", evalue, "-matrix", matrix, "-word_size" , word_size , "-seg" , seg , "-outfmt" , outfmt )
general.parameters <- paste( general.parameters, "-num_threads" , chipster.threads.max )
optional.parameters <- paste(" ")


if (outfmt < 5) {
	optional.parameters <- paste(optional.parameters, " -html -num_descriptions ", num_hits, "-num_alignments", num_hits )
}

if (outfmt > 5) {
	optional.parameters <- paste(optional.parameters, " -max_target_seqs ", num_hits )
}

#Check text formatted parameters
if (nchar(gapopen) > 2 ) {
	gapopen <- paste("default")
}
if (nchar(gapopen) < 1 ) {
	gapopen <- paste("default")
}

if (nchar(gapextend) > 2 ) {
	gapextend <- paste("default")
}
if (nchar(gapextend) < 1 ) {
	gapextend <- paste("default")
}

if (nchar(query_loc) < 3 ) {
	penalty <- paste("full length")
}

if ( gapopen != "default" ) {
	if ( gapextend == "default" ){
		stop(paste("If you define gap opening penalty you must also define gap extension pelalty"))
	}
}

if ( gapextend != "default" ) {
	if ( gapopen == "default" ){
		stop(paste("If you define gap opening penalty you must also define gap extension pelalty"))
	}
}

if ( gapopen !=  "default"){
	if ( gapextend != "default" ){
		optional.parameters <- paste(optional.parameters, " -gapopen ", gapopen, "-gapextend" , gapextend )	
	}
}

if ( query_loc != "full length"){
	optional.parameters <- paste(optional.parameters, " -query_loc ", query_loc )
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
