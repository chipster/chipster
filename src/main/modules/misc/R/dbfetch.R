# TOOL dbfetch.R: "Retrieve entries from public databases" (Tool to retrieve data from public databases based on the entry ID or name)
# OUTPUT OPTIONAL data.txt
# OUTPUT OPTIONAL data.fasta
# OUTPUT OPTIONAL data.html
# OUTPUT OPTIONAL dbfetch.log
# OUTPUT OPTIONAL sra_reads_1.fastq
# OUTPUT OPTIONAL sra_reads_2.fastq
# PARAMETER entry_id: "Name or ID of the sequence or database entry" TYPE STRING DEFAULT "entry" (Give the name of the database entry to be retrieved) 
# PARAMETER db: "Database" TYPE [sra: "SRA", uniprotkb: "UniProt", embl: "EMBL", ensemblgene: "Ensembl Gene", ensemblgenomesgene: "Ensembl Genomes Gene", ensemblgenomesranscrip: "Ensembl Genomes Transcript", ensemblranscrip: "Ensembl Transcript", interpro: "InterPro", medline: "MEDLINE", pdb: "PDB", refseqn : "RefSeq nucleoide", refseqp : "RefSeq protein", taxonomy: "Taxonomy", tracearchive: "Trace Archive" ] DEFAULT uniprotkb (Database to be used.)

# KM 8.11. 2013
if ( db == "sra"){
  sra.path <- file.path(chipster.tools.path, "sratoolkit", "bin")
  sra.binary <- file.path(sra.path, "fastq-dump")
  command.full <- paste(sra.binary, '--split-files', entry_id, ' 1> dbfetch.tmp 2> dbfetch2.tmp' )
  echo.command <- paste('echo "',command.full, ' "> dbfetch.log' )
  system(echo.command)
  system(command.full)
  system('cat dbfetch.tmp >> dbfetch.log')
  system('cat dbfetch2.tmp >> dbfetch.log')
#  system('ls -l >> dbfetch.log')
  system("mv *RR*_1.fastq sra_reads_1.fastq")
  system("mv *RR*_2.fastq sra_reads_2.fastq")  
#  system('ls -l >> dbfetch.log')
} else {  
  emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

  emboss.usa <- paste("dbfetch" ,db ,entry_id ,sep=":" )

  emboss.binary <- file.path(emboss.path, "textget")
  if ( db == "refseqp" ){
	  emboss.binary <- file.path(emboss.path, "entret")
	  emboss.usa <- paste("refseqp::dbfetch" ,db ,entry_id ,sep=":" )
  }

  command.full <- paste(emboss.binary, emboss.usa, '-outfile data.txt' )
  system(command.full)

  file.type <- system("file data.txt", intern = TRUE )

  if ( file.type == "data.txt: HTML document text"){
	  system("mv data.txt data.html")
  }
	
}