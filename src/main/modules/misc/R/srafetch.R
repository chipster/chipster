# TOOL srafetch.R: "Retrieve reads from SRA database" (Retrieve reads in FASTQ format from SRA databases based on the entry ID or name)
# OUTPUT OPTIONAL srafetch.log
# OUTPUT OPTIONAL sra_reads_1.fastq
# OUTPUT OPTIONAL sra_reads_2.fastq
# PARAMETER entry_id: "Name or SRR ID of the SRA entry" TYPE STRING DEFAULT "entry" (Give the SRR id of the SRA entry to be retrieved. For example: SRR000021) 
# PARAMETER dump: "Sequences to dump" TYPE [all: "All", aligned: "Only aligned sequences", unaligned: "Only unaligned sequences"  ] DEFAULT all (Define the reads to be retrieved from the SRA entry)

# KM 8.11. 2014

dump.par<-(" ")

if ( dump == "aligned"){
	dump.par <- ("--aligned")
}
if ( dump == "unaligned"){
	dump.par <- ("--unaligned")
}


sra.path <- file.path(chipster.tools.path, "sratoolkit", "bin")
sra.binary <- file.path(sra.path, "fastq-dump")
command.full <- paste(sra.binary, dump.par, '--split-files', entry_id, ' 1> srafetch.tmp 2> srafetch2.tmp' )
echo.command <- paste('echo "',command.full, ' "> srafetch.log' )
system(echo.command)
system(command.full)
system('cat srafetch.tmp >> srafetch.log')
system('cat srafetch2.tmp >> srafetch.log')
#  system('ls -l >> dbfetch.log')
system("mv *RR*_1.fastq sra_reads_1.fastq")
system("mv *RR*_2.fastq sra_reads_2.fastq")  
#  system('ls -l >> dbfetch.log')




