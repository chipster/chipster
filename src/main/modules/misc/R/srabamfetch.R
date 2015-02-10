# TOOL srabamfetch.R: "Retrieve read alignments from SRA database" (Tool to retrieve alignments in BAM format from SRA databases based on the entry ID or name)
# OUTPUT OPTIONAL srafetch.log
# OUTPUT OPTIONAL sra_aln.bam
# OUTPUT OPTIONAL sra_aln.bam.bai
# PARAMETER entry_id: "SRR ID of the SRA entry" TYPE STRING DEFAULT "entry" (Give the SRR id of the SRA entry to be retrieved. For example: SRR000021) 
# PARAMETER unaligned: "Output unaligned reads" TYPE [no: "No", yes: "Yes"] DEFAULT no (Output unaligned reads along with aligned reads)
# PARAMETER primary: "Output only primary alignments" TYPE [no: "No", yes: "Yes"] DEFAULT no ( Output only primary alignments)


# KM 8.11. 2014

dump.par<-(" ")


if ( unaligned == "yes"){
	dump.par <- paste(dump.par, "--unaligned")
}

if ( primary == "yes"){
	dump.par <- paste(dump.par, "--primary")
}


   	
sra.path <- file.path(chipster.tools.path, "sratoolkit", "bin")
sra.binary <- file.path(sra.path, "sam-dump")
command.full <- paste(sra.binary, dump.par, entry_id, ' > sra_tmp.sam 2> srafetch2.tmp' )
echo.command <- paste('echo "',command.full, ' "> srafetch.log' )
system(echo.command)
system(command.full)

sam.numrows <- system("ls -s sra_tmp.sam | cut -f1 -d ' '", intern = TRUE )
	
system('cat srafetch2.tmp >> srafetch.log')

if ( sam.numrows == "0"){	
	stop(paste('CHIPSTER-NOTE: The entry ID: ', entry_id,'retuned no hits from the SRA database' ))	
}



# samtools binary
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))
	
# convert sam to bam
system(paste(samtools.binary, "view -bS sra_tmp.sam -o sra_tmp.bam"))
	
# sort bam
system(paste(samtools.binary, "sort sra_tmp.bam sra_tmp.sorted"))
	
# index bam
system(paste(samtools.binary, "index sra_tmp.sorted.bam"))
	
# rename result files
system("mv sra_tmp.sorted.bam sra_aln.bam")
system("mv sra_tmp.sorted.bam.bai sra_aln.bam.bai")
system('ls -l >> srafetch.log')
	
	
