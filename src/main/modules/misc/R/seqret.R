# TOOL seqret.R: "Sequence format conversion" (Convert sequence and alignment files into different formats with the seqret EMBOSS command)
# INPUT input.txt: "Query sequences" TYPE GENERIC
# OUTPUT OPTIONAL converted.fasta
# OUTPUT OPTIONAL converted.ncbi
# OUTPUT OPTIONAL converted.gifasta
# OUTPUT OPTIONAL converted.pir
# OUTPUT OPTIONAL converted.text
# OUTPUT OPTIONAL converted.aln
# OUTPUT OPTIONAL converted.selex
# OUTPUT OPTIONAL converted.pfam
# OUTPUT OPTIONAL converted.phylip
# OUTPUT OPTIONAL converted.hennig86
# OUTPUT OPTIONAL converted.mega
# OUTPUT OPTIONAL converted.meganon
# OUTPUT OPTIONAL converted.nexus
# OUTPUT OPTIONAL converted.fastq
# PARAMETER osformat: "Output sequence format" TYPE [ fasta: "Standard Pearson FASTA format, but with the accession number included after the identifier if available.",  ncbi: "NCBI style FASTA format with the database name, entry name and accession number separated by pipe | characters.", gifasta: "NCBI fasta format with NCBI-style IDs using GI number", pir: "NBRF/PIR entry format, as used in the PIR database sequence files.", text: "Plain sequence, no annotation or heading.", aln: "Clustalw multiple alignment format", selex: "Selex format", pfam: "Stockholm pfam format", phylip: "Phylip interleaved format", hennig86: "Hennig86 output format", mega: "Mega interleaved output format", meganon: "Mega non-interleaved output format", nexus: "Nexus/paup interleaved format", fastq: "FASTQ short read format with phred quality"] DEFAULT fasta (Output format type)

# K.M 28.10.2013

source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("input.txt")

outfile <- paste("converted.",osformat , sep="" )

emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

#check sequece file type
inputfile.to.check <- ("input.txt")
sfcheck.binary <- file.path(chipster.module.path ,"/shell/sfcheck.sh")
sfcheck.command <- paste(sfcheck.binary, emboss.path, inputfile.to.check )
str.filetype <- system(sfcheck.command, intern = TRUE )

if ( str.filetype == "Not an EMBOSS compatible sequence file"){
	stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
}

#count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount -filter input.txt")
str.queryseq <- system(seqcount.exe, intern = TRUE )
num.queryseq <- as.integer(str.queryseq)


if (num.queryseq > 50000){
	stop(paste('CHIPSTER-NOTE: Too many query sequences. Maximun is 50000 but your file contains ', num.queryseq ))
}

# command settings
ecommand <- paste("seqret")
emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")
emboss.binary <- file.path(emboss.path,ecommand)
command.full <- paste(emboss.binary, "-sequence input.txt -osformat", osformat,"-auto -outseq ", outfile, " 2> log.txt" )
system(command.full)
system("ls -l > log.txt")


