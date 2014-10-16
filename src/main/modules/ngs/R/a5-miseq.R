# TOOL a5-miseq.R: "A5 assembly pipeline for microbial genomes" (A5-miseq_ is a pipeline for assembling DNA sequence data generated on the Illumina sequencing platform. There are many situations where A5-miseq is not the right tool for the job. In order to produce accurate results, A5-miseq requires Illumina data with certain characteristics. A5-miseq will likely not work well with Illumina reads shorter than around 80nt, or reads where the base qualities are low in all or most reads before 60nt. A5-miseq assumes it is assembling homozygous haploid genomes. Use a different assembler for metagenomes and heterozygous diploid or polyploid organisms. Use a different assembler if a tool like FastQC reports your data quality is dubious. )
# INPUT reads1.fastq: "First reads file" TYPE GENERIC (First reads file)
# INPUT reads2.fastq: "Second reads file." TYPE GENERIC (Second reads file.)
# OUTPUT OPTIONAL a5_assembly.final.scaffolds.fastq: "Final a5 scaffolds." (Final a5 scaffolds.)
# OUTPUT OPTIONAL a5_assembly.final.scaffolds.fasta: "Final a5 scaffolds." (Final a5 scaffolds.)
# OUTPUT OPTIONAL a5_assembly.final.scaffolds.qvl: "Final a5 scaffolds." (Final a5 scaffolds.)
# OUTPUT OPTIONAL a5_assembly.contigs.fasta: a5_assembly.contigs.fasta (a5_assembly.contigs.fasta)
# OUTPUT OPTIONAL a5_assembly.assembly_stats.tsv: a5_assembly.assembly_stats.tsv (a5_assembly.assembly_stats.tsv)
# OUTPUT OPTIONAL a5.log
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)

# KM 18.7. 2014

options(scipen=999)
emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")

source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("reads1.fastq")
unzipIfGZipFile("reads2.fastq")



##check sequece file type
#sfcheck.binary <- file.path(chipster.module.path ,"/shell/sfcheck.sh")
#sfcheck.command <- paste(sfcheck.binary, emboss.path, "reads1.fastq" )
#str.filetype <- system(sfcheck.command, intern = TRUE )
#
#if ( str.filetype == "Not an EMBOSS compatible sequence file"){
#	stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
#}
#
#sfcheck.command <- paste(sfcheck.binary, emboss.path, "reads2.fastq" )
#str.filetype <- system(sfcheck.command, intern = TRUE )
#
#if ( str.filetype == "Not an EMBOSS compatible sequence file"){
#	stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
#}


a5.binary <- file.path("/opt/chipster/tools2/a5/a5_miseq_linux_20140604/bin/a5_pipeline.pl")
a5.parameters <- paste('reads1.fastq reads2.fastq a5_assembly ')


command.full <- paste(a5.binary, a5.parameters, ' >> a5.log 2>&1' )
echo.command <- paste('echo "',command.full, ' "> a5.log' )
system(echo.command)

system(command.full)
system ("mv a5_assembly.assembly_stats.csv a5_assembly.assembly_stats.tsv")
system ("ls -l >> a5.log ")

if ( save_log == "no") {
	system ("rm -f a5.log")
}
