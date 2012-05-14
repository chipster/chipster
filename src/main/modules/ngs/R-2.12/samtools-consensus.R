# TOOL samtools-consensus.R: "Create consensus sequence from BAM" (Create consensus sequence from a BAM file. Please note that in addition to BAM file you have to provide an index file for it. You can create the index file using the tool Index BAM. You also need to provide the reference sequence in FASTA format. This tool is based on the SAMtools package.)
# INPUT aln.bam: "BAM file to us as input" TYPE GENERIC
# INPUT aln.bam.bai: "BAM index file to us as input" TYPE GENERIC
# INPUT ref.fa: "Reference genome in FASTA format" TYPE GENERIC
# OUTPUT OPTIONAL consensus.fasta
# OUTPUT OPTIONAL output.log



# KM 26.4.2012

# samtools path
samtools.binary <- c(file.path(chipster.tools.path, "samtools","samtools"))
bcftools.binary <- c(file.path (chipster.tools.path, "samtools","bcftools","bcftools"))
vcfutils.binary <- c(file.path (chipster.tools.path, "samtools","bcftools","vcfutils.pl"))


system("echo Running command: > output.log")
samtools.command <- paste( samtools.binary, "mpileup -uf ref.fa aln.bam |",bcftools.binary ," view -cg - |",vcfutils.binary, "vcf2fq | awk 'BEGIN {a = 1}{if ($1 == \"+\") a = 0}{if (a == 1) print $0}' | sed s/\"@\"/\">\"/  > consensus.fasta")
echo.command <- paste("echo '", samtools.command, "'>> output.log")

system(echo.command)
system(samtools.command)
system("ls -l  >> output.log")


