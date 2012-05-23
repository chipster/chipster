# TOOL samtools-snp-indel-single.R: "Call SNPs and short INDELs for one diploid individual" (This tool is based on the SAMtools package.)
# INPUT alignment.bam TYPE GENERIC 
# OUTPUT var.flt.vcf
# PARAMETER ref: "Reference sequence" TYPE [hg19.fa: "Human (hg19\)"] (Reference sequence)
# PARAMETER C: "Downgrading coefficient" TYPE INTEGER DEFAULT 0 (Coefficient)
# PARAMETER D: "Maximum read depth" TYPE INTEGER DEFAULT 100 (Maximum read depth. Should be adjusted to about twice the average read depth.)
#
# AMS 23.05.2012

# binaries
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))
bcftools.binary <- c(file.path(chipster.tools.path, "samtools", "bcftools", "bcftools"))
vcfutils.binary <- c(file.path(chipster.tools.path, "samtools", "bcftools", "vcfutils.pl"))

# reference sequence
path.bowtie <- c(file.path(chipster.tools.path, "bowtie"))
path.ref.seq <- c(file.path(path.bowtie, "indexes", ref))

# commands
command1 <- paste(samtools.binary, "mpileup -ugf", path.ref.seq, "-C", C, "| bcftools view -bvcg - > var.raw.bcf"))
command2 <- paste(bcftools.binary, "view var.raw.bcf |", vcfutils.binary, "varFilter -D", D, "> var.flt.vcf"))

# run
system(command1)
system(command2)




