# TOOL samtools-snp-indel-single.R: "Call SNPs and short INDELs for one diploid individual" (Call SNPs and short INDELs for one diploid individual. This tool is based on the SAMtools package.)
# INPUT alignment.bam TYPE GENERIC 
# OUTPUT var.flt.txt
# PARAMETER ref: "Reference sequence" TYPE [hg19.fa: "Human (hg19\)", mm9.fa: "Mouse (mm9\)", rn4.fa: "Rat (rn4\)", e_coli.fa: "E. coli"] DEFAULT hg19.fa (Reference sequence)
# PARAMETER c: "Downgrading coefficient" TYPE INTEGER DEFAULT 0 (Coefficient for downgrading mapping quality for reads containing excessive mismatches. Given a read with a phred-scaled probability q of being generated from the mapped position, the new mapping quality is about sqrt( (INT-q\)\/INT\)\*INT. A zero value disables this functionality. If enabled, the recommended value for BWA is 50.)
# PARAMETER d: "Maximum read depth" TYPE INTEGER DEFAULT 100 (Maximum read depth. Should be adjusted to about twice the average read depth.)

# binaries
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))
bcftools.binary <- c(file.path(chipster.tools.path, "samtools", "bcftools", "bcftools"))
vcfutils.binary <- c(file.path(chipster.tools.path, "samtools", "bcftools", "vcfutils.pl"))

# reference sequence
path.bowtie <- c(file.path(chipster.tools.path, "bowtie"))
path.ref.seq <- c(file.path(path.bowtie, "indexes", ref))

# commands
command1 <- paste(samtools.binary, "mpileup -ugf", path.ref.seq, "alignment.bam", "-C", c, "|", bcftools.binary, "view -bvcg - > var.raw.bcf")
command2 <- paste(bcftools.binary, "view var.raw.bcf |", vcfutils.binary, "varFilter -D", d, "> var.flt.txt")

# run
#stop(paste('CHIPSTER-NOTE: ', command1))
system(command1)
system(command2)
#system(command3)



