# TOOL samtools-snp-indel-single.R: "Call SNPs and short INDELs for one diploid individual" (Call SNPs and short INDELs for one diploid individual. You can provide your own reference sequence in FASTA format or choose one of the provided reference genomes. This tool is based on the SAMtools package.)
# INPUT alignment.bam: "Sorted BAM file" TYPE GENERIC 
# INPUT OPTIONAL ownref.fa: "Reference sequence FASTA" TYPE GENERIC
# OUTPUT var.flt.vcf
# PARAMETER ref: "Reference sequence" TYPE [hg19.fa: "Human (hg19\)", mm9.fa: "Mouse (mm9\)", rn4.fa: "Rat (rn4\)"] DEFAULT hg19.fa (Reference sequence)
# PARAMETER chr: "Chromosome names in my BAM file look like" TYPE [yes: "chr1", no: "1"] DEFAULT no (Chromosome names must match in the BAM file and in the reference sequence. Check your BAM and choose accordingly. This only applies to provided reference genomes.)
# PARAMETER b: "Disable probabilistic realignment for computation of base alignment quality" TYPE [yes, no] DEFAULT no
# PARAMETER c: "Downgrading coefficient" TYPE INTEGER DEFAULT 0 (Coefficient for downgrading mapping quality for reads containing excessive mismatches. Given a read with a phred-scaled probability q of being generated from the mapped position, the new mapping quality is about sqrt( (INT-q\)\/INT\)\*INT. A zero value disables this functionality. If enabled, the recommended value for BWA is 50.)
# PARAMETER q: "Minimum mapping quality for an alignment to be used" TYPE INTEGER DEFAULT 0
# PARAMETER cq: "Minimum base quality for a base to be considered" TYPE INTEGER DEFAULT 13
# PARAMETER r: "Call variants only for a certain region" TYPE STRING DEAFAULT "all sites"
# PARAMETER d: "Output per sample read depth" TYPE [yes, no] DEFAULT no 
# PARAMETER s: "Output per sample strand bias P-value" TYPE [yes, no] DEFAULT no
# PARAMETER i: "Call only SNPs, indels not considered" TYPE [yes, no] DEFAULT no
# PARAMETER p: "Call indels only from the specified platform" TYPE STRING DEAFAULT "all"
# PARAMETER n: "Skip sites where reference is ambiguous" TYPE [yes, no] DEFAULT no (Skip sites where the REF field is not A/C/G/T)
# PARAMETER e: "Use maximum likelihood to score variants, and perform HWE analysis" TYPE [yes, no] DEFAULT no (Perform max-likelihood inference only, including estimating the site allele frequency, testing Hardy-Weinberg equlibrium and testing associations with LRT.)
# PARAMETER v: "Output variant sites only" TYPE [yes, no] DEFAULT no (Output variant sites only. Variants are called using Bayesian inference.)
# PARAMETER d: "Maximum read depth" TYPE INTEGER DEFAULT 100 (Maximum read depth. Should be adjusted to about twice the average read depth.)

# AMS 30.5.2012

# binaries
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))
bcftools.binary <- c(file.path(chipster.tools.path, "samtools", "bcftools", "bcftools"))
vcfutils.binary <- c(file.path(chipster.tools.path, "samtools", "bcftools", "vcfutils.pl"))

# path to internal reference sequences
path.refseqs <- c(file.path(chipster.tools.path, "bowtie", "indexes"))

# check which reference sequence to use: own/internal
ref.seq <- ""
input_files <- dir()
is_own <- (length(grep("ownref.fa", input_files))>0)
if (is_own) {
	ref.seq <- paste("ownref.fa")	
} else {
	if (chr == "yes"){
		ref.seq <- c(file.path(path.refseqs, ref))
	}else{
		ref.seq <- c(file.path(path.refseqs, "nochr", ref))
	}
}

# mpileup otions
mpileup.options <- paste("-uf")

if (b == "yes") {
	mpileup.options <- paste("-B")
}
mpileup.options <- paste("-C", c)
mpileup.options <- paste("-q", q)
mpileup.options <- paste("-Q", cq)
mpileup.options <- paste("-r",r)
if (d == "yes") {
	mpileup.options <- paste("-D")
}
if (s == "yes") {
	mpileup.options <- paste("-S")
}
if (i == "yes") {
	mpileup.options <- paste("-I")
}
mpileup.options <- paste("-P", p)




# bcftools options
bcftools.options <- paste("bcg")
if (n == "yes") {
	bcftools.options <- paste("-N")
}
if (e == "yes") {
	bcftools.options <- paste("-e")
}
if (v == "yes") {
	bcftools.options <- paste("-v")
}

# vcfutils otions
vcfutils.options <- paste("varFilter -D", d)


# commands
command1 <- paste(samtools.binary, "mpileup -ugf", ref.seq, "alignment.bam", "-C", c, "|", bcftools.binary, bcftools.options, "- > var.raw.bcf")
command2 <- paste(bcftools.binary, "view var.raw.bcf |", vcfutils.binary, vcfutils.options, "> var.flt.vcf")

# run
#stop(paste('CHIPSTER-NOTE: ', command1))
system(command1)
system(command2)
