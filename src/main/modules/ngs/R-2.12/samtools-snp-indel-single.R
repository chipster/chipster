# TOOL samtools-snp-indel-single.R: "Call SNPs and short INDELs for one diploid individual" (Call SNPs and short INDELs for one diploid individual. You can provide your own reference sequence in FASTA format or choose one of the provided reference genomes. This tool is based on the SAMtools package.)
# INPUT alignment.bam: "Sorted BAM file" TYPE GENERIC 
# INPUT OPTIONAL ownref.fa: "Reference sequence FASTA" TYPE GENERIC
# OUTPUT var.flt.vcf
# PARAMETER ref: "Reference sequence" TYPE [hg19.fa: "Human (hg19\)", mm9.fa: "Mouse (mm9\)", rn4.fa: "Rat (rn4\)"] DEFAULT hg19.fa (Reference sequence)
# PARAMETER chr: "Chromosome names in my BAM file look like" TYPE [yes: "chr1", no: "1"] DEFAULT no (Chromosome names must match in the BAM file and in the reference sequence. Check your BAM and choose accordingly. This only applies to provided reference genomes.)
# PARAMETER mpileup.ub: "Disable probabilistic realignment for computation of BAC" TYPE [yes, no] DEFAULT no (Disable probabilistic realignment for the computation of base alignment quality (BAQ\). BAQ is the Phred-scaled probability of a read base being misaligned. Applying this option greatly helps to reduce false SNPs caused by misalignments.)
# PARAMETER mpileup.uc: "Downgrading coefficient" TYPE INTEGER DEFAULT 0 (Coefficient for downgrading mapping quality for reads containing excessive mismatches. Given a read with a phred-scaled probability q of being generated from the mapped position, the new mapping quality is about sqrt( (INT-q\)\/INT\)\*INT. A zero value disables this functionality. If enabled, the recommended value for BWA is 50.)
# PARAMETER mpileup.q: "Minimum mapping quality for an alignment to be used" TYPE INTEGER DEFAULT 0 (Minimum mapping quality for an alignment to be used. Deafult is 0)
# PARAMETER mpileup.uq: "Minimum base quality for a base to be considered" TYPE INTEGER DEFAULT 13 (Minimum base quality for a base to be considered. Default is 13)
# PARAMETER mpileup.r: "Call variants only for a certain region" TYPE STRING DEFAULT all (Only generate pileup in defined region)
# PARAMETER mpileup.ud: "Output per sample read depth" TYPE [yes, no] DEFAULT no (Output per-sample read depth)
# PARAMETER mpileup.us: "Output per sample strand bias P-value" TYPE [yes, no] DEFAULT no (Output per-sample Phred-scaled strand bias P-value)
# PARAMETER mpileup.ui: "Call only SNPs, indels not considered" TYPE [yes, no] DEFAULT no (Do not perform INDEL calling)
# PARAMETER mpileup.up: "Call indels only from platform" TYPE [all, CAPILLARY, HELICOS, IONTORRENT, ILLUMINA, IONTORRENT, LS454, PACBIO, SOLID] DEFAULT all (Platform/technology used to produce the reads.)
# PARAMETER bfctools.un: "Skip sites where reference is ambiguous" TYPE [yes, no] DEFAULT no (Skip sites where the REF field is not A/C/G/T)
# PARAMETER bfctools.e: "Use maximum likelihood to score variants, and perform HWE analysis" TYPE [yes, no] DEFAULT no (Perform max-likelihood inference only, including estimating the site allele frequency, testing Hardy-Weinberg equlibrium and testing associations with LRT.)
# PARAMETER bfctools.v: "Output variant sites only" TYPE [yes, no] DEFAULT no (Output variant sites only. Variants are called using Bayesian inference.)
# PARAMETER vcfutils.ud: "Maximum read depth" TYPE INTEGER DEFAULT 100 (Maximum read depth. Should be adjusted to about twice the average read depth.)

# AMS 30.5.2012
# AMS 17.08.2012: Additional parameters as suggested by Jarno

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
mpileup.options <- paste("mpileup -u")
if (mpileup.ub == "yes") {
	mpileup.options <- paste(mpileup.options, "-B")
}
mpileup.options <- paste(mpileup.options, "-C", mpileup.uc)
mpileup.options <- paste(mpileup.options, "-q", mpileup.q)
mpileup.options <- paste(mpileup.options, "-Q", mpileup.uq)
if (mpileup.r != "all"){
	mpileup.options <- paste(mpileup.options, "-r",mpileup.r)
}
if (mpileup.ud == "yes") {
	mpileup.options <- paste(mpileup.options, "-D")
}
if (mpileup.us == "yes") {
	mpileup.options <- paste(mpileup.options, "-S")
}
if (mpileup.ui == "yes") {
	mpileup.options <- paste(mpileup.options, "-I")
}
if (mpileup.up != "all"){
	mpileup.options <- paste(mpileup.options, "-P", mpileup.up)
}


# bcftools options
bcftools.options <- paste("view -bcg")
if (bfctools.un == "yes") {
	bcftools.options <- paste(bcftools.options, "-N")
}
if (bfctools.e == "yes") {
	bcftools.options <- paste(bcftools.options, "-e")
}
if (bfctools.v == "yes") {
	bcftools.options <- paste(bcftools.options, "-v")
}

# vcfutils otions
vcfutils.options <- paste("varFilter -D", vcfutils.ud)


# commands
command1 <- paste(samtools.binary, mpileup.options, "-f", ref.seq, "alignment.bam", "|", bcftools.binary, bcftools.options, "- > var.raw.bcf")
command2 <- paste(bcftools.binary, "view var.raw.bcf |", vcfutils.binary, vcfutils.options, "> var.flt.vcf")

# run
#stop(paste('CHIPSTER-NOTE: ', command1))
system(command1)
system(command2)
