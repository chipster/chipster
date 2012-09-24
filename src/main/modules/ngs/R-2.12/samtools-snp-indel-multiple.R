# TOOL samtools-snp-indel-multiple.R: "Call SNPs and short INDELs for multiple diploid individuals" (Call SNPs and short INDELs for multiple diploid individuals. You can provide your own reference sequence in FASTA format or choose one of the provided reference genomes. This tool is based on the SAMtools package.)
# INPUT alignment{...}.bam: "Sorted BAM file" TYPE BAM 
# INPUT OPTIONAL ownref.fa: "Reference sequence FASTA" TYPE FASTA
# OUTPUT var.flt.vcf
# PARAMETER ref: "Reference sequence" TYPE [hg19.fa: "Human (hg19\)", mm9.fa: "Mouse (mm9\)", rn4.fa: "Rat (rn4\)", e_coli.fa: "E. coli"] DEFAULT hg19.fa (Reference sequence)
# PARAMETER chr: "Chromosome names in my BAM file look like" TYPE [yes: "chr1", no: "1"] DEFAULT no (Chromosome names must match in the BAM file and in the reference sequence. Check your BAM and choose accordingly. This only applies to provided reference genomes.)
# PARAMETER p: "Call indels only from platform" TYPE [CAPILLARY, HELICOS, IONTORRENT, ILLUMINA, IONTORRENT, LS454, PACBIO, SOLID] DEFAULT ILLUMINA (Platform/technology used to produce the reads.)
# PARAMETER c: "Downgrading coefficient" TYPE INTEGER DEFAULT 0 (Coefficient for downgrading mapping quality for reads containing excessive mismatches. Given a read with a phred-scaled probability q of being generated from the mapped position, the new mapping quality is about sqrt( (INT-q\)\/INT\)\*INT. A zero value disables this functionality. If enabled, the recommended value for BWA is 50.)
# PARAMETER d: "Maximum read depth" TYPE INTEGER DEFAULT 2000 (Maximum read depth. Should be adjusted to about twice the average read depth.)

# AMS 15.6.2012

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


# commands
command1 <- paste(samtools.binary, "mpileup", "-P", p, "-ugf", ref.seq, "*.bam", "-C", c, "|", bcftools.binary, "view -bvcg - > var.raw.bcf")
command2 <- paste(bcftools.binary, "view var.raw.bcf |", vcfutils.binary, "varFilter -D", d, "> var.flt.vcf")

# run
#stop(paste('CHIPSTER-NOTE: ', command1))
system(command1)
system(command2)
