# TOOL samtools-snp-indel-single.R: "Call SNPs and short INDELs" (Calls SNPs and short INDELs for one or more diploid individuals. Alignment records are grouped by sample identifiers in @RG header lines. If sample identifiers are absent, each input file is regarded as one sample. You can provide your own reference sequence in FASTA format or choose one of the provided reference genomes. This tool is based on the SAMtools package.)
# INPUT alignment{...}.bam: "Sorted BAM files" TYPE BAM 
# INPUT OPTIONAL ownref.fa: "Reference sequence FASTA" TYPE GENERIC
# OUTPUT variants.vcf
# PARAMETER organism: "Reference sequence" TYPE [other, Arabidopsis_thaliana.TAIR10.26, Bos_taurus.UMD3.1, Canis_familiaris.BROADD2.67, Canis_familiaris.CanFam3.1, Drosophila_melanogaster.BDGP5, Drosophila_melanogaster.BDGP6, Felis_catus.Felis_catus_6.2, Gallus_gallus.Galgal4, Gasterosteus_aculeatus.BROADS1, Halorubrum_lacusprofundi_atcc_49239.GCA_000022205.1.26, Homo_sapiens.GRCh37.75, Homo_sapiens.GRCh38, Homo_sapiens.NCBI36.54, Medicago_truncatula.GCA_000219495.2.26, Mus_musculus.GRCm38, Mus_musculus.NCBIM37.67, Ovis_aries.Oar_v3.1, Populus_trichocarpa.JGI2.0.26, Rattus_norvegicus.RGSC3.4.69, Rattus_norvegicus.Rnor_5.0, Schizosaccharomyces_pombe.ASM294v2.26, Sus_scrofa.Sscrofa10.2, Vitis_vinifera.IGGP_12x.26, Yersinia_enterocolitica_subsp_palearctica_y11.GCA_000253175.1.25] DEFAULT other (Reference sequence.)
# PARAMETER chr: "Chromosome names in my BAM file look like" TYPE [chr1, 1] DEFAULT 1 (Chromosome names must match in the BAM file and in the reference sequence. Check your BAM and choose accordingly. This only applies to provided reference genomes.)
# PARAMETER mpileup.r: "Call variants only for a certain region" TYPE STRING DEFAULT all (Only generate pileup in defined region. Region given as chromosome:start-end, e.g. 20:131505-131550.)
# PARAMETER OPTIONAL mpileup.ub: "Disable probabilistic realignment for computation of BAC" TYPE [yes, no] DEFAULT no (Disable probabilistic realignment for the computation of base alignment quality (BAQ\). BAQ is the Phred-scaled probability of a read base being misaligned. Applying this option greatly helps to reduce false SNPs caused by misalignments.)
# PARAMETER OPTIONAL mpileup.uc: "Downgrading coefficient" TYPE INTEGER DEFAULT 0 (Coefficient for downgrading mapping quality for reads containing excessive mismatches. The recommended value for BWA is 50. A zero value disables this functionality.)
# PARAMETER OPTIONAL mpileup.q: "Minimum mapping quality for an alignment to be used" TYPE INTEGER DEFAULT 0 (Minimum mapping quality for an alignment to be used. Default is 0.)
# PARAMETER OPTIONAL mpileup.uq: "Minimum base quality for a base to be considered" TYPE INTEGER DEFAULT 13 (Minimum base quality for a base to be considered. Default is 13.)
# PARAMETER OPTIONAL mpileup.ud: "Output per sample read depth" TYPE [yes, no] DEFAULT no (Output per-sample read depth.)
# PARAMETER OPTIONAL mpileup.us: "Output per sample strand bias P-value" TYPE [yes, no] DEFAULT no (Output per-sample Phred-scaled strand bias P-value.)
# PARAMETER OPTIONAL mpileup.ui: "Call only SNPs, indels not considered" TYPE [yes, no] DEFAULT no (Do not perform INDEL calling.)
# PARAMETER OPTIONAL bfctools.un: "Skip sites where reference is ambiguous" TYPE [yes, no] DEFAULT no (Skip sites where the REF field is not A/C/G/T.)
# PARAMETER OPTIONAL bfctools.e: "Use maximum likelihood to score variants, and perform HWE analysis" TYPE [yes, no] DEFAULT no (Perform max-likelihood inference only, including estimating the site allele frequency, testing Hardy-Weinberg equlibrium and testing associations with LRT.)
# PARAMETER OPTIONAL vcfutils.d: "Minimum read depth" TYPE INTEGER DEFAULT 2 (Minimum read depth.)
# PARAMETER OPTIONAL vcfutils.ud: "Maximum read depth" TYPE INTEGER DEFAULT 100 (Maximum read depth. Should be adjusted to about twice the average read depth.)

# AMS 30.5.2012
# AMS 17.08.2012: Additional parameters as suggested by Jarno
# AMS 24.6.2014: Changed handling of fasta file
# AMS 04.07.2014 New genome/gtf/index locations & names

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("ownref.fa")

# binaries
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))
bcftools.binary <- c(file.path(chipster.tools.path, "samtools", "bcftools", "bcftools"))
vcfutils.binary <- c(file.path(chipster.tools.path, "samtools", "bcftools", "vcfutils.pl"))

# If user provided fasta we use it, else use internal fasta
if (organism == "other"){
	# If user has provided a FASTA, we use it
	if (file.exists("ownref.fa")){
		refseq <- paste("ownref.fa")
	}else{
		stop(paste('CHIPSTER-NOTE: ', "You need to provide a FASTA file or choose one of the provided reference genomes."))
	}
}else{
	# If not, we use the internal one.
	internal.fa <- file.path(chipster.tools.path, "genomes", "fasta", paste(organism,".fa",sep="",collapse=""))
	# If chromosome names in BAM have chr, we make a temporary copy of fasta with chr names, otherwise we use it as is.
	if(chr == "chr1"){
		source(file.path(chipster.common.path, "seq-utils.R"))
		addChrToFasta(internal.fa, "internal_chr.fa") 
		refseq <- paste("internal_chr.fa")
	}else{
		refseq <- paste(internal.fa)
	}
}	

# mpileup otions
mpileup.options <- paste("mpileup -u")
if (mpileup.r != "all"){
	mpileup.options <- paste(mpileup.options, "-r", mpileup.r)
}
if (mpileup.ub == "yes") {
	mpileup.options <- paste(mpileup.options, "-B")
}
mpileup.options <- paste(mpileup.options, "-C", mpileup.uc)
mpileup.options <- paste(mpileup.options, "-q", mpileup.q)
mpileup.options <- paste(mpileup.options, "-Q", mpileup.uq)
if (mpileup.ud == "yes") {
	mpileup.options <- paste(mpileup.options, "-D")
}
if (mpileup.us == "yes") {
	mpileup.options <- paste(mpileup.options, "-S")
}
if (mpileup.ui == "yes") {
	mpileup.options <- paste(mpileup.options, "-I")
}


# bcftools options
bcftools.options <- paste("view -bcg")
if (bfctools.un == "yes") {
	bcftools.options <- paste(bcftools.options, "-N")
}
if (bfctools.e == "yes") {
	bcftools.options <- paste(bcftools.options, "-e")
}


# vcfutils otions
vcfutils.options <- paste("varFilter", "-d", vcfutils.d,"-D", vcfutils.ud)


# commands
command1 <- paste(samtools.binary, mpileup.options, "-f", refseq, "*.bam", "|", bcftools.binary, bcftools.options, "- > var.raw.bcf")
command2 <- paste(bcftools.binary, "view var.raw.bcf |", vcfutils.binary, vcfutils.options, "> variants.vcf")

# run
#stop(paste('CHIPSTER-NOTE: ', command1))
system(command1)
system(command2)
