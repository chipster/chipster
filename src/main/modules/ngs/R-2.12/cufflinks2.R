# TOOL cufflinks2.R: "Assemble reads into transcripts using Cufflinks" (Given aligned RNA-seq reads as a BAM file, Cufflinks assembles the alignments into a parsimonious set of transcripts. It then estimates the relative abundances of these transcripts based on how many reads support each one. It is recommended to create the input BAM files using the TopHat aligner. You can merge the resulting GTF files of several samples using the Cuffmerge tool, and use the merged GTF file in differential expression analysis using Cuffdiff.)
# INPUT alignment.bam: "BAM file" TYPE BAM
# INPUT OPTIONAL annotation.gtf: "GTF for reference annotation based transcript assembly" TYPE GTF
# OUTPUT OPTIONAL genes.fpkm_tracking.tsv  
# OUTPUT OPTIONAL isoforms.fpkm_tracking.tsv
# OUTPUT OPTIONAL skipped.gtf
# OUTPUT OPTIONAL transcripts.gtf
# PARAMETER OPTIONAL normalize: "Upper-quartile normalization " TYPE [yes, no] DEFAULT no (Upper quartile normalization can improve robustness of differential expression calls for less abundant genes and transcripts. It excludes very abundant genes when normalizing expression values for the number of reads in each sample by using the upper quartile of the number of fragments mapping to individual loci.)
# PARAMETER OPTIONAL mmread: "Enable multi-mapped read correction" TYPE [yes, no] DEFAULT no (By default, Cufflinks will uniformly divide each multi-mapped read to all of the positions it maps to. If multi-mapped read correction is enabled, Cufflinks will re-estimate the transcript abundances dividing each multi-mapped read probabilistically based on the initial abundance estimation, the inferred fragment length and fragment bias, if bias correction is enabled.)
# PARAMETER OPTIONAL bias: "Bias correction for stranded data" TYPE [yes, no] DEFAULT no (Cufflinks can detect sequence-specific bias and correct for it in abundance estimation. Note that bias correction works only if your data was produced with a strand specific protocol.)
# PARAMETER OPTIONAL genome: "Genome used for bias correction" TYPE [hg19: "Human genome (hg19\)", mm9: "Mouse genome (mm9\)", mm10: "Mouse genome (mm10\)", rn4: "Rat genome (rn4\)"] DEFAULT hg19 (Genome used for bias correction.)
# PARAMETER chr: "Chromosome names in my BAM file look like" TYPE [chr1: "chr1", 1: "1"] DEFAULT chr1 (Chromosome names must match in the BAM file and in the reference annotation. Check your BAM and choose accordingly.)

# AMS 21.11.2012
# AMS 02.07.2013 Added chr1/1 option

# binary
cufflinks.binary <- c(file.path(chipster.tools.path, "cufflinks2", "cufflinks"))

# options
cufflinks.options <- ""
if (normalize == "yes") {
	cufflinks.options <- paste(cufflinks.options, "-N")
}
if (mmread == "yes") {
	cufflinks.options <- paste(cufflinks.options, "-u")
}
if (bias == "yes") {
	if (genome == "hg19"){
		genomefile <- "hg19.fa"
	}
	if (genome == "mm9"){
		genomefile <- "mm9.fa"
	}
	if (genome == "mm10"){
		genomefile <- "mm10.fa"
	}
	if (genome == "rn4"){
		genomefile <- "rn4.fa"
	}
	if (chr == "chr1"){
		genomefile <- c(file.path(chipster.tools.path, "genomes", "fasta", genomefile))
	}else{
		genomefile <- c(file.path(chipster.tools.path, "genomes", "fasta", "nochr", genomefile))
	}
	cufflinks.options <- paste(cufflinks.options, "-b", genomefile)
}
if (file.exists("annotation.gtf")){
	cufflinks.options <- paste(cufflinks.options, "-g", "annotation.gtf")
}	

# command
command <- paste(cufflinks.binary, cufflinks.options, "-q", "-o tmp", "alignment.bam")

# run
system(command)

# Rename files
if (file.exists("tmp/genes.fpkm_tracking") && file.info("tmp/genes.fpkm_tracking")$size > 0) {
	system("mv tmp/genes.fpkm_tracking genes.fpkm_tracking.tsv")
}
if (file.exists("tmp/isoforms.fpkm_tracking") && file.info("tmp/isoforms.fpkm_tracking")$size > 0) {
	system("mv tmp/isoforms.fpkm_tracking isoforms.fpkm_tracking.tsv")
}
if (file.exists("tmp/skipped.gtf") && file.info("tmp/skipped.gtf")$size > 0) {
	system("mv tmp/skipped.gtf skipped.gtf")
}
if (file.exists("tmp/transcripts.gtf") && file.info("tmp/transcripts.gtf")$size > 0) {
	system("mv tmp/transcripts.gtf transcripts.gtf")
}

