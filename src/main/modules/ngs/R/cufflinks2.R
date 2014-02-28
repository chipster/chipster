# TOOL cufflinks2.R: "Assemble reads into transcripts using Cufflinks" (Given aligned RNA-seq reads as a BAM file, Cufflinks assembles the alignments into a parsimonious set of transcripts. It then estimates the relative abundances of these transcripts based on how many reads support each one. It is recommended to create the input BAM files using the TopHat aligner. You can merge the resulting GTF files of several samples using the Cuffmerge tool, and use the merged GTF file in differential expression analysis using Cuffdiff.)
# INPUT alignment.bam: "BAM file" TYPE BAM
# INPUT OPTIONAL annotation.gtf: "GTF for reference annotation based transcript assembly" TYPE GTF
# OUTPUT OPTIONAL genes.fpkm_tracking.tsv  
# OUTPUT OPTIONAL isoforms.fpkm_tracking.tsv
# OUTPUT OPTIONAL skipped.gtf
# OUTPUT OPTIONAL transcripts.gtf
# PARAMETER chr: "Chromosome names in my BAM file look like" TYPE [chr1: "chr1", 1: "1"] DEFAULT chr1 (Chromosome names must match in the BAM file and in the reference annotation. Check your BAM and choose accordingly.)
# PARAMETER OPTIONAL lg: "Estimate expression of known isoforms and do not assembl novel transcripts" TYPE [yes, no] DEFAULT no (Reference annotation (a GFF/GTF file\) is used to estimate isoform expression. Program will not assemble novel transcripts, and the it will ignore alignments not structurally compatible with any reference transcript. You can supply your own GTF file or use one of the provided annotations.)
# PARAMETER OPTIONAL ug: "Do reference annotation based transcript assembly" TYPE [yes, no] DEFAULT no (Cufflinks will use the supplied reference annotation (a GFF/GTF file\) to guide RABT assembly. Reference transcripts will be tiled with faux-reads to provide additional information in assembly. Output will include all reference transcripts as well as any novel genes and isoforms that are assembled. You can supply your own GTF file or use one of the provided annotations.)
# PARAMETER OPTIONAL internalgtf: "Annotation GTF" TYPE [hg19: "Human (hg19\)", mm9: "Mouse (mm9\)", mm10: "Mouse (mm10\)", rn4: "Rat (rn4\)"] DEFAULT hg19 (You can use own GTF file or one of those provided on the server.)
# PARAMETER OPTIONAL mmread: "Enable multi-mapped read correction" TYPE [yes, no] DEFAULT no (By default, Cufflinks will uniformly divide each multi-mapped read to all of the positions it maps to. If multi-mapped read correction is enabled, Cufflinks will re-estimate the transcript abundances dividing each multi-mapped read probabilistically based on the initial abundance estimation, the inferred fragment length and fragment bias, if bias correction is enabled.)
# PARAMETER OPTIONAL bias: "Correct for sequence-specific bias" TYPE [yes, no] DEFAULT no (Cufflinks can detect sequence-specific bias and correct for it in abundance estimation.)
# PARAMETER OPTIONAL genome: "Genome used for bias correction" TYPE [hg19: "Human genome (hg19\)", mm9: "Mouse genome (mm9\)", mm10: "Mouse genome (mm10\)", rn4: "Rat genome (rn4\)"] DEFAULT hg19 (Genome used for bias correction.)

# AMS 21.11.2012
# AMS 2.07.2013 Added chr1/1 option
# EK 3.11.2013 Renamed bias correction parameter, removed quartile normalization parameter
# AMS 11.11.2013 Added thread support
# AMS 24.2.2014 Add G option, change g option, add support for internal gtf:s for both

# binary
cufflinks.binary <- c(file.path(chipster.tools.path, "cufflinks2", "cufflinks"))

# options
cufflinks.options <- ""
cufflinks.options <- paste(cufflinks.options, "-p", chipster.threads.max)
# if (normalize == "yes") {
# 	cufflinks.options <- paste(cufflinks.options, "-N")
# }
if (file.exists("annotation.gtf")){
	annotation.file <- "annotation.gtf"
}else{
	if (internalgtf == "hg19") {
		if (chr == 1){
			annotation.file <- "Homo_sapiens.GRCh37.68.gtf"
		}else {
			annotation.file <- "Homo_sapiens.GRCh37.68.chr.gtf"
		}		
	}
	if (internalgtf == "mm9") {
		if (chr == 1){
			annotation.file <- "Mus_musculus.NCBIM37.62.gtf"
		}else {
			annotation.file <- "Mus_musculus.NCBIM37.62.chr.gtf"
		}
	}
	if (internalgtf == "mm10") {
		if (chr == 1){
			annotation.file <- "Mus_musculus.GRCm38.68.gtf"
		}else{
			annotation.file <- "Mus_musculus.GRCm38.68.chr.gtf"
		}
	}
	if (internalgtf == "rn4") {
		if (chr == 1){
			annotation.file <- "Rattus_norvegicus.RGSC3.4.68.gtf"
		}else{
			annotation.file <- "Rattus_norvegicus.RGSC3.4.68.chr.gtf"
		}
	}
	annotation.file <- c(file.path(chipster.tools.path, "genomes", "gtf", annotation.file))
}
if (lg == "yes"){
	cufflinks.options <-paste(cufflinks.options, "-g", annotation.file)
}
if (ug == "yes"){
	cufflinks.options <-paste(cufflinks.options, "-G", annotation.file)
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

