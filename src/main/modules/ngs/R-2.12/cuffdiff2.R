# TOOL cuffdiff2.R: "Differential expression using Cuffdiff" (Given GTF and BAM files, Cuffdiff performs differential expression analysis of genes and transcripts using the Cufflinks algorithm. It is currently not possible to use replicate samples in Chipster, so you have to merge all samples belonging to the same experiment group into a single BAM file.)
# INPUT treatment1.bam: "Treatment BAM" TYPE BAM
# INPUT control1.bam: "Control BAM" TYPE BAM
# INPUT OPTIONAL annotation.gtf: "Annotation GTF" TYPE GTF
# OUTPUT OPTIONAL de-genes-cufflinks.tsv
# OUTPUT OPTIONAL de-isoforms-cufflinks.tsv
# OUTPUT OPTIONAL cufflinks-log.txt
# OUTPUT OPTIONAL de-genes-cufflinks.bed
# OUTPUT OPTIONAL de-isoforms-cufflinks.bed
# OUTPUT OPTIONAL cds.count_tracking.tsv
# OUTPUT OPTIONAL cds.diff.tsv
# OUTPUT OPTIONAL cds.fpkm_tracking.tsv
# OUTPUT OPTIONAL cds.read_group_tracking.tsv
# OUTPUT OPTIONAL cds_exp.diff.tsv
# OUTPUT OPTIONAL gene_exp.diff.tsv
# OUTPUT OPTIONAL genes.count_tracking.tsv
# OUTPUT OPTIONAL genes.fpkm_tracking.tsv
# OUTPUT OPTIONAL genes.read_group_tracking.tsv
# OUTPUT OPTIONAL isoform_exp.diff.tsv
# OUTPUT OPTIONAL isoforms.count_tracking.tsv
# OUTPUT OPTIONAL isoforms.fpkm_tracking.tsv
# OUTPUT OPTIONAL isoforms.read_group_tracking.tsv
# OUTPUT OPTIONAL promoters.diff.tsv
# OUTPUT OPTIONAL read_groups.info.txt
# OUTPUT OPTIONAL run.info.txt
# OUTPUT OPTIONAL splicing.diff.tsv
# OUTPUT OPTIONAL tss_group_exp.diff.tsv
# OUTPUT OPTIONAL tss_groups.count_tracking.tsv
# OUTPUT OPTIONAL tss_groups.fpkm_tracking.tsv
# OUTPUT OPTIONAL tss_groups.read_group_tracking.tsv
# PARAMETER output.type: "Output type" TYPE [concise, complete] DEFAULT concise (Cuffdiff produces a large number of output files (over 20\). You can choose to see the complete output or just concise processed output.)
# PARAMETER internalgtf: "Annotation GTF" TYPE [hg19: "Human (hg19\)", mm9: "Mouse (mm9\)", mm10: "Mouse (mm10\)", rn4: "Rat (rn4\)"] DEFAULT hg19 (You can use own GTF file or one of those provided on the server.)
# PARAMETER chr: "Chromosome names in my BAM file look like" TYPE [chr1: "chr1", 1: "1"] DEFAULT chr1 (Chromosome names must match in the BAM file and in the reference annotation. Check your BAM and choose accordingly.)
# PARAMETER OPTIONAL fdr: "Allowed false discovery rate" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (FDR-adjusted p-values (q-values\) are calculated. The concise output files include only those genes or transcripts which have a q-value lower than the given FDR. The value of the Significant-column is adjusted accordingly (yes/no\) in all output files.) 
# PARAMETER OPTIONAL mmread: "Enable multi-mapped read correction" TYPE [yes, no] DEFAULT no (By default, Cufflinks will uniformly divide each multi-mapped read to all of the positions it maps to. If multi-mapped read correction is enabled, Cufflinks will re-estimate the transcript abundances dividing each multi-mapped read probabilistically based on the initial abundance estimation, the inferred fragment length and fragment bias, if bias correction is enabled.)
# PARAMETER OPTIONAL bias: "Bias correction for stranded data" TYPE [yes, no] DEFAULT no (Cuffdiff can detect sequence-specific bias and correct for it in abundance estimation. Note that bias correction works only if your data was produced with a strand specific protocol.)
# PARAMETER OPTIONAL genome: "Genome used for bias correction" TYPE [hg19: "Human genome (hg19\)", mm9: "Mouse genome (mm9\)", mm10: "Mouse genome (mm10\)", rn4: "Rat genome (rn4\)"] DEFAULT hg19 (Genome used for bias correction.)

# AMS 21.1.2013
# check column renaming when replicates are enabled
# This parameter is not yet functional due to a bug in Cuffdiff itself. PARAMETER OPTIONAL normalize: "Upper-quartile normalization " TYPE [yes, no] DEFAULT yes (Upper quartile normalization can improve robustness of differential expression calls for less abundant genes and transcripts. It excludes very abundant genes when normalizing expression values for the number of reads in each sample by using the upper quartile of the number of fragments mapping to individual loci.)
# AMS 02.07.2013 Added chr1/1 option, fixed handling of errors when no results are found.

# binary
cuffdiff.binary <- c(file.path(chipster.tools.path, "cufflinks2", "cuffdiff"))

# options
cuffdiff.options <- ""
cuffdiff.options <- paste(cuffdiff.options, "-FDR", fdr)
# if (normalize == "yes") {
#	cuffdiff.options <- paste(cuffdiff.options, "-N")
# }
if (mmread == "yes") {
	cuffdiff.options <- paste(cuffdiff.options, "-u")
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
	cuffdiff.options <- paste(cuffdiff.options, "-b", genomefile)
}
if (file.exists("annotation.gtf")){
	cuffdiff.options <- paste(cuffdiff.options, "annotation.gtf")
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
	cuffdiff.options <- paste(cuffdiff.options, annotation.file)
}

# command
command <- paste(cuffdiff.binary, "-q", "-o tmp", cuffdiff.options, "treatment1.bam", "control1.bam")

# run
#stop(paste('CHIPSTER-NOTE: ', command))
system(command)

# The following code copied from dea-cufflinks.R
#
source(file.path(chipster.common.path, "bed-utils.R")) # bed sort

# Only do post-processing if file exists
if (file.exists("tmp/gene_exp.diff") && file.info("tmp/gene_exp.diff")$size > 0){
#if (file.exists("tmp/gene_exp.diff")){
	# DE genes
	# Extract chromosome locations and add in the first three table columns
	dat <- read.table(file="tmp/gene_exp.diff", header=T, sep="\t")
	regions_list <- as.character(dat$locus)
	chr_list <- character(length(regions_list))
	start_list <- numeric(length(regions_list))
	end_list <- numeric(length(regions_list))
	for (count in 1:length(regions_list)) {
		chr_list[count] <- unlist(strsplit (regions_list[count], split=":")) [1]
		start_list[count] <- unlist(strsplit(unlist(strsplit(regions_list[count], split=":"))[2], split="-")) [1]
		end_list[count] <- unlist(strsplit(unlist(strsplit(regions_list[count], split=":"))[2], split="-")) [2]	
	}
	dat2 <- data.frame(chr=chr_list, start=start_list, end=end_list, dat)

	# Rename gene to symbol for compatibility with venn diagram
	# colnames (dat2) [5] <- "ID"
	colnames (dat2) [6] <- "symbol"
	colnames (dat2) [11] <- "FPKM_1"
	colnames (dat2) [12] <- "FPKM_2"
	colnames (dat2) [13] <- "log2_FC"

	# Filter the gene output based on status and significant columns
	dat2 <- dat2[dat2$status=="OK",]
	results_list <- dat2
	results_list <- results_list[results_list$significant=="yes",]
	
	# Write output only if signicant results exist
	if (dim(results_list)[1] > 0) {
	
		# order according to increasing q-value
		results_list <- results_list[order(results_list$q_value, decreasing=FALSE),]
		number_genes <- dim (results_list) [1]
		row_names <- 1:number_genes
		rownames(results_list) <- row_names
		write.table(results_list, file="de-genes-cufflinks.tsv", sep="\t", row.names=TRUE, col.names=T, quote=F)
	
		# Also output a BED file for visualization and region matching tools
		bed_output <- results_list[,c("chr","start","end","gene_id","log2_FC")]
		# sort according to chromosome location
		write.table(bed_output, file="sortme.bed", sep="\t", row.names=F, col.names=T, quote=F)
		bed <- read.table(file="sortme.bed", skip=1, sep="\t") # assume file has 1 line header
		colnames(bed)[1:2] <- c("chr", "start")  # these named columns are required for sorting 
		sorted.bed <- sort.bed(bed)
		write.table(sorted.bed, file="de-genes-cufflinks.bed", sep="\t", row.names=F, col.names=F, quote=F)
	}	

	# Report numbers to the log file
	number_genes_tested <- dim(dat)[1]
	sink(file="cufflinks-log.txt")
	if (dim(results_list)[1] > 0) {
		number_filtered <- number_genes_tested-dim(results_list)[1]
		number_significant <- dim(results_list)[1]
		cat("GENE TEST SUMMARY\n")
		cat("In total,", number_genes_tested, "genes were tested for differential expression.\n")
		cat("Of these,", number_filtered, "didn't fulfill the technical criteria for testing or the significance cut-off specified.\n")
		cat(number_significant, "genes were found to be statistically significantly differentially expressed.\n")	
	} else {
		cat("GENE TEST SUMMARY\n")
		cat("Out of the", number_genes_tested, "genes tested, there were no statistically significantly differentially expressed ones found.")
	}
	sink()
}else{
	write("Cuffdiff did not complete succesfully. Please check your input files.", "cufflinks-log.txt", append=F)
}

# Only do post-processing if file exists
if (file.exists("tmp/isoform_exp.diff") && file.info("tmp/isoform_exp.diff")$size > 0) {
	# DE isoforms
	# Extract chromosome locations and add in the first three table columns
	dat <- read.table(file="tmp/isoform_exp.diff", header=T, sep="\t")
	regions_list <- as.character(dat$locus)
	chr_list <- character(length(regions_list))
	start_list <- numeric(length(regions_list))
	end_list <- numeric(length(regions_list))
	for (count in 1:length(regions_list)) {
		chr_list[count] <- unlist(strsplit (regions_list[count], split=":")) [1]
		start_list[count] <- unlist(strsplit(unlist(strsplit(regions_list[count], split=":"))[2], split="-")) [1]
		end_list[count] <- unlist(strsplit(unlist(strsplit(regions_list[count], split=":"))[2], split="-")) [2]	
	}
	dat2 <- data.frame(chr=chr_list, start=start_list, end=end_list, dat)
	
	# Rename gene to symbol for compability with venn diagram
	# colnames (dat2) [5] <- "ID"
	colnames (dat2) [6] <- "symbol"
	colnames (dat2) [11] <- "FPKM_1"
	colnames (dat2) [12] <- "FPKM_2"
	colnames (dat2) [13] <- "log2_FC"

	# Filter the gene output based on status and significant columns
	dat2 <- dat2[dat2$status=="OK",]
	results_list <- dat2
	results_list <- results_list[results_list$significant=="yes",]
	
	# Write output only if signicant results exist
	if (dim(results_list)[1] > 0) {
	
		# order according to increasing q-value
		results_list <- results_list[order(results_list$q_value, decreasing=FALSE),]
		number_genes <- dim (results_list) [1]
		row_names <- 1:number_genes
		rownames(results_list) <- row_names
		write.table(results_list, file="de-isoforms-cufflinks.tsv", sep="\t", row.names=TRUE, col.names=T, quote=F)
	
		# Also output a BED file for visualization and region matching tools
		bed_output <- results_list[,c("chr","start","end","gene_id","log2_FC")]
		write.table(bed_output, file="", sep="\t", row.names=F, col.names=T, quote=F)
		# sort according to chromosome location
		write.table(bed_output, file="sortme.bed", sep="\t", row.names=F, col.names=T, quote=F)
		bed <- read.table(file="sortme.bed", skip=1, sep="\t") # assume file has 1 line header
		colnames(bed)[1:2] <- c("chr", "start")  # these named columns are required for sorting 
		sorted.bed <- sort.bed(bed)
		write.table(sorted.bed, file="de-isoforms-cufflinks.bed", sep="\t", row.names=F, col.names=F, quote=F)
	}
	
	# Report numbers to the log file
	sink(file="cufflinks-log.txt", append=TRUE)
	number_genes_tested <- dim(dat)[1]
	if (dim(results_list)[1] > 0) {
		number_filtered <- number_genes_tested-dim(results_list)[1]
		number_significant <- dim(results_list)[1]
		cat("\n\nTRANSCRIPT ISOFORMS TEST SUMMARY\n")
		cat("In total,", number_genes_tested, "transcript isoforms were tested for differential expression.\n")
		cat("Of these,", number_filtered, "didn't fulfill the technical criteria for testing or the significance cut-off specified.\n")
		cat(number_significant, "transcripts were found to be statistically significantly differentially expressed.\n")	
	} else {
		cat("\n\nTRANSCRIPT ISOFORMS TEST SUMMARY\n")
		cat("Out of the", number_genes_tested, "transcripts tested, there were no statistically significantly differentially expressed ones found.\n")
	}
	sink()
}

# Rename files
# Non-processed cuffdiff output shown only when complete output selected
if (output.type == "complete"){
	if (file.exists("tmp/cds.count_tracking") && file.info("tmp/cds.count_tracking")$size > 12) {
		system("mv tmp/cds.count_tracking cds.count_tracking.tsv")
	}
	if (file.exists("tmp/cds.diff") && file.info("tmp/cds.diff")$size > 115) {
		system("mv tmp/cds.diff cds.diff.tsv")
	}
	if (file.exists("tmp/cds.fpkm_tracking") && file.info("tmp/cds.fpkm_tracking")$size > 91) {
		system("mv tmp/cds.fpkm_tracking cds.fpkm_tracking.tsv")
	}
	if (file.exists("tmp/cds.read_group_tracking") && file.info("tmp/cds.read_group_tracking")$size > 115) {
		system("mv tmp/cds.read_group_tracking cds.read_group_tracking.tsv")
	}
	if (file.exists("tmp/cds_exp.diff") && file.info("tmp/cds_exp.diff")$size > 124) {
		system("mv tmp/cds_exp.diff cds_exp.diff.tsv")
	}
	if (file.exists("tmp/gene_exp.diff") && file.info("tmp/gene_exp.diff")$size > 124) {
		system("mv tmp/gene_exp.diff gene_exp.diff.tsv")
	}
	if (file.exists("tmp/genes.count_tracking") && file.info("tmp/genes.count_tracking")$size > 184) {
		system("mv tmp/genes.count_tracking genes.count_tracking.tsv")
	}
	if (file.exists("tmp/genes.fpkm_tracking") && file.info("tmp/genes.fpkm_tracking")$size > 171) {
		system("mv tmp/genes.fpkm_tracking genes.fpkm_tracking.tsv")
	}
	if (file.exists("tmp/genes.read_group_tracking") && file.info("tmp/genes.read_group_tracking")$size > 115) {
		system("mv tmp/genes.read_group_tracking genes.read_group_tracking.tsv")
	}
	if (file.exists("tmp/isoform_exp.diff") && file.info("tmp/isoform_exp.diff")$size > 124) {
		system("mv tmp/isoform_exp.diff isoform_exp.diff.tsv")
	}
	if (file.exists("tmp/isoforms.count_tracking") && file.info("tmp/isoforms.count_tracking")$size > 184) {
		system("mv tmp/isoforms.count_tracking isoforms.count_tracking.tsv")
	}
	if (file.exists("tmp/isoforms.fpkm_tracking") && file.info("tmp/isoforms.fpkm_tracking")$size > 171) {
		system("mv tmp/isoforms.fpkm_tracking isoforms.fpkm_tracking.tsv")
	}
	if (file.exists("tmp/isoforms.read_group_tracking") && file.info("tmp/isoforms.read_group_tracking")$size > 115) {
		system("mv tmp/isoforms.read_group_tracking isoforms.read_group_tracking.tsv")
	}
	if (file.exists("tmp/promoters.diff") && file.info("tmp/promoters.diff")$size > 115) {
		system("mv tmp/promoters.diff promoters.diff.tsv")
	}
	if (file.exists("tmp/read_groups.info") && file.info("tmp/read_groups.info")$size > 0) {
		system("mv tmp/read_groups.info read_groups.info.txt")
	}
	if (file.exists("tmp/run.info") && file.info("tmp/run.info")$size > 0) {
		system("mv tmp/run.info run.info.txt")
	}
	if (file.exists("tmp/splicing.diff") && file.info("tmp/splicing.diff")$size > 115) {
		system("mv tmp/splicing.diff splicing.diff.tsv")
	}
	if (file.exists("tmp/tss_group_exp.diff") && file.info("tmp/tss_group_exp.diff")$size > 124) {
		system("mv tmp/tss_group_exp.diff tss_group_exp.diff.tsv")
	}
	if (file.exists("tmp/tss_groups.count_tracking") && file.info("tmp/tss_groups.count_tracking")$size > 12) {
		system("mv tmp/tss_groups.count_tracking tss_groups.count_tracking.tsv")
	}
	if (file.exists("tmp/tss_groups.fpkm_tracking") && file.info("tmp/tss_groups.fpkm_tracking")$size > 91) {
		system("mv tmp/tss_groups.fpkm_tracking tss_groups.fpkm_tracking.tsv")
	}
	if (file.exists("tmp/tss_groups.read_group_tracking") && file.info("tmp/tss_groups.read_group_tracking")$size > 115) {
		system("mv tmp/tss_groups.read_group_tracking tss_groups.read_group_tracking.tsv")
	}
}