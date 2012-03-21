# TOOL dea-cufflinks.R: "Differential expression analysis using Cufflinks"  (This tool will perform an analysis for differentially expressed genes and isoforms using the Cufflinks algorithm. Note that only one filtering criteria should be applied for a given analysis run. When left at default settings, Cufflinks filters out unsuccessfully tested loci, as well as those with a Benjamini-Hochberg adjusted false discovery rate less than 0.05.)
# INPUT treatment.bam: "BAM data file for the treatment sample" TYPE GENERIC
# INPUT control.bam: "BAM data file for the control sample" TYPE GENERIC
# OUTPUT cufflinks-log.txt
# OUTPUT de-genes-cufflinks.tsv
# OUTPUT de-isoforms-cufflinks.tsv
# OUTPUT OPTIONAL de-genes-cufflinks.bed
# OUTPUT OPTIONAL de-isoforms-cufflinks.bed
# PARAMETER genome: "Genome" TYPE [hg19: "Human (hg19\)", mm9: "Mouse (mm9\)", rn4: "Rat (rn4\)"] DEFAULT mm9 (Genome that your reads were aligned against.)
# PARAMETER fold.change.threshold: "Fold change cutoff" TYPE DECIMAL FROM 0 TO 1000000 DEFAULT 0 (The cutoff for differential expression. Note that the fold changes are reported using base 2 logarithmic scale, so the cutoff for finding 2-fold regulated genes should be given as 1.)
# PARAMETER p.value.threshold: "P-value cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 1 (The cutoff for statistical significance. Since the p-values are not adjusted to account for multiple testing correction, the cutoff needs to be substantially more conservative than what is usually applied.)
# PARAMETER q.value.threshold: "Q-value cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 1 (The cutoff for statistical significance. Note that q-values are adjusted to account for multiple testing correction.)


############################################################
#                                                          #
# Analaysis workflow using Cufflinks for normalization and #
# statistical testing for finding differentially expressed #
# known genes and transcript isoforms.                     #
#                                                          #
# The tool assumes that all samples belonging to each      #
# experiment condition have been merged into one single    #
# BAM file.                                                #
#                                                          #
# MG, 21.6.2011                                            #
#                                                          #
############################################################

# Output that is yet to be supported
# OUTPUT de-cds.tsv
# OUTPUT de-tss.tsv
# OUTPUT de-splicing.tsv
# OUTPUT de-promoters.tsv

# Cufflinks tools setup
cufflinks.binary <- c(file.path(chipster.tools.path, "cufflinks", "cuffdiff"))
command.start <- cufflinks.binary

# Annotation file setup
annotation.path <- c(file.path(chipster.tools.path, "genomes"))
if (genome == "hg19") {
	annotation.file <- "homo_sapiens/annotations/Homo_sapiens.GRCh37.62.gtf"
}
if (genome == "mm9") {
	annotation.file <- "mus_musculus/annotations/Mus_musculus.NCBIM37.62.gtf"
}
if (genome == "rn4") {
	annotation.file <- "rattus_norvegicus/annotations/Rattus_norvegicus.RGSC3.4.62.gtf"
}
annotation.file <- c(file.path(chipster.tools.path, "genomes", annotation.file))

# Run differential expression analysis for known genes and transcript isoforms
cufflinks.parameters <- annotation.file
cufflinks.input.treatment <- "treatment.bam"
cufflinks.input.control <- "control.bam"
cufflinks.command <- paste(command.start, cufflinks.parameters, cufflinks.input.treatment, cufflinks.input.control)
system(cufflinks.command)

# Rename output files for Chipster
system ("mv gene_exp.diff de-genes.tsv")
system ("mv isoform_exp.diff de-isoforms.tsv")
# system ("mv cds_exp.diff de-cds.tsv")
# system ("mv promoters.diff de-promoters.tsv")
# system ("mv splicing.diff de-splicing.tsv")
# system ("mv tss_group_exp.diff de-tss.tsv")

# DE genes
# Extract chromosome locations and add in the first three table columns
dat <- read.table(file="de-genes.tsv", header=T, sep="\t")
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
colnames (dat2) [5] <- "ensembl_id"
colnames (dat2) [6] <- "symbol"
colnames (dat2) [13] <- "ln(fold_change)"

# Filter the gene output based on user defined cutoffs
dat2 <- dat2[dat2$status=="OK",]
results_list <- dat2
if (fold.change.threshold != 0 || p.value.threshold < 1 || q.value.threshold < 1) {
	if (fold.change.threshold != 0) {
		dat3 <- dat2 [dat2$ln.fold_change. >= fold.change.threshold,]
		dat4 <- dat2 [dat2$ln.fold_change. <= -fold.change.threshold,]
		results_list <- rbind (dat3,dat4)
	}
	if (p.value.threshold < 1) {
		results_list <- dat2 [dat2$p_value <= p.value.threshold,]
	}
	if (q.value.threshold < 1) {
		results_list <- dat2 [dat2$q_value <= q.value.threshold,]
	}
} else {
	results_list <- results_list[results_list$significant=="yes",]
}
# order according to increasing q-value
results_list <- results_list[order(results_list$q_value, decreasing=FALSE),]
write.table(results_list, file="de-genes-cufflinks.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# Also output a bed graph file for visualization and region matching tools
if (dim(results_list)[1] > 0) {
	bed_output <- results_list[,c("chr","start","end","symbol","ln(fold_change)")]
	# sort according to chromosome location
	bed_output <- bed_output[order(bed_output$chr, bed_output$start, bed_output$end, decreasing=FALSE),]
	# add chr to the chromosome name for genome browser compability
	bed_output[,1] <- paste("chr",bed_output[,1],sep="")
	write.table(bed_output, file="de-genes-cufflinks.bed", sep="\t", row.names=F, col.names=F, quote=F)
}

# Report numbers to the log file
if (dim(results_list)[1] > 0) {
	sink(file="cufflinks-log.txt")
	number_genes_tested <- dim(dat)[1]
	number_filtered <- number_genes_tested-dim(results_list)[1]
	number_significant <- dim(results_list)[1]
	cat("GENE TEST SUMMARY\n")
	cat("In total,", number_genes_tested, "genes were tested for differential expression.\n")
	cat("Of these,", number_filtered, "didn't fulfill the technical criteria for testing or the significance cut-off specified.\n")
	cat(number_significant, "genes were found to be statiscially significantly differentially expressed.")	
} else {
	cat("GENE TEST SUMMARY\n")
	cat("Out of the", number_genes_tested, "genes tested there were no statistically significantly differentially expressed ones found.")
}

# DE isoforms
# Extract chromosome locations and add in the first three table columns
dat <- read.table(file="de-isoforms.tsv", header=T, sep="\t")
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
colnames (dat2) [5] <- "ensembl_id"
colnames (dat2) [6] <- "symbol"
colnames (dat2) [13] <- "ln(fold_change)"

# Filter the isoforms output based on user defined cutoffs
dat2 <- dat2[dat2$status=="OK",]
results_list <- dat2
if (fold.change.threshold != 0 || p.value.threshold < 1 || q.value.threshold < 1) {
	if (fold.change.threshold != 0) {
		dat3 <- dat2 [dat2$ln.fold_change. >= fold.change.threshold,]
		dat4 <- dat2 [dat2$ln.fold_change. <= -fold.change.threshold,]
		results_list <- rbind (dat3,dat4)
	}
	if (p.value.threshold < 1) {
		results_list <- dat2 [dat2$p_value <= p.value.threshold,]
	}
	if (q.value.threshold < 1) {
		results_list <- dat2 [dat2$q_value <= q.value.threshold,]
	}
} else {
	results_list <- results_list[results_list$significant=="yes",]
}
# order according to increasing q-value
results_list <- results_list[order(results_list$q_value, decreasing=FALSE),]
write.table(results_list, file="de-isoforms-cufflinks.tsv", sep="\t", row.names=F, col.names=T, quote=F)

# Also output a bed graph file for visualization and region matching tools
if (dim(results_list)[1] > 0) {
	bed_output <- results_list[,c("chr","start","end","symbol","ln(fold_change)")]
	# sort according to chromosome location
	bed_output <- bed_output[order(bed_output$chr, bed_output$start, bed_output$end, decreasing=FALSE),]
	# add chr to the chromosome name for genome browser compability
	bed_output[,1] <- paste("chr",bed_output[,1],sep="")
	write.table(bed_output, file="de-isoforms-cufflinks.bed", sep="\t", row.names=F, col.names=F, quote=F)
}

# Report numbers to the log file
if (dim(results_list)[1] > 0) {
	number_genes_tested <- dim(dat)[1]
	number_filtered <- number_genes_tested-dim(results_list)[1]
	number_significant <- dim(results_list)[1]
	cat("\n\nTRANSCRIPT ISOFORMS TEST SUMMARY\n")
	cat("In total,", number_genes_tested, "transcript isoforms were tested for differential expression.\n")
	cat("Of these,", number_filtered, "didn't fulfill the technical criteria for testing or the significance cut-off specified.\n")
	cat(number_significant, "transcripts were found to be statiscially significantly differentially expressed.")	
} else {
	cat("\n\nTRANSCRIPT ISOFORMS TEST SUMMARY\n")
	cat("Out of the", number_genes_tested, "transcripts tested there were no statistically significantly differentially expressed ones found.")
}
sink()

# EOF

