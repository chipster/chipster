# TOOL dea-cufflinks.R: "Differential expression analysis using Cufflinks"  (This tool will perform an analysis for differentially expressed genes and isoforms using the Cufflinks algorithm.)
# INPUT treatment.bam: "BAM data file for the treatment sample" TYPE GENERIC
# INPUT control.bam: "BAM data file for the control sample" TYPE GENERIC
# OUTPUT cufflinks-log.txt
# OUTPUT de-genes.tsv
# OUTPUT de-isoforms.tsv
# PARAMETER genome: "Genome" TYPE [hg19: "Human (hg19\)", mm9: "Mouse (mm9\)", rn4: "Rat (rn4\)"] DEFAULT mm9 (Genome that your reads were aligned against.)
# PARAMETER fold.change.threshold: "Fold change cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0 (The cutoff for differential expression. Notice that the fold changes are reported using base 2 logarithmic scale, so the cutoff for finding 2-fold regulated genes should be given as 1.)
# PARAMETER p.value.threshold: "P-value cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 1 (The cutoff for statistical significance. Since the p-values are not adjusted to account for multiple testing correction the cutoff needs to be substantially more conservative than what is usually applied.)
# PARAMETER q.value.threshold: "Q-value cutoff" TYPE DECIMAL FROM 0 TO 1 DEFAULT 1 (The cutoff for statistical significance. Notice that q-values are adjusted to account for multiple testing correction.)


############################################################
#                                                          #
# Analaysis workflow using Cufflinks for normalization and #
# statistical testing for finding differentially expressed #
# sequence tags mapping to genes and transcript isoforms.  #
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
cufflinks.command <- paste(command.start, cufflinks.parameters, cufflinks.input.treatment, cufflinks.input.control, " > cufflinks-log.txt")
system(cufflinks.command)

# Rename output files for Chipster
system ("mv cds_exp.diff de-cds.tsv")
system ("mv gene_exp.diff de-genes.tsv")
system ("mv isoform_exp.diff de-isoforms.tsv")
system ("mv promoters.diff de-promoters.tsv")
system ("mv splicing.diff de-splicing.tsv")
system ("mv tss_group_exp.diff de-tss.tsv")

# Filter the output based on user defined cutoffs
dat <- read.table(file="de-genes.tsv", header=T, sep="\t")
if (fold.change.threshold != 0) {
	dat2 <- data [dat$ln(fold_change) >= fold.change.threshold,]
	dat3 <- data [dat$ln(fold_change) <= fold.change.threshold,]
	dat4 <- rbind (dat2,dat3)
	write.table(dat4, file="de-genes.tsv", sep="\t", row.names=F, col.names=T, quote=F)
}
if (p.value.threshold != 1) {
	dat2 <- data [dat$p_value <= p.value.threshold,]
	write.table(dat2, file="de-genes.tsv", sep="\t", row.names=F, col.names=T, quote=F)
	}
if (q.value.threshold != 1) {
	dat2 <- data [dat$q_value <= q.value.threshold,]
	write.table(dat2, file="de-genes.tsv", sep="\t", row.names=F, col.names=T, quote=F)
}

# EOF

