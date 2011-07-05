# TOOL dea-cufflinks.R: "Differential expression analysis using Cufflinks"  (This tool will perform an analysis for differentially expressed genes and isoforms using the Cufflinks algorithm.)
# INPUT treatment.bam: "BAM data file for the treatment sample" TYPE GENERIC
# INPUT control.bam: "BAM data file for the control sample" TYPE GENERIC
# OUTPUT cufflinks-log.txt
# OUTPUT de-cds.tsv
# OUTPUT de-genes.tsv
# OUTPUT de-isoforms.tsv
# OUTPUT de-promoters.tsv
# OUTPUT de-splicing.tsv
# OUTPUT de-tss.tsv
# PARAMETER genome: "Genome" TYPE [hg19: "Human (hg19\)", mm9: "Mouse (mm9\)", rn4: "Rat (rn4\)"] DEFAULT mm9 (Genome that your reads were aligned against.)


############################################################
#                                                          #
# Analaysis workflow using Cufflinks for normalization and #
# statistical testing for finding differentially expressed #
# sequence tags mapping to genes and transcript isoforms   #
#                                                          #
# MG, 21.6.2011                                            #
# development version, 1 sample vs 1 sample                #
#                                                          #
############################################################

# Cufflinks tools setup
cufflinks.binary <- c(file.path(chipster.tools.path, "cufflinks", "cuffdiff"))
command.start <- cufflinks.binary

# Annotation file setup
annotation.path <- c(file.path(chipster.tools.path, "genomes"))
if (genome == "hg19") {
	annotation.file <- "homo_sampiens/annotations/Homo_sapiens.GRCh37.62.gtf"
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
sink(file="cufflinks-log.txt")
system(cufflinks.command)
sink()

# Rename output files for Chipster
system ("mv cds_exp.diff de-cds.tsv")
system ("mv gene_exp.diff de-genes.tsv")
system ("mv isoform_exp.diff de-isoforms.tsv")
system ("mv promoters.diff de-promoters.tsv")
system ("mv splicing.diff de-splicing.tsv")
system ("mv tss_group_exp.diff de-tss.tsv")

# EOF

