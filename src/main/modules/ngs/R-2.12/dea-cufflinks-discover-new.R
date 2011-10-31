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
# PARAMETER genome: "Genome or transcriptome" TYPE [hg19: "Human genome (hg19\)", mm9: "Mouse genome (mm9\)", rn4: "Rat genome (rn4\)"] DEFAULT mm9 (Genome or transcriptome that your reads were aligned against.)


############################################################
#                                                          #
# Analaysis workflow using Cufflinks for normalization and #
# statistical testing for finding differentially expressed #
# sequence tags mapping to novel genes, transcript         #
# isoforms and splice junctions.                           #
#                                                          #
# The tool assumes that all samples belonging to each      #
# experiment condition have been merged into one single    #
# BAM file.                                                #
#                                                          #
# MG, 21.6.2011                                            #
#                                                          #
############################################################



# Map the reads to reference genome using TopHat
tophat -r 50 -o tophat_brain /seqdata/indexes/hg19 brain_1.fq brain_2.fq
tophat -r 50 -o tophat_liver /seqdata/indexes/hg19 liver_1.fq liver_2.fq
tophat -r 50 -o tophat_heart /seqdata/indexes/hg19 heart_1.fq heart_2.fq 

# Assemble genes and transcript info for each sample using Cufflinks
cufflinks -o cufflinks_brain tophat_brain/accepted_hits.bam
cufflinks -o cufflinks_liver tophat_liver/accepted_hits.bam
cufflinks -o cufflinks_heart tophat_liver/accepted_hits.bam

# Merge the assemblies into a single gtf file
cufflinks_brain/transcripts.gtf
cufflinks_liver/transcripts.gtf
cufflinks_heart/transcripts.gtf
cuffmerge -s /seqdata/fastafiles/hg19/hg19.fa assemblies.txt

# Optional, compare merged assembly to known genes from ENSEMBL annotations
cuffcompare -s /seqdata/fastafiles/hg19/hg19.fa -r known_annotation.gtf merged_asm/merged.gtf

# Run differential expression analysis with gene and transcript discovery
# Take the merged assembly from produced in step 3 of the discovery protocol and provide it to cuffdiff along with the BAM files from TopHat:
cuffdiff merged_asm/merged.gtf liver1.bam,liver2.bam brain1.bam,brain2.bam

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

# EOF

