# TOOL cuffdiff2.R: "Cuffdiff2" (Cuffdiff2)
# INPUT treatment1.bam: "Treatment BAM" TYPE GENERIC
# INPUT control1.bam: "Control BAM" TYPE GENERIC
# INPUT OPTIONAL annotation.gtf: "Annotation GTF" TYPE GENERIC
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
# PARAMETER normalize: "Normalize" TYPE [yes, no] DEFAULT no (Normalize.)
# PARAMETER bias: "Bias correction" TYPE [yes, no] DEFAULT no (Bias detection and correction.)
# PARAMETER genome: "Genome" TYPE [hg19: "Human genome (hg19\)", mm9: "Mouse genome (mm9\)", mm10: "Mouse genome (mm10\)", rn4: "Rat genome (rn4\)"] DEFAULT hg19 (Genome used for bias correction.)
# PARAMETER internalgtf: "Annotation GTF" TYPE [hg19: "Human (hg19\)", mm9: "Mouse (mm9\)", mm10: "Mouse (mm10\)", rn4: "Rat (rn4\)"] DEFAULT hg19 (You can use your own GTF or select one of the provided ones.)



# binary
cuffdiff.binary <- c(file.path(chipster.tools.path, "cufflinks-2.0.2.Linux_x86_64", "cuffdiff"))

# options
cuffdiff.options <- ""
if (normalize == "yes") {
	cuffdiff.options <- paste(cuffdiff.options, "--upper-quartile-norm --total-hits-norm")
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
	genomefile <- c(file.path(chipster.tools.path, "genomes", "fasta", "nochr", genomefile))
	cuffdiff.options <- paste(cuffdiff.options, "-b", genomefile)
}
if (file.exists("annotation.gtf")){
	cuffdiff.options <- paste(cuffdiff.options, "annotation.gtf")
}else{
	if (internalgtf == "hg19") {
		annotation.file <- "Homo_sapiens.GRCh37.68.gtf"
	}
	if (internalgtf == "mm9") {
		annotation.file <- "Mus_musculus.NCBIM37.62.gtf"
	}
	if (internalgtf == "mm10") {
		annotation.file <- "Mus_musculus.GRCm38.68.gtf"
	}
	if (internalgtf == "rn4") {
		annotation.file <- "Rattus_norvegicus.RGSC3.4.68.gtf"
	}
	annotation.file <- c(file.path(chipster.tools.path, "genomes", "gtf", annotation.file))
	cuffdiff.options <- paste(cuffdiff.options, annotation.file)
}


# command
command <- paste(cuffdiff.binary, "-q", "-o tmp", cuffdiff.options, "treatment1.bam", "control1.bam")
# run
#stop(paste('CHIPSTER-NOTE: ', command))
system(command)

# Rename files
if (file.info("tmp/cds.count_tracking")$size > 12) {
	system("mv tmp/cds.count_tracking cds.count_tracking.tsv")
}
if (file.info("tmp/cds.diff")$size > 115) {
	system("mv tmp/cds.diff cds.diff.tsv")
}
if (file.info("tmp/cds.fpkm_tracking")$size > 91) {
	system("mv tmp/cds.fpkm_tracking cds.fpkm_tracking.tsv")
}
if (file.info("tmp/cds.read_group_tracking")$size > 115) {
	system("mv tmp/cds.read_group_tracking cds.read_group_tracking.tsv")
}
if (file.info("tmp/cds_exp.diff")$size > 124) {
	system("mv tmp/cds_exp.diff cds_exp.diff.tsv")
}
if (file.info("tmp/gene_exp.diff")$size > 124) {
	system("mv tmp/gene_exp.diff gene_exp.diff.tsv")
}
if (file.info("tmp/genes.count_tracking")$size > 184) {
	system("mv tmp/genes.count_tracking genes.count_tracking.tsv")
}
if (file.info("tmp/genes.fpkm_tracking")$size > 171) {
	system("mv tmp/genes.fpkm_tracking genes.fpkm_tracking.tsv")
}
if (file.info("tmp/genes.read_group_tracking")$size > 115) {
	system("mv tmp/genes.read_group_tracking genes.read_group_tracking.tsv")
}
if (file.info("tmp/isoform_exp.diff")$size > 124) {
	system("mv tmp/isoform_exp.diff isoform_exp.diff.tsv")
}
if (file.info("tmp/isoforms.count_tracking")$size > 184) {
	system("mv tmp/isoforms.count_tracking isoforms.count_tracking.tsv")
}
if (file.info("tmp/isoforms.fpkm_tracking")$size > 171) {
	system("mv tmp/isoforms.fpkm_tracking isoforms.fpkm_tracking.tsv")
}
if (file.info("tmp/isoforms.read_group_tracking")$size > 115) {
	system("mv tmp/isoforms.read_group_tracking isoforms.read_group_tracking.tsv")
}
if (file.info("tmp/promoters.diff")$size > 115) {
	system("mv tmp/promoters.diff promoters.diff.tsv")
}
if (file.info("tmp/read_groups.info")$size > 0) {
	system("mv tmp/read_groups.info read_groups.info.txt")
}
if (file.info("tmp/run.info")$size > 0) {
	system("mv tmp/run.info run.info.txt")
}
if (file.info("tmp/splicing.diff")$size > 115) {
	system("mv tmp/splicing.diff splicing.diff.tsv")
}
if (file.info("tmp/tss_group_exp.diff")$size > 124) {
	system("mv tmp/tss_group_exp.diff tss_group_exp.diff.tsv")
}
if (file.info("tmp/tss_groups.count_tracking")$size > 12) {
	system("mv tmp/tss_groups.count_tracking tss_groups.count_tracking.tsv")
}
if (file.info("tmp/tss_groups.fpkm_tracking")$size > 91) {
	system("mv tmp/tss_groups.fpkm_tracking tss_groups.fpkm_tracking.tsv")
}
if (file.info("tmp/tss_groups.read_group_tracking")$size > 115) {
	system("mv tmp/tss_groups.read_group_tracking tss_groups.read_group_tracking.tsv")
}



#binary <- "ls -l > list.txt"