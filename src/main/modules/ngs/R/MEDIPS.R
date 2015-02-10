# TOOL MEDIPS.R: "MEDIPS - methylation analysis" (Methylation analysis for sequencing data. This tool can be used to analyze a single condition, or also of treatment vs. control, in which case the differences between the two conditions are tested using edgeR and resulting p-values returned. Additionally, an optional input sequencing data file can be provided as well.)
# INPUT treatment.bam: "BAM data file for the treatment" TYPE GENERIC
# INPUT OPTIONAL control.bam: "Optional BAM data file for the control" TYPE GENERIC
# INPUT OPTIONAL input.bam: "Optional BAM data file for the input" TYPE GENERIC
# OUTPUT methylation.tsv: "Enrichment data"
# OUTPUT saturation.pdf: "Saturation plot"
# OUTPUT coverage.pdf: "Coverage plot"
# OUTPUT calibration.pdf: "Calibration plot"
# OUTPUT OPTIONAL treatment.bed: "Enrichment data for treatment"
# OUTPUT OPTIONAL control.bed: "Enrichment data for control"
# OUTPUT OPTIONAL input.bed: "Enrichment data for input"
# PARAMETER chr: "Chromosome names in my BAM file look like" TYPE [chr1, 1] DEFAULT 1 (Chromosome names must match in the BAM file and in the reference annotation. Check your BAM and choose accordingly.)
# PARAMETER genome: "Genome" TYPE [hg19: "human hg19", GRCh38: "human GRCh38", mm10: "mouse mm10", rn5: "rat rn5"] DEFAULT GRCh38 (Select the genome build)
# PARAMETER promoters.only: "Focus on promoter regions only" TYPE [yes, no] DEFAULT no (Should the analyses be restricted to promoter regions only)
# PARAMETER fragment.length: "Fragment length" TYPE [400, 800, 1600, 2400] DEFAULT 800 (Fragment length, used for calculating local CpGs)
# PARAMETER OPTIONAL coverage.resolution: "Coverage resolution" TYPE [25, 50, 100, 200] DEFAULT 50 (Targeted data resolution, in base pairs, when the genome-wide coverage is calculated)
# PARAMETER OPTIONAL smoothing.extension: "Smoothing extension" TYPE [200, 400, 800, 1200] DEFAULT 400 (The amount of data smoothing, in base pairs, by extending the reads)
# PARAMETER OPTIONAL promoter.upstream: "Number of promoter base pairs upstream" TYPE [1000, 2000, 5000] DEFAULT 1000 (How much upstream, in base pairs, from the transcription start site does the promoter extend)
# PARAMETER OPTIONAL promoter.downstream: "Number of promoter base pairs downstream" TYPE [250, 500, 750, 1000] DEFAULT 500 (How much downstream, in base pairs, from the transcription start site does the promoter extend)
# PARAMETER OPTIONAL save.bed: "Save results also as a BED file" TYPE [yes, no] DEFAULT yes (Should the normalized data, as RPM, be saved as a BED file)

# Parameters for testing purposes only
# genome <- "hg19"
# fragment.length <- "800"
# coverage.resolution <- "50"
# smoothing.extension <- "400"
# promoters.only <- "no"
# promoter.upstream <- "1000"
# promoter.downstream <- "500"
# save.bed <- "yes"

# Add chr to chromosome names in BAM if necessary. Optional inputs are processed only if present.
tfile <- "treatment.bam"
if (chr == "1") {
	source(file.path(chipster.common.path, "bam-utils.R"))
	addChrToBAM("treatment.bam", "treatment_chr.bam")
	system("mv treatment_chr.bam treatment.bam")
}
cfile <- "control.bam"
use.control <- file.exists(cfile)
if (use.control){
	if (chr == "1") {
		addChrToBAM("control.bam", "control_chr.bam")
		cfile <- "control_chr.bam"
	}
}
ifile <- "input.bam"
use.input <- file.exists(ifile)
if (use.input){
	if (chr == "1") {
		addChrToBAM("input.bam", "input_chr.bam")
		ifile <- "input_chr.bam"
	}
}

# Processing the parameters
if (genome == "hg19") {
  genome <- "BSgenome.Hsapiens.UCSC.hg19"
  pgenome <- "hg19"  
  library(BSgenome.Hsapiens.UCSC.hg19)
}
if (genome == "GRCh38") {
	genome <- "BSgenome.Hsapiens.NCBI.GRCh38"
	pgenome <- "hg38"  
	library(BSgenome.Hsapiens.NCBI.GRCh38)
}
if (genome == "mm10") {
	genome <- "BSgenome.Mmusculus.UCSC.mm10"
	pgenome <- "mm10"  
	library(BSgenome.Mmusculus.UCSC.mm10)
}
if (genome == "rn5") {
	genome <- "BSgenome.Rnorvegicus.UCSC.rn5"
	pgenome <- "rn5"  
	library(BSgenome.Rnorvegicus.UCSC.rn5)
}

fragment.length <- as.integer(fragment.length)
coverage.resolution <- as.integer(coverage.resolution)
smoothing.extension <- as.integer(smoothing.extension)
promoter.upstream <- as.integer(promoter.upstream)
promoter.downstream <- as.integer(promoter.downstream)

	
# Load library to memory
library(MEDIPS)

# saturation analysis
pdf("saturation.pdf", width=0, height=0, paper="a4r")
sr.treatment <- MEDIPS.saturation(file=tfile, BSgenome=genome, extend=smoothing.extension, window_size=coverage.resolution)
MEDIPS.plotSaturation(sr.treatment, main="Saturation analysis for treatment")
if (use.control) {
  sr.control <- MEDIPS.saturation(file=cfile, BSgenome=genome, extend=smoothing.extension, window_size=coverage.resolution)
  MEDIPS.plotSaturation(sr.control, main="Saturation analysis for control")
}
if (use.input) {
  sr.input <- MEDIPS.saturation(file=ifile, BSgenome=genome, extend=smoothing.extension, window_size=coverage.resolution)
  MEDIPS.plotSaturation(sr.input, main="Saturation analysis for input")
}
dev.off()

# sequence pattern coverage
pdf("coverage.pdf", width=0, height=0, paper="a4r")
cr.treatment <- MEDIPS.seqCoverage(file=tfile, pattern="CG", BSgenome=genome, extend=smoothing.extension)
MEDIPS.plotSeqCoverage(cr.treatment, type="pie")
title("\n\ntreatment")
MEDIPS.plotSeqCoverage(cr.treatment, type="hist")
title("\n\ntreatment")
if (use.control) {
  cr.control <- MEDIPS.seqCoverage(file=cfile, pattern="CG", BSgenome=genome, extend=smoothing.extension)
  MEDIPS.plotSeqCoverage(cr.control, type="pie")
  title("\n\ncontrol")
  MEDIPS.plotSeqCoverage(cr.control, type="hist")
  title("\n\ncontrol")
}
if (use.input) {
  cr.input <- MEDIPS.seqCoverage(file=ifile, pattern="CG", BSgenome=genome, extend=smoothing.extension)
  MEDIPS.plotSeqCoverage(cr.input, type="pie")
  title("\n\ninput")
  MEDIPS.plotSeqCoverage(cr.input, type="hist")
  title("\n\ninput")
}
dev.off()

# Reads the data
treatment <- MEDIPS.createSet(file=tfile, BSgenome=genome,  extend=as.numeric(smoothing.extension), window_size=as.numeric(coverage.resolution))
if (use.control) {
  control <- MEDIPS.createSet(file=cfile, BSgenome=genome, extend=as.numeric(smoothing.extension),  window_size=as.numeric(coverage.resolution))
}
if (use.input) {
  input <- MEDIPS.createSet(file=ifile, BSgenome=genome,  extend=as.numeric(smoothing.extension), window_size=as.numeric(coverage.resolution))
}

# Coupling vector
cs <- MEDIPS.couplingVector(pattern="CG", refObj=treatment)

# calibration
library(png)
tmpfiles <- c(tempfile(), tempfile())
bitmap(tmpfiles[1], width=11.7, height=8.3, units='in', res=300)
if (use.input) {
  MEDIPS.plotCalibrationPlot(MSet=treatment, ISet=input, CSet=cs, main=paste("Calibration plot for treatment, ", treatment@chr_names[1]), plot_chr=treatment@chr_names[1])
} else {
  MEDIPS.plotCalibrationPlot(MSet=treatment, CSet=cs, main=paste("Calibration plot for treatment, ", treatment@chr_names[1]), plot_chr=treatment@chr_names[1])
}
dev.off()
if (use.control) {
  bitmap(tmpfiles[2], width=11.7, height=8.3, units='in', res=300)
  MEDIPS.plotCalibrationPlot(MSet=control, CSet=cs, main=paste("Calibration plot for control, ", control@chr_names[1]), plot_chr=control@chr_names[1])
  dev.off()
}

pdf("calibration.pdf", width=0, height=0, paper="a4r")
par(mai=c(0, 0, 0, 0))
plot.new()
rasterImage(readPNG(tmpfiles[1]), 0, 0, 1, 1)
if (use.control) {
  plot.new()
  rasterImage(readPNG(tmpfiles[2]), 0, 0, 1, 1)
}
dev.off()

if (use.control) {
	if (use.input) {
		mr.edgeR <- MEDIPS.meth(MSet1=control, MSet2=treatment, ISet1=input, CSet=cs, diff.method="edgeR")
	} else {
		mr.edgeR <- MEDIPS.meth(MSet1=control, MSet2=treatment, CSet=cs, diff.method="edgeR")
	}
} else {
	if (use.input) {
		mr.edgeR <- MEDIPS.meth(MSet1=treatment, ISet1=input, CSet=cs, diff.method="edgeR")
	} else {
		mr.edgeR <- MEDIPS.meth(MSet1=treatment, CSet=cs, diff.method="edgeR")
	}
}
colnames(mr.edgeR) <- sub("\\.bam\\.", ".", colnames(mr.edgeR))
mr.edgeR <- mr.edgeR[, grep("\\.mean$", colnames(mr.edgeR), invert=TRUE)]

if (promoters.only == "yes") {
	library(rtracklayer)
	session <- browserSession()
	genome(session) <- pgenome
	query <- ucscTableQuery(session, "refGene")
	refGene <- getTable(query)
	minus <- refGene$strand == "-"
	sta <- refGene$txStart - promoter.upstream
	sta[minus] <- refGene$txEnd[minus] - promoter.downstream
	sto <- refGene$txStart + promoter.downstream
	sto[minus] <- refGene$txEnd[minus] + promoter.upstream
	rois <- data.frame(chromosome=as.character(refGene$chrom), start=sta, end=sto, gene=as.character(refGene$name2), stringsAsFactors=FALSE)
	mr.edgeR <- MEDIPS.selectROIs(results=mr.edgeR, rois=rois)
	colnames(mr.edgeR)[colnames(mr.edgeR) == "ROI"] <- "symbol"
}

# order the result st
chr <- mr.edgeR$chr
chr <- sub("^CHR", "", toupper(chr))
chr[chr == "X"] <- 23
chr[chr == "Y"] <- 24
chr[chr == "M" | chr == "MT"] <- 25
chr <- as.integer(chr)
mr.edgeR <- mr.edgeR[order(chr, mr.edgeR$start), ]

# remove duplicate rows when same coverage window overlaps multiple promoters
conc <- function(x) {
	x <- unique(x)
	paste(x, collapse=";")
}
position <- sprintf("%s:%i-%i", mr.edgeR$chr, mr.edgeR$start, mr.edgeR$stop)
duplicates <- unique(position[duplicated(position)])
for (duplicate in duplicates)
	mr.edgeR[position == duplicate, 'symbol'] <- conc(mr.edgeR[position == duplicate, 'symbol'])
mr.edgeR <- mr.edgeR[!duplicated(position), ]
position <- position[!duplicated(position)]
rownames(mr.edgeR) <- position

for (i in colnames(mr.edgeR))
	if (is.numeric(mr.edgeR[, i]))
		mr.edgeR[, i] <- round(mr.edgeR[, i], digits=3)

options(scipen=10)

if (save.bed == "yes") {
	if ("symbol" %in% colnames(mr.edgeR)) {
		gene <- mr.edgeR$symbol
	} else {
		gene <- position
	}
	bed <- data.frame(chromosome=mr.edgeR$chr, start=mr.edgeR$start, end=mr.edgeR$stop, name=gene, score=mr.edgeR$treatment.rpkm, strand=rep("+", nrow(mr.edgeR)), stringsAsFactors=FALSE)
	bed$start <- bed$start - 1
	write.table(bed, "treatment.bed", quote=FALSE, sep="\t", na="", row.names=FALSE, col.names=FALSE)
	if (use.control) {
		bed$score <- mr.edgeR$control.rpkm
		write.table(bed, "control.bed", quote=FALSE, sep="\t", na="", row.names=FALSE, col.names=FALSE)
	}
	if (use.input) {
		bed$score <- mr.edgeR$input.rpkm
		write.table(bed, "input.bed", quote=FALSE, sep="\t", na="", row.names=FALSE, col.names=FALSE)
	}
}

write.table(mr.edgeR, file="methylation.tsv", quote=FALSE, sep="\t", na="")

# EOF