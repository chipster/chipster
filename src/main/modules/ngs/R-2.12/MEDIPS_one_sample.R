# TOOL MEDIPS_one_sample.R: "MEDIPS - methylation analysis, treatment only" (Methylation analysis for sequencing data. Analysis for a single condition or file.)
# INPUT MEDIPS-input.tsv: "Converted BAM data files in MEDIPS format" TYPE GENERIC
# OUTPUT OPTIONAL saturationplot.png: "Saturation plot"
# OUTPUT OPTIONAL coverageplot.pdf: "Coverage plot"
# OUTPUT OPTIONAL calibrationplot.pdf: "Calibration plot"
# OUTPUT OPTIONAL CpGdensities.pdf: "CpG density plot"
# OUTPUT OPTIONAL RPMsignal.pdf: "RPM plot"
# OUTPUT OPTIONAL AMSsignal.pdf: "AMS plot"
# OUTPUT OPTIONAL rois.tsv: "Promoter regions"
# OUTPUT OPTIONAL methylation.bed: "BED file"
# OUTPUT methylation.tsv: "Enrichment data"
# PARAMETER species TYPE [human] DEFAULT human (Select the species)
# PARAMETER promoters.only TYPE [yes, no] DEFAULT no (Should the analyses be restricted to promoter regions only)
# PARAMETER fragment.length TYPE [400,800,1600,2400] DEFAULT 800 (Fragment length, used for calculating local CpGs)
# PARAMETER OPTIONAL coverage.resolution TYPE [25,50,100,200] DEFAULT  50 (Targeted data resolution, in base pairs, when the genome-wide coverage is calculated)
# PARAMETER OPTIONAL smoothing.extension TYPE [200,400,800,1200] DEFAULT 400 (The amount of data smoothing, in base pairs, by extending the reads)
# PARAMETER OPTIONAL save.bed TYPE [yes, no] DEFAULT yes (Should the normalized data, as RPM, be saved as a BED file)
# PARAMETER OPTIONAL promoter.upstream TYPE [1000,2000,5000] DEFAULT 1000 (How much upstream, in base pairs, from the transcription start site does the promoter extend)
# PARAMETER OPTIONAL promoter.downstream TYPE [250,500,750,1000] DEFAULT 500 (How much downstream, in base pairs, from the transcription start site does the promoter extend)
# PARAMETER OPTIONAL image.width TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted images in pixels)
# PARAMETER OPTIONAL image.height TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted images in pixels)


# Parameters for testing purposes only
#species<-"human"
#promoters.only<-"no"
#coverage.resolution<-"100"
#smoothing.extension<-"1200"
#fragment.length<-2400
#save.bed<-"yes"
#promoter.upstream<-"1000"
#promoter.downstream<-"500"
#image.width<-600
#image.height<-600

# Processing of the parameters
if(species=="human") {
   genome<-c("BSgenome.Hsapiens.UCSC.hg19")
   pgenome<-"hg19"  
   library(BSgenome.Hsapiens.UCSC.hg19)
}
w<-image.width
h<-image.height

# Input files
files<-dir(pattern="MEDIPS-input")





# Reading the MEDIPS input file
library(MEDIPS)
dat<-MEDIPS.readAlignedSequences(BSgenome=genome, file=files)

# Creating the genome vector
dat<-MEDIPS.genomeVector(data=dat, bin_size=as.numeric(coverage.resolution), extend=as.numeric(smoothing.extension))

# Cleaning up to save some memory
gc()

# Saturation analysis
sr<-MEDIPS.saturationAnalysis(data=dat, bin_size=as.numeric(coverage.resolution), extend=as.numeric(smoothing.extension), no_iterations=10, no_random_iterations=1)
bitmap(file="saturationplot.png", width=w/72, height=h/72)
par(mar=c(5,5,5,5))
MEDIPS.plotSaturation(sr)
dev.off()

# Pattern Positions
dat<-MEDIPS.getPositions(data=dat, pattern="CG")

# Coverage analysis
cr<-MEDIPS.coverageAnalysis(data=dat, extend=as.numeric(smoothing.extension), no_iterations=10)
pdf(file="coverageplot.pdf", width=w/72, height=h/72)
par(mar=c(5,5,5,5))
MEDIPS.plotCoverage(cr)
dev.off()

# Coupling vector
dat<-MEDIPS.couplingVector(data=dat, fragmentLength=as.numeric(fragment.length), func="count")

# Calibration curve plot
dat<-MEDIPS.calibrationCurve(data=dat)
pdf(file="calibrationplot.pdf", width=w/72, height=h/72)
par(mar=c(5,5,5,5))
MEDIPS.plotCalibrationPlot(dat)
dev.off()

# Normalization
dat<-MEDIPS.normalize(data=dat)

# Save the BED file
if(save.bed=="yes") {
   source(file.path(chipster.common.path, "bed-utils.R"))
   MEDIPS.exportWIG(file="output.bed", data=dat, raw=F, descr=files)
   library(rtracklayer)
   wig<-import("output.bed", format="wig")
   export(wig, "output.bed", format="bed")
   wig<-read.table("output.bed", header=F, sep="\t")
   wig<-wig[wig[,5]>0,]
   wig <- sort.bed(wig)
   write.table(wig, "methylation.bed", col.names=F, row.names=F, sep="\t", quote=F)
}

if(promoters.only=="yes") {
   library(rtracklayer)
   session <- browserSession()
   genome(session) <- pgenome
   query <- ucscTableQuery(session, "refGene")
   refGene<-getTable(query)
   sta<-refGene$txStart
   sta[refGene$strand=="+"]<-sta[refGene$strand=="+"]-1000
   sta[refGene$strand=="-"]<-sta[refGene$strand=="-"]-500
   sto<-refGene$txStart
   sto[refGene$strand=="+"]<-sto[refGene$strand=="+"]+500
   sto[refGene$strand=="-"]<-sto[refGene$strand=="-"]+1000
   rois<-data.frame(refGene$chrom, sta, sto, make.unique(as.character(refGene$name)))
   rois<-rois[rois[,1] %in% unique(dat@genome_chr),]
   write.table(rois, "rois.tsv", sep="\t", quote=F, col.names=F, row.names=F)
   frames<-MEDIPS.methylProfiling(data1=dat, ROI_file="rois.tsv", math=mean, select=2)
   write.table(frames, "methylation.tsv", sep="\t", quote=F, col.names=T, row.names=F)
}

if(promoters.only=="no") {
   frames<-MEDIPS.methylProfiling(data1=dat, frame_size=as.numeric(smoothing.extension), step=as.numeric(smoothing.extension)/2, math=mean, select=2)
   write.table(frames, "methylation.tsv", sep="\t", quote=F, col.names=T, row.names=F)
}

# Some plots of the enrichment results
pdf(file="CpGdensities.pdf", width=w/72, height=h/72)
par(mar=c(5,5,5,5))
hist(frames$coupling, breaks=100, main="CpG densities", xlab="CpG coupling factors") 
dev.off()

pdf(file="RPMsignal.pdf", width=w/72, height=h/72)
par(mar=c(5,5,5,5))
hist(frames$rpm_A[frames$rpm_A!=0], breaks=100, main="RPM signals", xlab="reads/bin") 
dev.off()

pdf(file="AMSsignal.pdf", width=w/72, height=h/72)
par(mar=c(5,5,5,5))
hist(frames$ams_A[frames$ams_A!=0], breaks=100, main="AMS signals", xlab="absolute methylation score (ams)") 
dev.off()
