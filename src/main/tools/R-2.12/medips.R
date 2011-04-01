# TOOL medips.R: "MEDIPS - methylation analysis" (Methylation analysis for sequencing data.)
# INPUT MEDIPS-input{...}.tsv: "Converted BAM data files" GENERIC
# INPUT phenodata.tsv: "Phenodata describing the experiment" GENERIC 
# OUTPUT saturationplot.png: "Saturation plot"
# OUTPUT coverageplot.png: "Coverage plot"
# OUTPUT calibrationplot.png: "Calibration plot"
# OUTPUT CpGdensities.png: "CpG density plot"
# OUTPUT RPMsignal.png: "RPM plot"
# OUTPUT AMSsignal.png: "AMS plot"
# OUTPUT rois.tsv: "Promoter regions"
# OUTPUT output.WIG: "Wiggle file"
# OUTPUT output.tsv: "Enrichment data"
# PARAMETER species [human] DEFAULT human (Select the species)
# PARAMETER promoters.only [yes, no] DEFAULT no (Should the analyses be restricted to promoter regions only)
# PARAMETER coverage.resolution [25,50,100,200] DEFAULT  50 (Targeted data resolution, in base pairs, when the genome-wide coverage is calculated)
# PARAMETER smoothing.extension [200,400,800,1200] DEFAULT 400 (The amount of data smoothing, in base pairs, by extending the reads)
# PARAMETER fragment.length [400,800,1600,2400] DEFAULT 800 (The amount of data smoothing, in base pairs, by extending the reads)
# PARAMETER save.wiggle [yes, no] DEFAULT yes (Should the normalized data, as rpm, be saved as a wiggle file)
# PARAMETER promoter.upstream [1000,2000,5000] DEFAULT 1000 (How much upstream, in base pairs, from the transcription start site does the promoter extend)
# PARAMETER promoter.downstream [250,500,750,1000] DEFAULT 500 (How much downstream, in base pairs, from the transcription start site does the promoter extend)
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Parameters for testing purposes only
#species<-"human"
#promoters.only<-"no"
#coverage.resolution<-"100"
#smoothing.extension<-"1200"
#fragment.length<-2400
#save.wiggle<-"yes"
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

# Reading phenodata
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")

# Analysis type
analtype<-length(unique(phenodata$group ))
repl<-table(phenodata$group)
repltype<-length(unique(repl))

if(analtype==1 & repltype==1) {
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
bitmap(file="coverageplot.png", width=w/72, height=h/72)
par(mar=c(5,5,5,5))
MEDIPS.plotCoverage(cr)
dev.off()

# Coupling vector
dat<-MEDIPS.couplingVector(data=dat, fragmentLength=fragment.length, func="count")

# Calibration curve plot
dat<-MEDIPS.calibrationCurve(data=dat)
bitmap(file="calibrationplot.png", width=w/72, height=h/72)
par(mar=c(5,5,5,5))
MEDIPS.plotCalibrationPlot(dat)
dev.off()

# Normalization
dat<-MEDIPS.normalize(data=dat)

# Save the wiggle file
if(save.wiggle=="yes") {
   MEDIPS.exportWIG(file="output.WIG", data=dat, raw=T, descr=files)
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
   write.table(frames, "output.tsv", sep="\t", quote=F, col.names=T, row.names=F)
}

if(promoters.only=="no") {
   frames<-MEDIPS.methylProfiling(data1=dat, frame_size=as.numeric(smoothing.extension), step=as.numeric(smoothing.extension)/2, math=mean, select=2)
   write.table(frames, "output.tsv", sep="\t", quote=F, col.names=T, row.names=F)
}

# Some plots of the enrichment results
bitmap(file="CpGdensities.png", width=w/72, height=h/72)
par(mar=c(5,5,5,5))
hist(frames$coupling, breaks=100, main="CpG densities", xlab="CpG coupling factors") 
dev.off()

bitmap(file="RPMsignal.png", width=w/72, height=h/72)
par(mar=c(5,5,5,5))
hist(frames$rpm_A[frames$rpm_A!=0], breaks=100, main="RPM signals", xlab="reads/bin") 
dev.off()

bitmap(file="AMSsignal.png", width=w/72, height=h/72)
par(mar=c(5,5,5,5))
hist(frames$ams_A[frames$ams_A!=0], breaks=100, main="AMS signals", xlab="absolute methylation score (ams)") 
dev.off()

# Do not delete the following curly bracket - it ends the if-clause for analysis and replication types
}


