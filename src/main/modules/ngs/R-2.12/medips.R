# TOOL medips.R: "MEDIPS - methylation analysis" (Methylation analysis for sequencing data.)
# INPUT MEDIPS-input1.tsv: "Converted BAM data files" TYPE GENERIC
# INPUT META phenodata.tsv: "Phenodata describing the experiment" TYPE GENERIC 
# OUTPUT OPTIONAL saturationplot.pdf: "Saturation plot"
# OUTPUT OPTIONAL coverageplot.pdf: "Coverage plot"
# OUTPUT OPTIONAL calibrationplot.png: "Calibration plot"
# OUTPUT OPTIONAL CpGdensities.pdf: "CpG density plot"
# OUTPUT OPTIONAL RPMsignal.pdf: "RPM plot"
# OUTPUT OPTIONAL AMSsignal.pdf: "AMS plot"
# OUTPUT OPTIONAL rois.tsv: "Promoter regions"
# OUTPUT OPTIONAL output.bed: "BED file"
# OUTPUT output.tsv: "Enrichment data"
# PARAMETER species TYPE [human] DEFAULT human (Select the species)
# PARAMETER promoters.only TYPE [yes, no] DEFAULT no (Should the analyses be restricted to promoter regions only)
# PARAMETER coverage.resolution TYPE [25,50,100,200] DEFAULT  50 (Targeted data resolution, in base pairs, when the genome-wide coverage is calculated)
# PARAMETER smoothing.extension TYPE [200,400,800,1200] DEFAULT 400 (The amount of data smoothing, in base pairs, by extending the reads)
# PARAMETER fragment.length TYPE [400,800,1600,2400] DEFAULT 800 (Length of the fragment for calculation of local CpGs)
# PARAMETER save.bed TYPE [yes, no] DEFAULT yes (Should the normalized data, as rpm, be saved as a BED file)
# PARAMETER promoter.upstream TYPE [1000,2000,5000] DEFAULT 1000 (How much upstream, in base pairs, from the transcription start site does the promoter extend)
# PARAMETER promoter.downstream TYPE [250,500,750,1000] DEFAULT 500 (How much downstream, in base pairs, from the transcription start site does the promoter extend)
# PARAMETER image.width TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


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
	genome<-c("BSgenome.Hsapiens.UCSC.hg18")
	pgenome<-"hg19"  
	library(BSgenome.Hsapiens.UCSC.hg18)
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
if(any(repl>1)) {
	repltype<-c(2)
} else {
	repltype<-c(1)
}

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
	pdf(file="saturationplot.pdf", width=w/72, height=h/72)
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
	bitmap(file="calibrationplot.png", width=w/72, height=h/72)
	par(mar=c(5,5,5,5))
	MEDIPS.plotCalibrationPlot(dat)
	dev.off()
	
	# Normalization
	dat<-MEDIPS.normalize(data=dat)
	
	# Save the bed file
	if(save.bed=="yes") {
		MEDIPS.exportWIG(file="output.bed", data=dat, raw=T, descr=files)
		library(rtracklayer)
		wig<-import("output.bed", format="wig")
		export(wig, "output.bed", format="bed")
		wig<-read.table("output.bed", header=F, sep="\t")
		wig<-wig[wig[,5]>0,]
		wig <- sort.bed(wig)
		write.table(wig, "output.bed", col.names=F, row.names=F, sep="\t", quote=F)
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
	
	# Do not delete the following curly bracket - it ends the if-clause for analysis and replication types
}





if(analtype==2 & repltype==1) {
	library(MEDIPS)
	
	# Assign the files to control and treatment based on the group column from phenodata
	cfile<-as.character(phenodata$sample[phenodata$group==head(sort(phenodata$group), n=1)])
	tfile<-as.character(phenodata$sample[phenodata$group==tail(sort(phenodata$group), n=1)])
	
	# Reads the data
	control<-MEDIPS.readAlignedSequences(BSgenome=genome, file=cfile)
	treatment<-MEDIPS.readAlignedSequences(BSgenome=genome, file=tfile)
	
	# Creating the genome vector
	control<-MEDIPS.genomeVector(data=control, bin_size=as.numeric(coverage.resolution), extend=as.numeric(smoothing.extension))
	treatment<-MEDIPS.genomeVector(data=treatment, bin_size=as.numeric(coverage.resolution), extend=as.numeric(smoothing.extension))
	
	# Pattern Positions
	control<-MEDIPS.getPositions(data=control, pattern="CG")
	treatment<-MEDIPS.getPositions(data=treatment, pattern="CG")
	
	# Coupling vector
	control<-MEDIPS.couplingVector(data=control, fragmentLength=as.numeric(fragment.length), func="count")
	treatment<-MEDIPS.couplingVector(data=treatment, fragmentLength=as.numeric(fragment.length), func="count")
	
	# Calibration curve
	control<-MEDIPS.calibrationCurve(data=control)
	treatment<-MEDIPS.calibrationCurve(data=treatment)
	
	# Normalization
	control<-MEDIPS.normalize(data=control)
	treatment<-MEDIPS.normalize(data=treatment)
	
	# Comparing experiments
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
		rois<-rois[rois[,1] %in% unique(c(control@genome_chr), treatment@genome_chr),]
		write.table(rois, "rois.tsv", sep="\t", quote=F, col.names=F, row.names=F)
		diff.meth<-MEDIPS.methylProfiling(data1=control, data2=treatment, select=2, frame_size=400, ROI_file="rois.tsv")
		write.table(diff.meth, "output.tsv", sep="\t", quote=F, col.names=T, row.names=F)
	}
	
	if(promoters.only=="no") {
		diff.meth<-MEDIPS.methylProfiling(data1=control, data2=treatment, select=2, frame_size=400)
		write.table(diff.meth, file="output.tsv", sep="\t", quote=F, col.names=T, row.names=F)
	}
	
	# Do not delete the following curly bracket - it ends the if-clause for analysis and replication types
}





if( (analtype==1 | analtype==2) & repltype==2) {
	library(MEDIPS)
	
	# Assign the files to control and treatment based on the group column from phenodata
	cfiles<-as.character(phenodata$sample[phenodata$group==head(sort(phenodata$group), n=1)])
	tfiles<-as.character(phenodata$sample[phenodata$group==tail(sort(phenodata$group), n=1)])
	
	# Combining the replicates
	for(i in 1:length(cfiles)) {
		d<-read.table(cfiles, header=F, sep="\t")
		write.table(d, "group1.tsv", col.names=F, row.names=F, sep="\t", quote=F, append=TRUE)
	}
	for(i in 1:length(cfiles)) {
		d<-read.table(tfiles, header=F, sep="\t")
		write.table(d, "group2.tsv", col.names=F, row.names=F, sep="\t", quote=F, append=TRUE)
	}
	
	# Reads the data
	control<-MEDIPS.readAlignedSequences(BSgenome=genome, file="group1.tsv")
	treatment<-MEDIPS.readAlignedSequences(BSgenome=genome, file="group2.tsv")
	
	# Creating the genome vector
	control<-MEDIPS.genomeVector(data=control, bin_size=as.numeric(coverage.resolution), extend=as.numeric(smoothing.extension))
	treatment<-MEDIPS.genomeVector(data=treatment, bin_size=as.numeric(coverage.resolution), extend=as.numeric(smoothing.extension))
	
	# Pattern Positions
	control<-MEDIPS.getPositions(data=control, pattern="CG")
	treatment<-MEDIPS.getPositions(data=treatment, pattern="CG")
	
	# Coupling vector
	control<-MEDIPS.couplingVector(data=control, fragmentLength=as.numeric(fragment.length), func="count")
	treatment<-MEDIPS.couplingVector(data=treatment, fragmentLength=as.numeric(fragment.length), func="count")
	
	# Calibration curve
	control<-MEDIPS.calibrationCurve(data=control)
	treatment<-MEDIPS.calibrationCurve(data=treatment)
	
	# Normalization
	control<-MEDIPS.normalize(data=control)
	treatment<-MEDIPS.normalize(data=treatment)
	
	# Comparing experiments
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
		rois<-rois[rois[,1] %in% unique(c(control@genome_chr), treatment@genome_chr),]
		write.table(rois, "rois.tsv", sep="\t", quote=F, col.names=F, row.names=F)
		diff.meth<-MEDIPS.methylProfiling(data1=control, data2=treatment, select=2, frame_size=400, ROI_file="rois.tsv")
		write.table(diff.meth, "output.tsv", sep="\t", quote=F, col.names=T, row.names=F)
	}
	
	if(promoters.only=="no") {
		diff.meth<-MEDIPS.methylProfiling(data1=control, data2=treatment, select=2, frame_size=400)
		write.table(diff.meth, file="output.tsv", sep="\t", quote=F, col.names=T, row.names=F)
	}
	
	# Do not delete the following curly bracket - it ends the if-clause for analysis and replication types
}


