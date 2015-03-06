# TOOL qc-affy-exon-rle-nuse.R: "Affymetrix exon/gene arrays - using RLE and NUSE" (Affymetrix quality control using NUSE and RLE. This tool should be run on RAW data, i.e., CEL-files, for exon arrays. The chip type and the summary feature have to be specified.)
# INPUT microarray{...}.cel: microarray{...}.cel TYPE AFFY 
# OUTPUT rle-plot.png: rle-plot.png 
# OUTPUT nuse-plot.png: nuse-plot.png 
# PARAMETER image.width: "Image width" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height: "Image height" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)
# PARAMETER chiptype: "Chiptype" TYPE [empty: empty, human-exon: human-exon, mouse-exon: mouse-exon, rat-exon: rat-exon, human-1.0-ST: HuGene-1.0-ST, human-1.1-ST: HuGene-1.1-ST, human-2.0-ST: HuGene-2.0-ST, human-2.1-ST: HuGene-2.1-ST, human-hta20: human-hta20, human-prime: human-prime, mouse-1.0-ST: MoGene-1.0-ST, mouse-1.1-ST: MoGene-1.1-ST, mouse-2.0-ST: MoGene-2.0-ST, mouse-2.1-ST: MoGene-2.1-ST, rat-1.0-ST: RaGene-1.0-ST, rat-1.1-ST: RaGene-1.1-ST, rat-2.0-ST: RaGene-2.0-ST, rat-2.1-ST: RaGene-2.1-ST, zebra_fish-1.0-ST: zebra_fish_gene-1.0-ST, zebra_fish-1.1-ST: zebra_fish_gene-1.1-ST, arabidopsis-1.0-ST-entrez: arabidopsis_gene-1.0-ST-entrez, arabidopsis-1.1-ST-entrez: arabidopsis_gene-1.1-ST-entrez, arabidopsis-1.0-ST-tair: arabidopsis_gene-1.0-ST-tair, arabidopsis-1.1-ST-tair: arabidopsis_gene-1.1-ST-tair, fly-1.0-ST: FlyGene-1.0-ST, fly-1.1-ST: FlyGene-1.1-ST, celegans-1.0-ST: celegans_gene-1.0-ST, celegans-1.1-ST: celegans_gene-1.1-ST, dog-1.0-ST: DogGene-1.0-ST, dog-1.1-ST: DogGene-1.1-ST, rice-1.0-ST-entrez: rice_gene-1.0-ST, rice-1.1-ST-entrez: rice_gene-1.1-ST, oligo: oligo] DEFAULT empty (The use of empty and oligo envokes the default behaviour of the AffyPLM and oligo packages, respectively. If default annotation packages are found, AffyPLM and oligo will use them)
# PARAMETER summary.feature: "Summary feature" TYPE [gene: gene, exon: exon] DEFAULT gene (Output summary type for Exon arrays)


# Affymetrix quality control
# MG 12.1.2010
# MK: 12.06.2013 added possibility to analyse custom chips
# MK: 20.02.2014 added support for oligo-package
# MK: 07.04.2014 added support for new gene chips

# Loading the libraries
library(affy)
library(affyPLM)
library(oligo)

# Renaming variables
w<-image.width
h<-image.height

headdetails <- affyio::read.celfile.header(as.character(list.celfiles()[1]))

# Reading in data
if(chiptype == "oligo") {
	if(length(grep("HuEx", headdetails$cdfName, ignore.case = T)) > 0) {
		data.raw <- read.celfiles(filenames=list.celfiles(), pkgname="pd.huex.1.0.st.v2")
	} else {
		data.raw <- read.celfiles(filenames=list.celfiles())		
	}
} else {
	#Try if cel-files can be read using ReadAffy
	dat <- try(ReadAffy())

	#If ReadAffy not working (i.e. gene st2, exon arrays), fake affy-batch object
	#Could be done commeting out line 81-85 from ReadAffy
	if(class(dat)=="try-error") {
		if(length(grep("HuEx", headdetails$cdfName, ignore.case = T)) > 0) {
			dat <- read.celfiles(filenames=list.celfiles(), pkgname="pd.huex.1.0.st.v2")
		} else {
			dat <- read.celfiles(filenames=list.celfiles())		
		}
		dim.intensity <- headdetails[[2]]

		class(dat) <- "AffyBatch"
		dat@nrow <- dim.intensity[2]
		dat@ncol <- dim.intensity[1]
	}
}

# Set up proper cdf package information
if(chiptype=="empty") {
	chiptype<-dat@annotation
	a <- try(library(paste(chiptype,"cdf",sep=""), character.only=T))
	if(class(a) == "try-error") {
		stop("You need to specify the chiptype. Please run the script again.")
	}
}
if(chiptype=="human-exon" & summary.feature=="exon") {
	dat@cdfName<-"huex10stv2hsensecdf"
	dat@annotation<-"huex10stv2hsensecdf"
}
if(chiptype=="mouse-exon" & summary.feature=="exon") {
	dat@cdfName<-"moex10stv1mmensecdf"
	dat@annotation<-"moex10stv1mmensecdf"
}
if(chiptype=="rat-exon" & summary.feature=="exon") {
	dat@cdfName<-"raex10stv1rnensecdf"
	dat@annotation<-"raex10stv1rnensecdf"
}

if(chiptype=="human-exon" & summary.feature=="gene") {
	dat@cdfName<-"huex10stv2hsentrezgcdf"
	dat@annotation<-"huex10stv2hsentrezg.db"
}
if(chiptype=="mouse-exon" & summary.feature=="gene") {
	dat@cdfName<-"moex10stv1mmentrezgcdf"
	dat@annotation<-"moex10stv1mmentrezg.db"
}
if(chiptype=="rat-exon" & summary.feature=="gene") {
	dat@cdfName<-"raex10stv1rnentrezgcdf"
	dat@annotation<-"raex10stv1rnentrezg.db"
}

if(chiptype=="human-1.0-ST") {
	dat@cdfName<-"hugene10sthsentrezgcdf"
	dat@annotation<-"hugene10sthsentrezgcdf"
}
if(chiptype=="human-1.1-ST") {
	dat@cdfName<-"hugene11sthsentrezgcdf"
	dat@annotation<-"hugene11sthsentrezgcdf"
}
if(chiptype=="human-2.0-ST") {
	dat@cdfName<-"hugene20sthsentrezgcdf"
	dat@annotation<-"hugene20sthsentrezgcdf"
}
if(chiptype=="human-2.1-ST") {
	dat@cdfName<-"hugene21sthsentrezgcdf"
	dat@annotation<-"hugene21sthsentrezgcdf"
}
if(chiptype=="human-hta20") {
	dat@cdfName<-"hta20hsentrezgcdf"
	dat@annotation<-"hta20hsentrezgcdf"
}
if(chiptype=="mouse-1.0-ST") {
	dat@cdfName<-"mogene10stmmentrezgcdf"
	dat@annotation<-"mogene10stmmentrezgcdf"
}
if(chiptype=="mouse-1.1-ST") {
	dat@cdfName<-"mogene11stmmentrezgcdf"
	dat@annotation<-"mogene11stmmentrezgcdf"
}
if(chiptype=="mouse-2.0-ST") {
	dat@cdfName<-"mogene20stmmentrezgcdf"
	dat@annotation<-"mogene20stmmentrezgcdf"
}
if(chiptype=="mouse-2.1-ST") {
	dat@cdfName<-"mogene21stmmentrezgcdf"
	dat@annotation<-"mogene21stmmentrezgcdf"
}
if(chiptype=="rat-1.0-ST") {
	dat@cdfName<-"ragene10strnentrezgcdf"
	dat@annotation<-"ragene10strnentrezgcdf"
}
if(chiptype=="rat-1.1-ST") {
	dat@cdfName<-"ragene11strnentrezgcdf"
	dat@annotation<-"ragene11strnentrezgcdf"
}
if(chiptype=="rat-2.0-ST") {
	dat@cdfName<-"ragene20strnentrezgcdf"
	dat@annotation<-"ragene20strnentrezgcdf"
}
if(chiptype=="rat-2.1-ST") {
	dat@cdfName<-"ragene21strnentrezgcdf"
	dat@annotation<-"ragene21strnentrezgcdf"
}
if(chiptype=="zebra_fish-1.0-ST") {
	dat@cdfName<-"zebgene10stdrentrezgcdf"
	dat@annotation<-"zebgene10stdrentrezgcdf"
}
if(chiptype=="zebra_fish-1.1-ST") {
	dat@cdfName<-"zebgene11stdrentrezgcdf"
	dat@annotation<-"zebgene11stdrentrezgcdf"
}
if(chiptype=="arabidopsis-1.0-ST-entrez") {
	dat@cdfName<-"aragene10statentrezgcdf"
	dat@annotation<-"aragene10statentrezgcdf"
}
if(chiptype=="arabidopsis-1.1-ST-entrez") {
	dat@cdfName<-"aragene11statentrezgcdf"
	dat@annotation<-"aragene11statentrezgcdf"
}
if(chiptype=="arabidopsis-1.0-ST-tair") {
	dat@cdfName<-"aragene10stattairgcdf"
	dat@annotation<-"aragene10stattairgcdf"
}
if(chiptype=="arabidopsis-1.1-ST-tair") {
	dat@cdfName<-"aragene11stattairgcdf"
	dat@annotation<-"aragene11stattairgcdf"
}

if(chiptype=="human-prime") {
	dat@cdfName<-"primeviewhsentrezgcdf"
	dat@annotation<-"primeviewhsentrezgcdf"
}
if(chiptype=="fly-1.0-ST") {
	dat@cdfName<-"drogene10stdmentrezgcdf"
	dat@annotation<-"drogene10stdmentrezgcdf"
}
if(chiptype=="fly-1.1-ST") {
	dat@cdfName<-"drogene11stdmentrezgcdf"
	dat@annotation<-"drogene11stdmentrezgcdf"
}
if(chiptype=="celegans-1.0-ST") {
	dat@cdfName<-"elegene10stceentrezgcdf"
	dat@annotation<-"elegene10stceentrezgcdf"
}
if(chiptype=="celegans-1.1-ST") {
	dat@cdfName<-"elegene11stceentrezgcdf"
	dat@annotation<-"elegene11stceentrezgcdf"

}
if(chiptype=="dog-1.0-ST") {
	dat@cdfName<-"cangene10stcfentrezgcdf"
	dat@annotation<-"cangene10stcfentrezgcdf"
}
if(chiptype=="dog-1.1-ST") {
	dat@cdfName<-"cangene11stcfentrezgcdf"
	dat@annotation<-"cangene11stcfentrezgcdf"
}
if(chiptype=="rice-1.0-ST-entrez") {
	dat@cdfName<-"ricegene10stosentrezgcdf"
	dat@annotation<-"ricegene10stosentrezgcdf"
}
if(chiptype=="rice-1.1-ST-entrez") {
	dat@cdfName<-"ricegene11stosentrezgcdf"
	dat@annotation<-"ricegene11stosentrezgcdf"
}

# Calculating quality control values
if(chiptype == "oligo") {
	aqc <- fitProbeLevelModel(data.raw, target="core")

	# Plotting the QC-values
	par(mar=c(7, 4, 4, 2) + 0.1)
	title <- paste("RLE\n(", summary.feature, " level)", sep="")
	bitmap(file="rle-plot.png", width=w/72, height=h/72)
	RLE(aqc, main=title, las=2)
	dev.off()

	par(mar=c(7, 4, 4, 2) + 0.1)
	title <- paste("NUSE\n(", summary.feature, " level)", sep="")
	bitmap(file="nuse-plot.png", width=w/72, height=h/72)
	NUSE(aqc, main=title, las=2)
	dev.off()
} else {
	aqc<-affyPLM::fitPLM(dat)

	# Plotting the QC-values
	par(mar=c(7, 4, 4, 2) + 0.1)
	title <- paste("RLE\n(", summary.feature, " level)", sep="")
	bitmap(file="rle-plot.png", width=w/72, height=h/72)
	Mbox(aqc, main=title, las=2)
	dev.off()

	par(mar=c(7, 4, 4, 2) + 0.1)
	title <- paste("NUSE\n(", summary.feature, " level)", sep="")
	bitmap(file="nuse-plot.png", width=w/72, height=h/72)
	boxplot(aqc, main=title, las=2)
	dev.off()
}
