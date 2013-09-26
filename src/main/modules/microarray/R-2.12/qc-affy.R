# TOOL qc-affy.R: "Affymetrix basic" (Affymetrix quality control for RNA degradation and general quality parameters, such as scaling factor. This tool should be run on RAW data, i.e., CEL-files.)
# INPUT microarray{...}.cel: microarray{...}.cel TYPE AFFY 
# OUTPUT RNA-degradation-plot.pdf: RNA-degradation-plot.pdf 
# OUTPUT simpleaffy-plot.pdf: simpleaffy-plot.pdf 
# OUTPUT spike-in-plot.pdf: spike-in-plot.pdf 
# PARAMETER image.width: image.width TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height: image.height TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Affymetrix quality control
# JTT:  9.6.2006
# MG: 31.12.2009
# MK: 20.05.2013

#image.width<-600
#image.height<-600

# Loading the libraries
library(affy)
library(simpleaffy)

# Renaming variables
w<-image.width
h<-image.height

# Reading in data
dat<-ReadAffy()

# Setting up coloring and line types for plotting
library (RColorBrewer)	
number_samples <- dim(pData(dat))[1]
number_colors <- number_samples
number_columns <- ceiling(number_samples/16)
color_range <- brewer.pal (n=7, name="Set1")
color_range <- c(color_range,"black")
if (number_samples >= 8) {
	color_array <- rep(color_range, 100)
}
if (number_samples < 8) {
	color_array <- rep(color_range[1:number_samples], 100)
}
line_type <- c(rep(1,8),rep(2,8),rep(3,8),rep(4,8),rep(5,8),rep(6,8))
line_type <- rep(line_type, 100)
plot_symbol <- 1:8
if (number_samples >= 8) {
	plot_symbol <- rep(plot_symbol, 100)
}
if (number_samples < 8) {
	plot_symbol <- rep(plot_symbol[1:number_samples], 100)
}

if(length(which(is.na(rownames(mm(dat))))) / length(rownames(mm(dat))) > 0) {
	stop("CHIPSTER-NOTE: The simpleaffy quality analysis tool has not been designed for PM-only arrays. Please use another quality assessment method");
}

aqc<-try(qc(dat))
if(class(aqc)=="try-error") {
	stop("CHIPSTER-NOTE: Your array type is not supported by the simpleaffy quality analysis tool. Please use another quality assessment method");
}
	
# Calculating quality control values
#aqc<-qc(dat)

# Plotting the QC-values
pdf(file="simpleaffy-plot.pdf", width=w/72, height=h/72)
plot(aqc)
dev.off()

# Checking the RNA degradation
deg<-AffyRNAdeg(dat)

# Define modified function for plotting RNA degradation plot
plot_RNA_deg <-function (rna.deg.obj, transform = "shift.scale", cols = NULL,
		line_type = 1, ...) 
{
	if (!is.element(transform, c("shift.scale", "shift.only", 
					"neither"))) 
		stop("Tranform must be 'shift.scale','shift.only', or 'neither'")
	mns <- rna.deg.obj$means.by.number
	if (is.null(cols)) 
		cols = rep(4, dim(mns)[1])
	ylab = "Mean Intensity"
	if (transform == "shift.scale") {
		sds <- rna.deg.obj$ses
		mn <- mns[, 1]
		mns <- sweep(mns, 1, mn)
		mns <- mns/(sds)
		mns <- sweep(mns, 1, 1:(dim(mns)[1]), "+")
		ylab <- paste(ylab, ": shifted and scaled")
	}
	else if (transform == "shift.only") {
		mn <- mns[, 1]
		mns <- sweep(mns, 1, mn)
		mns <- sweep(mns, 1, 1:(dim(mns)[1]), "+")
		ylab <- paste(ylab, ": shifted")
	}
	plot(-2, -1, pch = "", xlim = range(-1, (dim(mns)[2])), ylim = range(min(as.vector(mns)) - 
							1, max(as.vector(mns)) + 1), xlab = "5' <-----> 3'\n Probe Number ", 
			ylab = ylab, axes = FALSE, main = "RNA degradation plot", 
			...)
	axis(1)
	axis(2)
	for (i in 1:dim(mns)[1]) lines(0:((dim(mns)[2] - 1)), mns[i, 
				], col = cols[i], lty=line_type[i])
}
# Saving the degradation result into a file
pdf(file="RNA-degradation-plot.pdf", width=w/72, height=h/72)
plot_RNA_deg(deg, col=color_array, line_type=line_type)
legend(legend=sampleNames(dat), x="topleft", lty=line_type, cex=0.75, 
	col=color_array, ncol=number_columns)
dev.off()

# Assess and plot spike-in performance
par(cex=0.75, cex.main=1.5, mar=c(5,5,4,2))
concentration <- log(c(1.5, 5, 25, 100))
x_values <- array(concentration, c(4, number_samples))
x_values <- t(x_values)
y_values <- spikeInProbes(aqc)
pdf(file="spike-in-plot.pdf", width=w/72, height=h/72)
plot(x_values, y_values, col=color_array, main="Spike-in performance",
		xlab="log10 (concentration in pM)", ylab="log2 (expression)",
		pch=plot_symbol, type="p")
legend(legend=sampleNames(dat), x="topleft",
		col=color_array, lty=line_type, ncol=number_columns, cex=0.75,
		pch=plot_symbol)

for (loop_count in 1:length(dat))
{
	y_values <- spikeInProbes(aqc) [loop_count,]
	lm_spike <- lm(y_values~concentration)
	slope <- coef(lm_spike) [2]
	intercept <- coef(lm_spike) [1]
	abline(intercept, slope, col=color_array[loop_count], lty=line_type[loop_count])
}
dev.off()
