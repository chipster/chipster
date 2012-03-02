# TOOL dichotomize-continious-data.R: "Dichotomize continious data" (Dichotomize the data into values 0 and 1 based on a single cut-off value.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT dichotomized.tsv: dichotomized.tsv
# PARAMETER what.mode: "what mode" TYPE [absolute: absolute, relative: relative] DEFAULT relative (Determines whether to apply an absolute or relative level cut-off. If relative level is used the cut-off represents a percentile of the data distribution.)
# PARAMETER cut.off: cut-off TYPE DECIMAL FROM 0 TO 1000000 DEFAULT 50 (The absolute or relative percentage cut-off defining the data to be assigned the value of 1. The remainder of the data will be assigned the value of 0.)

# MG, 23.9.2011

# Parameter settings (default) for testing purposes
# what.mode <- "absolute"
# cut.off <- 0.2

# Loading the libraries
library(genefilter)

# Loads the normalized data
file <- c("normalized.tsv")
dat <- read.table(file, header=T, sep="\t", row.names=1)

# Extract the data and keep the rest for a while
dat2 <-dat[,grep("chip", names(dat))]
dat2 <- as.matrix(dat2)
the.rest <- dat[,-grep("chip",names(dat))]

# Sanity checks
if (what.mode == "relative" && cut.off > 100) {
	stop("CHIPSTER-NOTE: In relative mode the cut-off cannot by definition be outside the range 0-100")
}
if (ncol(dat2) == 1) {
	stop("CHIPSTER-NOTE: You need to have at least two chips to apply filtering by genes!")
}

# Figure out what value the cut-off percentile corresponds to and dichotomize
dat3 <- as.vector(dat2)
if (what.mode == "relative") {
	cutoff.value <- quantile(dat3, probs=cut.off/100)
} else {
	cutoff.value <- cut.off
}
dat4 <- ifelse((dat3>cutoff.value),1,0) 
dat5 <- matrix(dat4, ncol=dim(dat2)[2], nrow=dim(dat2)[1])
colnames(dat5) <- colnames(dat2)
rownames(dat5) <- rownames(dat2)


# Join the transformed data and the rest
if (ncol(the.rest)==0) {
	output.table <- dat5
} else {
	output.table <- cbind (as.data.frame(dat5), the.rest)
}

# Saving the results
write.table(data.frame(output.table), file=("dichotomized.tsv"), sep="\t", row.names=T, col.names=T, quote=F)

