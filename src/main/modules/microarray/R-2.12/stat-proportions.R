# TOOL stat-proportions.R: "Test proportions" (Tests for comparing proportions in dichotomized data.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT prop-test.tsv: prop-test.tsv 
# PARAMETER column: column TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER p.value.adjustment.method: p.value.adjustment.method TYPE [none: none, Bonferroni: Bonferroni, Holm: Holm, Hochberg: Hochberg, BH: BH, BY: BY, Storey-Q: Storey-Q] DEFAULT BH (Multiple testing correction method)
# PARAMETER p.value.threshold: p.value.threshold TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (P-value cut-off for significant results)

# MG, 29.9.2011

# Optional parameter
# PARAMETER test [Fisher, Chisquare, hypergeometric] DEFAULT Chisquare (Test type)

# Loads the libraries
library(multtest)
library(qvalue)

# Loads the normalized data
file <- c("normalized.tsv")
dat <- read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls <- dat[,grep("flag", names(dat))]
dat_2 <- dat[,grep("chip", names(dat))]

# Test needs a parameter "groups" that specifies the grouping of the samples
phenodata <- read.table("phenodata.tsv", header=T, sep="\t")
groups <- phenodata[,pmatch(column,colnames(phenodata))]
group_levels <- (unique(groups))

# Sanity checks
num_groups <- length(unique(groups))
if(num_groups<=1) {
	stop("CHIPSTER NOTE: You need to have at least two groups to run this analysis")
}

# Calculate number of 0 and 1 occurrencies for each group of samples
dat_0 <- matrix(nrow=nrow(dat_2), ncol=length(group_levels), NA)
dat_1 <- matrix(nrow=nrow(dat_2), ncol=length(group_levels), NA)
for(group_count in 1:length(group_levels)) {
	dat_0[,group_count] <- rowSums(data.frame(dat_2[,which(groups==group_levels[group_count])]))
	dat_1[,group_count] <- ncol(dat_2)-dat_0[,group_count] 
}

# Perform test of equal proportions per row of data
results <- matrix (nrow=nrow(dat_2), ncol=2)
for (row_number in 1:nrow(dat_2)) {
	counts_table <- matrix(nrow=length(group_levels),ncol=2,data=c(dat_0[row_number,],dat_1[row_number,]), byrow=FALSE)
	results[row_number,1] <- prop.test(counts_table)$statistic
	results[row_number,2] <- prop.test(counts_table)$p.value
}

# If NAN result replace with statistic 0 and p-value 1
results[is.nan(results[,1]),1] <- 0
results[is.nan(results[,2]),2] <- 1

# Extract raw p-values
raw_p <- results[,2]

# Perform multiple testing correction
if (p.value.adjustment.method != "none") {
	if (p.value.adjustment.method == "Storey-Q") {
		adj_p <- qvalue(raw_p)
		adj_p2 <- adj_p$qvalues
	} else {
		adj_p <- mt.rawp2adjp(raw_p, proc=as.character(p.value.adjustment.method))
		adj_p2 <- adj_p$adjp[order(adj_p$index),][,2]
	}
	adj_p2df <- as.data.frame(adj_p2)
	results <- cbind(results, adj_p2df)
	colnames(results) <- c("chi_square","p_value", "adjusted_p_value")
} else {
	colnames(results) <- c("chi_square","p_value")
}

# Extract only significant results
dat_out <- data.frame(dat,results)
if (p.value.adjustment.method != "none") {
	dat_out <- dat_out[dat_out$adjusted_p_value <= p.value.threshold,]   
} else {
	dat_out <- dat_out[dat_out$p_value <= p.value.threshold,]   
}

# Writing out data
write.table(dat_out, file="prop-test.tsv", sep="\t", row.names=T, col.names=T, quote=F)
