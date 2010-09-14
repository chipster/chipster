# ANALYSIS "aCGH tools (beta testing)"/"Calculate descriptive statistics for a column" (Calculates basic descriptive statistics for a given column.)
# INPUT GENERIC normalized.tsv
# OUTPUT statistics.tsv
# PARAMETER column COLUMN_SEL (The column to calculate statistics for.)

# calculate-descriptive-statistics-for-col.R
# Ilari Scheinin <firstname.lastname@helsinki.fi>
# 2010-09-05

# load input
dat <- read.table('normalized.tsv', header=TRUE, sep='\t', as.is=TRUE)

# calculate statistics and round off extra digits
output <- data.frame(column=column, minimum=min(dat[,column]), median=median(dat[,column]), mean=mean(dat[,column]), maximum=max(dat[,column]), standard.deviation=sd(dat[,column]), variance=var(dat[,column]))
output[,-1] <- round(output[,-1], digits=3)

# write output
write.table(output, file="statistics.tsv", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

# EOF