# ANALYSIS "Normalisation"/"Median normalisation (cDNA)" (Calculates the mean of the log ratios for one microarray. Produces the centered data by subtracting this mean from the log ratio of every gene.)
# INPUT CDNA microarray.txt OUTPUT results.txt
# PARAMETER Percentile INTEGER FROM 0 TO 100 DEFAULT 50 (the percentile to center with)
# PARAMETER Percentile2 PERCENT DEFAULT 50 (the percentile to center with)
# PARAMETER "Remove negative values" [yes, no] DEFAULT yes (remove negative values)
# PARAMETER "Replacing value" DECIMAL DEFAULT 0.01 (when replacing negative values, replace with this value)

t <- read.table("microarray.txt", header=T, sep="\t");
ratio <- (t['CH1I_MEAN'] - t['CH1B_MEDIAN'])/(t['CH2I_MEAN'] - t['CH2B_MEDIAN'])
norm.ratio <- ratio - median(ratio[,1])
result <- cbind(t['Accession'], norm.ratio)
write.table(result, file="results.txt", quote=FALSE, col.names=c('PROBESET_ID', 'EXPRS'), row.names=FALSE)
