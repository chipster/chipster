# TOOL extract-features-from-BED.R: "Extract features from BED" (This tool extracts the description column, usually gene identifiers, removes duplicate entries and creates a table containing only the extracted features. The output is compatible with downstream tools.)
# INPUT bedfile.bed "BED file with description column to extract" TYPE GENERIC 
# OUTPUT extracted-features.tsv: "Table listing the unique features found in the description column of the bed file." 

#######################################################
#                                                     #
# MG, 1.11.2011                                       #
#                                                     #
# Tool that extracts a column with genomic features   #
# from a BED file in order to make it suitable to     #
# downstream analysis.                                #
#                                                     #
#######################################################

file <- "bedfile.bed"
dat2 <- read.table(file=file, header=TRUE, sep="\t")
list_ids <- dat[,4]
unique_ids <- unique(list_ids)
output_table <- as.data.frame(unique_ids)
rownames(output_table) <- unique_ids
colnames(output_table) <- "id"
write.table(output_table, file="extracted-features.tsv", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

# EOF
