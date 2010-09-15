# TOOL "Statistics" / Test-multiple-file-output.R: "Takes a table of data, extracts the first 3 rows and returns them as individual 1-row tables" (Takes a table of data, extracts the first 3 rows and returns them as individual 1-row tables)
# INPUT table.tsv: "Data file" TYPE GENERIC
# OUTPUT result_{...}.tsv: "testing multiple file output"
# OUTPUT table_1.tsv: "also one at a time, first row"
# OUTPUT table_2.tsv: "also one at a time, second row"
# OUTPUT table_3.tsv: "also one at a time, thirs row"

#####################################################
#                                                   #
# MG, 24.6.2010                                     #
#                                                   #
# For debugging purposes only                       #
#                                                   #
#####################################################

# Read in data and convert to BED format
data_table <- read.table (file="table.tsv", sep="\t", header=T, row.names=1)

# Extract the first three rows into individual variables
first_row <- data_table[1,]
second_row <- data_table[2,]
third_row <- data_table[3,]

# Create the multiple files output
write.table(first_row, "result_1.tsv", sep="\t", row.names=T, col.names=T, quote=F)
write.table(second_row, "result_2.tsv", sep="\t", row.names=T, col.names=T, quote=F)
write.table(third_row, "result_3.tsv", sep="\t", row.names=T, col.names=T, quote=F)

# Create the individual files output
write.table(first_row, "table_1.tsv", sep="\t", row.names=T, col.names=T, quote=F)
write.table(second_row, "table_2.tsv", sep="\t", row.names=T, col.names=T, quote=F)
write.table(third_row, "table_3.tsv", sep="\t", row.names=T, col.names=T, quote=F)

# The end