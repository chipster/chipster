# ANALYSIS "Utilities"/"Parameter demo" (Demonstration of all parameter types.)
# INPUT CDNA microarray.txt OUTPUT results.txt
# PARAMETER value1 INTEGER FROM 0 TO 200 DEFAULT 10 (the first value of the result set)
# PARAMETER value2 DECIMAL FROM 0 TO 200 DEFAULT 20 (the second value of the result set)
# PARAMETER value3 DECIMAL FROM 0 TO 200 DEFAULT 30.2 (the third value of the result set)
# PARAMETER cutoff PERCENT DEFAULT 34 (cut-off value for selecting data points)
# PARAMETER genename STRING DEFAULT at_something (which gene we are interested in)
# PARAMETER key COLUMN_SEL (which column we use as a key)
# PARAMETER metakey METACOLUMN_SEL (which metadata column we use as a metakey)
# PARAMETER input INPUT_SEL (which input we want)
# PARAMETER method [linear, logarithmic, exponential] DEFAULT logarithmic (which method to apply)

write.table(c(value1, value2, value3), file="results.txt", quote=FALSE, col.names=c('EXPRS'), row.names=FALSE)
