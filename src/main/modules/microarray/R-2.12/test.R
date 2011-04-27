# TOOL test.R: "Parameter demo" (Demonstration of all parameter types.)
# INPUT microarray.txt: microarray.txt TYPE CDNA 
# OUTPUT results.txt: results.txt 
# PARAMETER value1: value1 TYPE INTEGER FROM 0 TO 200 DEFAULT 10 (the first value of the result set)
# PARAMETER value2: value2 TYPE DECIMAL FROM 0 TO 200 DEFAULT 20 (the second value of the result set)
# PARAMETER value3: value3 TYPE DECIMAL FROM 0 TO 200 DEFAULT 30.2 (the third value of the result set)
# PARAMETER cutoff: cutoff TYPE PERCENT DEFAULT 34 (cut-off value for selecting data points)
# PARAMETER genename: genename TYPE STRING DEFAULT at_something (which gene we are interested in)
# PARAMETER key: key TYPE COLUMN_SEL (which column we use as a key)
# PARAMETER metakey: metakey TYPE METACOLUMN_SEL (which metadata column we use as a metakey)
# PARAMETER input: input TYPE INPUT_SEL (which input we want)
# PARAMETER method: method TYPE [linear: linear, logarithmic: logarithmic, exponential: exponential] DEFAULT logarithmic (which method to apply)

write.table(c(value1, value2, value3), file="results.txt", quote=FALSE, col.names=c('EXPRS'), row.names=FALSE)
