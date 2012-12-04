# TOOL head.R: "Get first lines" (Get first lines of a dataset. Similar to Unix command head.)
# INPUT input-file: input-file TYPE GENERIC 
# OUTPUT first-lines.tsv: first-lines.tsv 
# PARAMETER count:  TYPE INTEGER FROM 0 (Number of lines to get)

system(paste("head -n", count, "input-file > first-lines.tsv"))