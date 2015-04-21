# TOOL join_tsv_files.R: "Join tables" (Joins TSV format tables into one file. The input files must be in same format.)
# INPUT table{...}.tsv: "TSV files" TYPE GENERIC
# OUTPUT merged.tsv 
# PARAMETER header: "Tables contain a header row" TYPE [yes, no] DEFAULT yes (Do the TSV files contain a header row?)
# PARAMETER remove.duplicates: "Remove duplicate rows" TYPE [yes, no] DEFAULT no (Remove exact duplicates of rows.)

source(file.path(chipster.common.path, "zip-utils.R"))
files <- list.files(pattern="*.tsv", full.names=FALSE, recursive=FALSE)
for (i in 1:length(files)) {
	unzipIfGZipFile(files[i])
}

if (header == "yes"){
	system("head -1 table001.tsv > merged_tmp.tsv")
	system("for file in $( ls table*.tsv ) ; do tail -n+2 ${file} ; done >> merged_tmp.tsv")
}else{
	system("cat *.tsv | grep -v ^$ > merged_tmp.tsv")
}

if (remove.duplicates == "yes"){
	system("cat merged_tmp.tsv|sort|uniq > merged_tmp2.tsv")
	system("mv merged_tmp2.tsv merged_tmp.tsv")
}

# Remove empty lines
system("tr -d '\r' < merged_tmp.tsv > t ; mv t merged_tmp.tsv")
system("grep -v ^$ merged_tmp.tsv > merged.tsv")





