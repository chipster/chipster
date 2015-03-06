# TOOL output_test.R: "Output test" (Output test.)
# INPUT bedfile.bed: "BED file with \"description\" column to be extracted" TYPE GENERIC 
# INPUT sequences.fa: "Your fasta file" TYPE GENERIC
# OUTPUT my-output{...}.tsv: "Multiple files" 
# PARAMETER column: Column TYPE COLUMN_SEL (Data \"column\" from first input)
# PARAMETER num: Num TYPE DECIMAL FROM 0 TO 10 DEFAULT 1 (Some \(test\) number)

# MK 18.11.2013

# check for valid accession

rnd_num <- round((runif(1, 1.0, 10)),0)

for(i in 1:rnd_num) {
  write.table(rnd_num, file=paste("some-other-output", i, ".tsv", sep=""))
}

for(i in 1:rnd_num) {
  write.table(rnd_num, file=paste("my-output", i, "blaah.tsv", sep=""))
}

# EOF
