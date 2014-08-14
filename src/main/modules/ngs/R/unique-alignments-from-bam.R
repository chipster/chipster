# TOOL unique-alignments-from-bam.R: "Retrieve unique alignments from BAM" (Retrieve unique alignments from BAM. Input BAM needs to be sorted. This tool is based on the SAMtools package.)
# INPUT alignment.bam: "BAM alignment file" TYPE GENERIC
# OUTPUT OPTIONAL unique_alignments.bam
# PARAMETER paired: "Does the alignment file contain paired-end data" TYPE [yes, no] DEFAULT no (Does the alignment data contain paired end or single end reads?)

# 2014.06.02 AMS 

# convert bam to sam, sort bam if the data is paired-end
samtools.binary <- file.path(chipster.tools.path, "samtools", "samtools")
filtersamreads.binary <- file.path(chipster.tools.path, "picard-tools", "FilterSamReads.jar")

if(paired == "yes") {
	#List of IDs having mapping quality >= 4 and having proper pairing (0x0002). Uniq command make sure that both pairs must fulfill filtering criteria
	keep_id.command <- paste(samtools.binary, "view -q4 -F4 -f2 alignment.bam | cut -f1 | sort | uniq -d > alignment_id_keep.txt")
	system(keep_id.command)
	filtersamreads.command <- paste("java -Xmx4096m -jar", filtersamreads.binary, "I=alignment.bam FILTER=includeReadList RLF=alignment_id_keep.txt O=alignment_keep.bam VALIDATION_STRINGENCY=LENIENT")
	system(filtersamreads.command) 
	samtools.view <- paste(samtools.binary, "view -h alignment_keep.bam | perl -p -e 's/\n$/\tNH:i:1\n/g'")
} else if(paired == "no") {
	samtools.view <- paste(samtools.binary, "view -h -q4 -F4 alignment.bam | perl -p -e 's/\n$/\tNH:i:1\n/g'")
}
# run
command <- paste(samtools.view, "|", samtools.binary, "view -bS - > unique_alignments.bam")
system(command)
