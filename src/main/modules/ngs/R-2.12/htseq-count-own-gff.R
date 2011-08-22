# TOOL htseq-count-own-gtf.R: "Count reads with HTSeq using your own feature files" (Given a BAM file, this tool calculates how many reads map to each gene or exon.
# This tool allows you count reads against your own list of genes or exons, which you have to supply in GTF format.
# If you would like to align reads against public GTF files supplied by Chipster, please use the tool \"Count reads with HTSeq\".)
# INPUT alignment.bam: "BAM alignment file" TYPE GENERIC
# INPUT features.gtf: "GTF feature file" TYPE GENERIC
# OUTPUT htseq-counts.txt
# OUTPUT OPTIONAL htseq-count-info.txt
# PARAMETER paired: "Alignment file contains paired-end data" TYPE [yes, no] DEFAULT no (Does the alignment data contain paired end or single end reads?)
# PARAMETER mode: "Mode to handle reads overlapping more than one feature" TYPE [union, intersection-strict, intersection-nonempty] DEFAULT union (How to deal with reads that overlap more than one gene or exon?)
# PARAMETER stranded: "Was a strand-specific protocol used for RNA-seq" TYPE [yes, no, reverse] DEFAULT no (If you select no, a read is considered overlapping with a feature regardless of whether it is mapped to the same or the opposite strand as the feature. If you select yes, the read has to be mapped to the same strand as the feature. You have to say no, if yours was not made with a strand-specific RNA-seq protocol, because otherwise half your reads will be lost.)
# PARAMETER minaqual: "Minimum alignment quality" TYPE INTEGER FROM 0 TO 100 DEFAULT 0 (Skip all reads with alignment quality lower than the given minimum value.)
# PARAMETER feature.type: "Which feature type to count" TYPE [exon, CDS] DEFAULT exon (Which feature type to use, all features of other type are ignored.)
# PARAMETER id.attribute: "Which feature ID to use" TYPE [gene_id, transcript_id, gene_name, transcript_name, protein_name] DEFAULT gene_id (GFF attribute to be used as feature ID. Several GFF lines with the same feature ID will be considered as parts of the same feature. The feature ID is used to identity the counts in the output table.)

# EK 22.8.2011

# bash wrapping
python.path <- paste(sep="", "PYTHONPATH=", file.path(chipster.tools.path, "lib", "python2.6", "site-packages"), ":$PYTHONPATH")
command.start <- paste("bash -c '", python.path, ";")
command.end <- "'"

# sort bam if needed
samtools.binary <- file.path(chipster.tools.path, "samtools", "samtools")
samtools.sort <- ifelse(paired == "yes", paste(samtools.binary, "sort -on alignment.bam sorted-by-name"), "cat alignment.bam")

# convert bam to sam
samtools.view <- paste(samtools.binary, "view -")

# htseq-count
htseq.binary <- c(file.path(chipster.tools.path, "htseq", "htseq-count"))
htseq <- paste(htseq.binary, "-q -m", mode, "-s", stranded, "-a", minaqual, "-t", feature.type, "-i", id.attribute, "-", "features.gtf > htseq-counts-out.txt")

# run
command <- paste(command.start, samtools.sort, " | ", samtools.view, " | ", htseq, command.end)
system(command)

# separate result file
system("head -n -5 htseq-counts-out.txt > htseq-counts.txt")
system("tail -n 5 htseq-counts-out.txt > htseq-count-info.txt")
