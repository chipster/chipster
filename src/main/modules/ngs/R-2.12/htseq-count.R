# TOOL htseq-count.R: "Map aligned reads to genes with HTSeq" (Calculates how many reads in a BAM file map to each gene. If you would like to map reads against your own GTF files, please use the tool \"Map aligned reads to genes with HTSeq using own GTF\". This tool is based on the HTSeq package. In order to use the output in edgeR, you need to select all samples and run the tool \"Utilities - Define NGS experiment\".)
# INPUT alignment.bam: "BAM alignment file" TYPE GENERIC
# OUTPUT htseq-counts.tsv 
# OUTPUT OPTIONAL htseq-count-info.txt
# PARAMETER paired: "Does the alignment file contain paired-end data" TYPE [yes, no] DEFAULT no (Does the alignment data contain paired end or single end reads?)
# PARAMETER stranded: "Was the data produced with a strand-specific RNA-seq protocol" TYPE [yes, no, reverse] DEFAULT no (If you select no, a read is considered overlapping with a feature regardless of whether it is mapped to the same or the opposite strand as the feature. If you select yes, the read has to be mapped to the same strand as the feature. You have to say no, if your was not made with a strand-specific RNA-seq protocol, because otherwise half your reads will be lost.)
# PARAMETER organism: "Organism" TYPE [Homo_sapiens.GRCh37.62.chr.gtf: "Human (hg19.62)", Mus_musculus.NCBIM37.62.chr.gtf: "Mouse (mm9.62)", Rattus_norvegicus.RGSC3.4.62.chr.gtf: "Rat (rn4.62)"] DEFAULT Homo_sapiens.GRCh37.62.chr.gtf (Which organism is your data from.)
# PARAMETER OPTIONAL mode: "Mode to handle reads overlapping more than one feature" TYPE [union, intersection-strict, intersection-nonempty] DEFAULT union (How to deal with reads that overlap more than one gene or exon?)
# PARAMETER OPTIONAL minaqual: "Minimum alignment quality" TYPE INTEGER FROM 0 TO 100 DEFAULT 0 (Skip all reads with alignment quality lower than the given minimum value.)
# PARAMETER OPTIONAL feature.type: "Feature type to count" TYPE [exon, CDS] DEFAULT exon (Which feature type to use, all features of other type are ignored.)
# PARAMETER OPTIONAL id.attribute: "Feature ID to use" TYPE [gene_id, transcript_id, gene_name, transcript_name, protein_name] DEFAULT gene_id (GFF attribute to be used as feature ID. Several GFF lines with the same feature ID will be considered as parts of the same feature. The feature ID is used to identity the counts in the output table.)

# TH and EK 18.1.2012
# EK changed to use Ensembl GTFs 17.4.2012

# bash wrapping
python.path <- paste(sep="", "PYTHONPATH=", file.path(chipster.tools.path, "lib", "python2.6", "site-packages"), ":$PYTHONPATH")
command.start <- paste("bash -c '", python.path, ";")
command.end <- "'"

# sort bam if the data is paired-end
samtools.binary <- file.path(chipster.tools.path, "samtools", "samtools")
samtools.sort <- ifelse(paired == "yes", paste(samtools.binary, "sort -on alignment.bam sorted-by-name"), "cat alignment.bam")

# convert bam to sam
samtools.view <- paste(samtools.binary, "view -")

# htseq-count
htseq.binary <- file.path(chipster.tools.path, "htseq", "htseq-count")
gtf <- file.path(chipster.tools.path, "genomes", organism)
htseq <- paste(htseq.binary, "-q -m", mode, "-s", stranded, "-a", minaqual, "-t", feature.type, "-i", id.attribute, "-", gtf, " > htseq-counts-out.txt")

# run
command <- paste(command.start, samtools.sort, " | ", samtools.view, " | ", htseq, command.end)
system(command)

# separate result file
system("head -n -5 htseq-counts-out.txt > htseq-counts.tsv")
system("tail -n 5 htseq-counts-out.txt > htseq-count-info.txt")

# bring in file to R environment for formating
file <- c("htseq-counts.tsv")
dat <- read.table(file, header=F, sep="\t")
names(dat) <- c("id", "count")

# write result table to output
write.table(dat, file="htseq-counts.tsv", col.names=T, quote=F, sep="\t", row.names=F)

# EOF


