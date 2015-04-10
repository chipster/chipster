# TOOL htseq-count.R: "Count aligned reads per genes with HTSeq" (Calculates how many reads in a BAM file map to each gene. If you would like to map reads against your own GTF files, please use the tool \"Count aligned reads per genes with HTSeq using own GTF\". This tool is based on the HTSeq package. In order to use the output in edgeR or DESeq, you need to select all samples and run the tool \"Utilities - Define NGS experiment\".)
# INPUT alignment.bam: "BAM alignment file" TYPE GENERIC
# OUTPUT htseq-counts.tsv
# OUTPUT OPTIONAL htseq-count-info.txt
# PARAMETER organism: "Reference organism" TYPE [Arabidopsis_thaliana.TAIR10.26, Bos_taurus.UMD3.1.79, Canis_familiaris.BROADD2.67, Canis_familiaris.CanFam3.1.79, Drosophila_melanogaster.BDGP5.78, Drosophila_melanogaster.BDGP6.79, Felis_catus.Felis_catus_6.2.79, Gallus_gallus.Galgal4.79, Gasterosteus_aculeatus.BROADS1.79, Halorubrum_lacusprofundi_atcc_49239.GCA_000022205.1.26, Homo_sapiens.GRCh37.75, Homo_sapiens.GRCh38.79, Homo_sapiens.NCBI36.54, Medicago_truncatula.GCA_000219495.2.26, Mus_musculus.GRCm38.79, Mus_musculus.NCBIM37.67, Ovis_aries.Oar_v3.1.79, Populus_trichocarpa.JGI2.0.26, Rattus_norvegicus.RGSC3.4.69, Rattus_norvegicus.Rnor_5.0.79, Schizosaccharomyces_pombe.ASM294v2.26, Sus_scrofa.Sscrofa10.2.79, Vitis_vinifera.IGGP_12x.26, Yersinia_enterocolitica_subsp_palearctica_y11.GCA_000253175.1.25] DEFAULT Homo_sapiens.GRCh38.79 (Which organism is your data from.)
# PARAMETER chr: "Chromosome names in the BAM file look like" TYPE [chr1, 1] DEFAULT 1 (Chromosome names must match in the BAM file and in the reference annotation. Check your BAM and choose accordingly.)
# PARAMETER paired: "Does the BAM file contain paired-end data" TYPE [yes, no] DEFAULT no (Does the alignment data contain paired end or single end reads?)
# PARAMETER stranded: "Was the data produced with a strand-specific protocol" TYPE [yes, no, reverse] DEFAULT no (Select no if your data was not produced with a strand-specific RNA-seq protocol, so that a read is considered overlapping with a feature regardless of whether it is mapped to the same or the opposite strand as the feature. If you select yes, the read has to be mapped to the same strand as the feature.)
# PARAMETER OPTIONAL mode: "Mode to handle reads overlapping more than one feature" TYPE [union, intersection-strict, intersection-nonempty] DEFAULT union (How to deal with reads that overlap more than one gene or exon?)
# PARAMETER OPTIONAL minaqual: "Minimum alignment quality" TYPE INTEGER FROM 0 TO 100 DEFAULT 10 (Skip all reads with alignment quality lower than the given minimum value.)
# PARAMETER OPTIONAL feature.type: "Feature type to count" TYPE [exon, CDS] DEFAULT exon (Which feature type to use, all features of other type are ignored.)
# PARAMETER OPTIONAL id.attribute: "Feature ID to use" TYPE [gene_id, transcript_id, gene_name, transcript_name, protein_name] DEFAULT gene_id (GFF attribute to be used as feature ID. Several GFF lines with the same feature ID will be considered as parts of the same feature. The feature ID is used to identify the counts in the output table.)
# PARAMETER OPTIONAL print.coord: "Add chromosomal coordinates to the count table" TYPE [yes, no] DEFAULT yes (If you select yes, chromosomal coordinates are added to the output file. Given are the minimum and maximum coordinates of features, e.g. exons, associated with a given identifier)

# 18.1.2012 TH and EK 
# 17.4.2012 EK changed to use Ensembl GTFs 
# 3.2.2013 AMS added chr/nochr option
# 6.5.2013 MK added chr-location information to the output
# 30.5.2013 EK changed the default for "add chromosomal coordinates" to no
# 21.5.2014 EK updated to use HTSeq 0.6.1
# 19.6.2014 AMS changed handling of GTFs
# AMS 04.07.2014 New genome/gtf/index locations & names

# bash wrapping
python.path <- paste(sep="", "PYTHONPATH=", file.path(chipster.tools.path, "lib", "python2.7", "site-packages"), ":$PYTHONPATH")
command.start <- paste("bash -c '", python.path, ";")
command.end <- "'"

# sort bam if the data is paired-end
samtools.binary <- file.path(chipster.tools.path, "samtools", "samtools")
if(paired == "yes"){
	system(paste(samtools.binary, "sort -n alignment.bam name-sorted"))
	bam<-"name-sorted.bam"
} else {
	bam<-"alignment.bam"
}

# htseq-count
if(print.coord == "no") {
	htseq.binary <- file.path(chipster.tools.path, "htseq", "htseq-count")
} else {
	htseq.binary <- file.path(chipster.tools.path, "htseq", "htseq-count_chr")
}


internal.gtf <- file.path(chipster.tools.path, "genomes", "gtf", paste(organism, ".gtf" ,sep="" ,collapse=""))
if(chr == "1"){
	annotation.file <- paste(internal.gtf)
}else{
	source(file.path(chipster.common.path, "gtf-utils.R"))
	addChrToGtf(internal.gtf, "internal_chr.gtf") 
	annotation.file <- paste("internal_chr.gtf")
}


htseq <- paste(htseq.binary, "-f bam -q -m", mode, "-s", stranded, "-a", minaqual, "-t", feature.type, "-i", id.attribute, bam, annotation.file, " > htseq-counts-out.txt")

htseq.command <- paste(command.start, htseq, command.end)
system(htseq.command)

# separate result file
system("head -n -5 htseq-counts-out.txt > htseq-counts.tsv")
system("tail -n 5 htseq-counts-out.txt > htseq-count-info.txt")

# bring in file to R environment for formating
file <- c("htseq-counts.tsv")
dat <- read.table(file, header=F, sep="\t")

if(print.coord == "no") {
	names(dat) <- c("id", "count")
} else {
	names(dat) <- c("id", "chr", "start", "end", "len", "strand", "count")
}

# write result table to output
write.table(dat, file="htseq-counts.tsv", col.names=T, quote=F, sep="\t", row.names=F)

# EOF


