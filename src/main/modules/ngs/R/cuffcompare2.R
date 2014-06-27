# TOOL cuffcompare2.R: "Compare assembly to reference using Cuffcompare" (Compare GTF of assembled transcripts to reference annotation. Isoforms are tagged as overlapping, matching or novel. You can supply your own reference GTF or use one from the server.)
# INPUT input.gtf: "Input GTF file" TYPE GTF
# INPUT OPTIONAL reference.gtf: "Reference GTF file" TYPE GTF
# OUTPUT OPTIONAL cuffcmp.refmap.tsv
# OUTPUT OPTIONAL cuffcmp.tmap.tsv
# OUTPUT OPTIONAL cuffcmp.combined.gtf
# OUTPUT OPTIONAL cuffcmp.loci.tsv
# OUTPUT OPTIONAL cuffcmp.stats.txt
# OUTPUT OPTIONAL cuffcmp.tracking.tsv
# PARAMETER chr: "Chromosome names in my GTF files look like" TYPE [chr1: "chr1", 1: "1"] DEFAULT 1 (If you are using the reference annotations provided in Chipster, check your GTF files and choose accordingly. This option is not used if you use your own reference GTF.)
# PARAMETER OPTIONAL organism: "Annotation GTF" TYPE [Homo_sapiens.GRCh37.75.gtf: "Human (hg19\)", Mus_musculus.GRCm38.75.gtf: "Mouse (mm10\)", Rattus_norvegicus.Rnor_5.0.75.gtf: "Rat (rn5\)", Schizosaccharomyces_pombe.ASM294v2.22.gtf: "Schizosaccharomyces pombe (ASM294v2.22\)"] DEFAULT Homo_sapiens.GRCh37.75.gtf (You can use own GTF file or one of those provided on the server.)
# PARAMETER OPTIONAL r: "Ignore non-overlapping reference transcripts" TYPE [yes, no] DEFAULT no (Ignore reference transcripts that are not overlapped by any transcript in any of the GTF files. Useful for ignoring annotated transcripts that are not present in your RNA-seq samples and thus adjusting the sensitivity calculation in the accuracy report.)

# AMS 24.2.2014
# AMS 2014.06.18 Changed the handling of GTF files

# binary
cuffcompare.binary <- c(file.path(chipster.tools.path, "cufflinks2", "cuffcompare"))


cuffcompare.options <- ""

if (file.exists("reference.gtf")){
	# If user has provided a GTF, we use it
	annotation.file <- "reference.gtf"
}else{
	# If not, we use the internal one.
	internal.gtf <- file.path(chipster.tools.path, "genomes", "gtf", organism)
	# If chromosome names in BAM have chr, we make a temporary copy of gtf with chr names, otherwise we use it as is.
	if(chr == "chr1"){
		source(file.path(chipster.common.path, "gtf-utils.R"))
		addChrToGtf(internal.gtf, "internal_chr.gtf") 
		annotation.file <- paste("internal_chr.gtf")
	}else{
		annotation.file <- paste(internal.gtf)
	}
}	
cuffcompare.options <-paste(cuffcompare.options, "-r", annotation.file)

if (r == "yes"){
	cuffcompare.options <-paste(cuffcompare.options, "-R")
}


# command
command <- paste(cuffcompare.binary, cuffcompare.options, "input.gtf")

# run
system(command)

# rename outputs
system("mv cuffcmp.input.gtf.refmap cuffcmp.refmap.tsv")
system("mv cuffcmp.input.gtf.tmap cuffcmp.tmap.tsv")
system("mv cuffcmp.combined cuffcmp.combined.gtf")
system("mv cuffcmp.loci cuffcmp.loci.tsv")
system("mv cuffcmp.stats cuffcmp.stats.txt")
system("mv cuffcmp.tracking cuffcmp.tracking.tsv")
