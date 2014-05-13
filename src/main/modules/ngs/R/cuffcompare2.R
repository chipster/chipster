# TOOL cuffcompare2.R: "Compare assembly to reference using Cuffcompare" (Compare GTF of assembled transcripts to reference annotation. Isoforms are tagged as overlapping, matching or novel. You can supply your own reference GTF or use one from the server.)
# INPUT input.gtf: "Input GTF file" TYPE GTF
# INPUT OPTIONAL reference.gtf: "Reference GTF file" TYPE GTF
# OUTPUT OPTIONAL cuffcmp.refmap.tsv
# OUTPUT OPTIONAL cuffcmp.tmap.tsv
# OUTPUT OPTIONAL cuffcmp.combined.gtf
# OUTPUT OPTIONAL cuffcmp.loci.tsv
# OUTPUT OPTIONAL cuffcmp.stats.txt
# OUTPUT OPTIONAL cuffcmp.tracking.tsv
# PARAMETER chr: "Chromosome names in my GTF files look like" TYPE [chr1: "chr1", 1: "1"] DEFAULT chr1 (If you are using the reference annotations provided in Chipster, check your GTF files and choose accordingly. This option is not used if you use your own reference GTF.)
# PARAMETER OPTIONAL internalgtf: "Annotation GTF" TYPE [hg19: "Human (hg19\)", mm9: "Mouse (mm9\)", mm10: "Mouse (mm10\)", rn4: "Rat (rn4\)", Schizosaccharomyces_pombe.ASM294v2.22: "Schizosaccharomyces pombe (ASM294v2.22\)"] DEFAULT hg19 (You can use own GTF file or one of those provided on the server.)
# PARAMETER OPTIONAL r: "Ignore non-overlapping reference transcripts" TYPE [yes, no] DEFAULT no (Ignore reference transcripts that are not overlapped by any transcript in any of the GTF files. Useful for ignoring annotated transcripts that are not present in your RNA-seq samples and thus adjusting the sensitivity calculation in the accuracy report.)

# AMS 24.2.2014

# binary
cuffcompare.binary <- c(file.path(chipster.tools.path, "cufflinks2", "cuffcompare"))

cuffcompare.options <- ""

if (file.exists("reference.gtf")){
	annotation.file <- "reference.gtf"
}else{
	if (internalgtf == "hg19") {
		if (chr == 1){
			annotation.file <- "Homo_sapiens.GRCh37.68.gtf"
		}else {
			annotation.file <- "Homo_sapiens.GRCh37.68.chr.gtf"
		}		
	}
	if (internalgtf == "mm9") {
		if (chr == 1){
			annotation.file <- "Mus_musculus.NCBIM37.62.gtf"
		}else {
			annotation.file <- "Mus_musculus.NCBIM37.62.chr.gtf"
		}
	}
	if (internalgtf == "mm10") {
		if (chr == 1){
			annotation.file <- "Mus_musculus.GRCm38.68.gtf"
		}else{
			annotation.file <- "Mus_musculus.GRCm38.68.chr.gtf"
		}
	}
	if (internalgtf == "rn4") {
		if (chr == 1){
			annotation.file <- "Rattus_norvegicus.RGSC3.4.68.gtf"
		}else{
			annotation.file <- "Rattus_norvegicus.RGSC3.4.68.chr.gtf"
		}
	}
	if (internalgtf =="Schizosaccharomyces_pombe.ASM294v2.22"){
		annotation.file <- "Schizosaccharomyces_pombe.ASM294v2.22.gtf"
	}
	annotation.file <- c(file.path(chipster.tools.path, "genomes", "gtf", annotation.file))
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
