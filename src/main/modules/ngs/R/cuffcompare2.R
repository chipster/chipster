# TOOL cuffcompare2.R: "Compare assembly to reference using Cuffcompare" (Compare GTF of assembled transcripts to reference annotation. Isoforms are tagged as overlapping, matching or novel. You can supply your own reference GTF or use one from the server.)
# INPUT input.gtf: "Input GTF file" TYPE GTF
# INPUT OPTIONAL reference.gtf: "Reference GTF file" TYPE GTF
# OUTPUT OPTIONAL cuffcmp.refmap.tsv
# OUTPUT OPTIONAL cuffcmp.tmap.tsv
# OUTPUT OPTIONAL cuffcmp.combined.gtf
# OUTPUT OPTIONAL cuffcmp.loci.tsv
# OUTPUT OPTIONAL cuffcmp.stats.txt
# OUTPUT OPTIONAL cuffcmp.tracking.tsv
# PARAMETER chr: "Chromosome names in my GTF files look like" TYPE [chr1, 1] DEFAULT 1 (If you are using the reference annotations provided in Chipster, check your GTF files and choose accordingly. This option is not used if you use your own reference GTF.)
# PARAMETER OPTIONAL organism: "Reference organism" TYPE [other, Arabidopsis_thaliana.TAIR10.26, Bos_taurus.UMD3.1.79, Canis_familiaris.BROADD2.67, Canis_familiaris.CanFam3.1.79, Drosophila_melanogaster.BDGP5.78, Drosophila_melanogaster.BDGP6.79, Felis_catus.Felis_catus_6.2.79, Gallus_gallus.Galgal4.79, Gasterosteus_aculeatus.BROADS1.79, Halorubrum_lacusprofundi_atcc_49239.GCA_000022205.1.26, Homo_sapiens.GRCh37.75, Homo_sapiens.GRCh38.79, Homo_sapiens.NCBI36.54, Medicago_truncatula.GCA_000219495.2.26, Mus_musculus.GRCm38.79, Mus_musculus.NCBIM37.67, Ovis_aries.Oar_v3.1.79, Populus_trichocarpa.JGI2.0.26, Rattus_norvegicus.RGSC3.4.69, Rattus_norvegicus.Rnor_5.0.79, Schizosaccharomyces_pombe.ASM294v2.26, Sus_scrofa.Sscrofa10.2.79, Vitis_vinifera.IGGP_12x.26, Yersinia_enterocolitica_subsp_palearctica_y11.GCA_000253175.1.25] DEFAULT other (You can use own GTF file or one of those provided on the server.)
# PARAMETER OPTIONAL r: "Ignore non-overlapping reference transcripts" TYPE [yes, no] DEFAULT no (Ignore reference transcripts that are not overlapped by any transcript in any of the GTF files. Useful for ignoring annotated transcripts that are not present in your RNA-seq samples and thus adjusting the sensitivity calculation in the accuracy report.)

# AMS 24.2.2014
# AMS 2014.06.18 Changed the handling of GTF files
# AMS 04.07.2014 New genome/gtf/index locations & names

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("input.gtf")
unzipIfGZipFile("reference.gtf")

# binary
cuffcompare.binary <- c(file.path(chipster.tools.path, "cufflinks2", "cuffcompare"))


cuffcompare.options <- ""

if (organism == "other"){
	# If user has provided a GTF, we use it
	if (file.exists("reference.gtf")){
		annotation.file <- "reference.gtf"
	}else{
		stop(paste('CHIPSTER-NOTE: ', "You need provide a GTF file if you are not using one of the provided ones."))
	}
}else{
	# If not, we use the internal one.
	internal.gtf <- file.path(chipster.tools.path, "genomes", "gtf", paste(organism, ".gtf" ,sep="" ,collapse=""))
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
