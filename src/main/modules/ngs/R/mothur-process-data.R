# TOOL mothur-process-data.R: "Mothur - process data" (This tool preprocesses the sequence files. Primers and adapters are removed, sequences are aligned against templates, chimeric sequences are removed, and finally, sequences are classified to OTUs.)
# INPUT all.fasta: "Input sequences in fasta format" TYPE GENERIC
# INPUT META phenodata.tsv: "Phenodata table phenodata.tsv" TYPE GENERIC
# INPUT all.qual: "Quality data" TYPE GENERIC
# OUTPUT OPTIONAL all.groups: all.groups
# OUTPUT OPTIONAL all.trim.unique.good.filter.unique.precluster.pick.silva.taxonomy: all.trim.unique.good.filter.unique.precluster.pick.silva.taxonomy
# OUTPUT OPTIONAL mothur-logfile.txt: mothur-logfile.txt
# PARAMETER primers: "primers" TYPE [forward: "forward-only", reverse: "reverse-only", both: "both", neither: "neither"] DEFAULT both (Primer direction)
# PARAMETER use.quality: "use.quality" TYPE [yes: "yes", no: "no"] DEFAULT yes (If the quality is turned off, the sequence quality file is not used at all, but other quality options are still enforced.)
# PARAMETER qaverage: "qaverage" TYPE INTEGER FROM 0 TO 40 DEFAULT 20 (Minimum average quality of the sequence. Sequences that have a lower average quality are dropped.)
# PARAMETER maxambig: "maxambig" TYPE INTEGER FROM 0 TO 10 DEFAULT 0 (Maximum number of ambiguous bases allowed in any sequence)
# PARAMETER maxhomop: "maxhomop" TYPE INTEGER FROM 0 TO 50 DEFAULT 8 (Maximum length of a homopolymere allowed in any sequence)
# PARAMETER minlength: "minlength" TYPE INTEGER FROM 0 TO 1000 DEFAULT 250 (Minimum length of an allowed sequence)
# PARAMETER pdiffs: "pdiffs" TYPE INTEGER FROM 0 TO 10 DEFAULT 0 (Maximum number of allowed differences to primer sequences)
# PARAMETER bdiffs: "bdiffs" TYPE INTEGER FROM 0 TO 10 DEFAULT 0 (Maximum number of allowed differences to barcode sequences)
# PARAMETER minlength.screen: "minlength.screen" TYPE INTEGER FROM 0 TO 1000 DEFAULT 0 ()
# PARAMETER maxlength.screen: "maxlength.screen" TYPE INTEGER FROM 0 TO 1000 DEFAULT 500 ()
# PARAMETER chimera.method: "chimera.method" TYPE [uchime: "uchime", chimeraslayer: "slayer", perseus: "perseus"] DEFAULT uchime ()

# JTT 2012-10-02
# KM all.groups tiedostoa ei löydi vaiheessa screen.seqs 
# Missä vika
# oligos-fitteräinti kommentoitu pois päältä.
# aikaisempi vaihe prepare data ei tee kunnollista phenodata-tiedostoa

# check out if the file is compressed and if so unzip it
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("all.fasta")
unzipIfGZipFile("all.qual")

# Set the path to the mothur executable
mothur.path<-c("/opt/chipster/tools/mothur/1.28.0")
mothur.exe<-c("mothur")
align.template.path<-c("/opt/chipster/tools/mothur/data")

# Reads the phenodata and prepares the oligos file for mothur
phenodata<-read.table("phenodata.tsv", header=T, sep="\t", comment.char="")
pdata<-data.frame(phenodata$oligo, phenodata$sequence)
pdata[,1]<-as.character(pdata[,1])
if(primers=="forward") {
	o<-which(pdata[,1]=="reverse")
	pdata[o,1]<-paste("#", pdata[o,1], sep="")
}
if(primers=="reverse") {
	o<-which(pdata[,1]=="forward")
	pdata[o,1]<-paste("#", pdata[o,1], sep="")
}
if(primers=="neither") {
	o<-which(pdata[,1]=="reverse")
	pdata[o,1]<-paste("#", pdata[o,1], sep="")
	o<-which(pdata[,1]=="forward")
	pdata[o,1]<-paste("#", pdata[o,1], sep="")
}

pdata$sample<-c(as.character(phenodata$sample))
write.table(pdata, "oligos", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

####################################
#
# Prepare the batch file for mothur
#
####################################

# Remove primers, linkers, etc., and control quality
unlink("batch.mth")
write("summary.seqs(fasta=all.fasta)", "batch.mth", append=T)
if(use.quality=="yes") {
	#write(paste("trim.seqs(fasta=all.fasta, qaverage=", qaverage, ", qfile=all.qual, maxambig=", maxambig, ", maxhomop=", maxhomop, ", minlength=", minlength, ", pdiffs=", pdiffs, ", bdiffs=", bdiffs, ")", sep=""), "batch.mth", append=T)
	write(paste("trim.seqs(fasta=all.fasta, oligos=oligos, qaverage=", qaverage, ", qfile=all.qual, maxambig=", maxambig, ", maxhomop=", maxhomop, ", minlength=", minlength, ", pdiffs=", pdiffs, ", bdiffs=", bdiffs, ")", sep=""), "batch.mth", append=T)
}
if(use.quality=="no") {
	write(paste("trim.seqs(fasta=all.fasta, oligos=oligos, maxambig=", maxambig, ", maxhomop=", maxhomop, ", minlength=", minlength, ", pdiffs=", pdiffs, ", bdiffs=", bdiffs, ")", sep=""), "batch.mth", append=T)
}
write("summary.seqs(fasta=all.trim.fasta)", "batch.mth", append=T)
write("unique.seqs(fasta=all.trim.fasta)", "batch.mth", append=T)
write("summary.seqs(fasta=all.trim.unique.fasta, name=all.trim.names)", "batch.mth", append=T)
# Align sequences
write(paste("align.seqs(candidate=all.trim.unique.fasta, template=", align.template.path, "/silva.bacteria.fasta)", sep=""), "batch.mth", append=T)
write("summary.seqs(fasta=all.trim.unique.align, name=all.trim.names)", "batch.mth", append=T)
write(paste("screen.seqs(fasta=all.trim.unique.align, name=all.trim.names, group=all.groups, minlength=", minlength.screen, ", maxlength=", maxlength.screen, ")", sep=""), "batch.mth", append=T)
write("summary.seqs(fasta=all.trim.unique.good.align, name=all.trim.good.names)", "batch.mth", append=T)

# Filter and precluster
write("filter.seqs(fasta=all.trim.unique.good.align, vertical=T, trump=.)", "batch.mth", append=T)
write("unique.seqs(fasta=all.trim.unique.good.filter.fasta, name=all.trim.good.names)", "batch.mth", append=T)
write("summary.seqs(fasta=current, name=current)", "batch.mth", append=T)
write("pre.cluster(fasta=all.trim.unique.good.filter.unique.fasta, name=all.trim.unique.good.filter.names, diffs=1)", "batch.mth", append=T)
write("summary.seqs(fasta=current, name=current)", "batch.mth", append=T)

# Remove chimeric sequences
if(chimera.method=="slayer") {
	write(paste("chimera.slayer(fasta=all.trim.unique.good.filter.unique.precluster.fasta, template=", align.template.path, "/silva.gold.align, blastlocation=", mothur.path, "\\blast\\bin\\", ")", sep=""), "batch.mth", append=T)
	write("remove.seqs(accnos=all.trim.unique.good.filter.unique.precluster.slayer.accnos, fasta=all.trim.unique.good.filter.unique.precluster.fasta, name=all.trim.unique.good.filter.unique.precluster.names)", "batch.mth", append=T)
}
if(chimera.method=="perseus") {
	write(paste("chimera.perseus(fasta=all.trim.unique.good.filter.unique.precluster.fasta)", sep=""), "batch.mth", append=T)
	write("remove.seqs(accnos=all.trim.unique.good.filter.unique.precluster.perseus.accnos, fasta=all.trim.unique.good.filter.unique.precluster.fasta, name=all.trim.unique.good.filter.unique.precluster.names)", "batch.mth", append=T)
}
if(chimera.method=="uchime") {
	write(paste("chimera.uchime(fasta=all.trim.unique.good.filter.unique.precluster.fasta, template=", align.template.path, "/silva.gold.align)", sep=""), "batch.mth", append=T)
}
write("summary.seqs(fasta=all.trim.unique.good.filter.unique.precluster.pick.fasta, name=all.trim.unique.good.filter.unique.precluster.pick.names)", "batch.mth", append=T) 

# Classifying sequences
write(paste("classify.seqs(fasta=all.trim.unique.good.filter.unique.precluster.pick.fasta, template=", align.template.path, "/silva.bacteria.fasta, taxonomy=", align.template.path, "/silva.bacteria.silva.tax, iters=1000)", sep=""), "batch.mth", append=T)
#remove.lineage(fasta=stool.trim.unique.good.filter.unique.precluster.pick.fasta, name=stool.trim.unique.good.filter.unique.precluster.pick.names, group=stool.good.pick.groups, taxonomy=stool.trim.unique.good.filter.unique.precluster.pick.rdp.taxonomy, taxon=Cyanobacteria)

# Ensure that mothur quits when the run is done
write("quit()", "batch.mth", append=T)

# Runs mothur 
system(paste(mothur.path, "/", mothur.exe, " batch.mth", sep=""))
system("cat *.logfile > mothur-logfile.txt ")
system("sleep 100")
# Renames the logfile
#file.rename(dir()[grep("logfile", dir())], "mothur-logfile.txt")