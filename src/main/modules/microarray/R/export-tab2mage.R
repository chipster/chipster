# TOOL export-tab2mage.R: "Export tab2mage format" (Writes out a text file in tab2mage format. This file is suitable for batch submission to ArrayExpress database. The scripts writes out a blank file, which you need to fill yourself. This file does not necessarily cover all experimental situations, and especially the hybridization section might need modifications. For more information, see tab2mage documentation at http: tab2mage.sourceforge.net docs index.html.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT tab2mage.txt: tab2mage.txt 


# Writes out a blank tab2mage description of the experiment
# JTT 7.8.2007

# Sets the data file name
file<-c("normalized.tsv")

# Reads the files
dat<-read.table(file, header=T, sep="\t", row.names=1)
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")


# Experiment section
write(file="tab2mage.txt", "Experiment section", append=T)
write(file="tab2mage.txt", paste("domain", "ebi.ac.uk", sep="\t"), append=T)
write(file="tab2mage.txt", paste("accession", "", sep="\t"), append=T)
write(file="tab2mage.txt", paste("quality_control", "", sep="\t"), append=T)
write(file="tab2mage.txt", paste("experiment_design_type", "", sep="\t"), append=T)
write(file="tab2mage.txt", paste("name", "", sep="\t"), append=T)
write(file="tab2mage.txt", paste("description", "", sep="\t"), append=T)
write(file="tab2mage.txt", paste("release_date", Sys.Date(), sep="\t"), append=T)
write(file="tab2mage.txt", paste("submission_date", Sys.Date(), sep="\t"), append=T)
write(file="tab2mage.txt", paste("submitter", "", sep="\t"), append=T)
write(file="tab2mage.txt", paste("organization", "", sep="\t"), append=T)
write(file="tab2mage.txt", paste("publication_title", "", sep="\t"), append=T)
write(file="tab2mage.txt", paste("authors", "", sep="\t"), append=T)
write(file="tab2mage.txt", paste("journal", "", sep="\t"), append=T)
write(file="tab2mage.txt", paste("volume", "", sep="\t"), append=T)
write(file="tab2mage.txt", paste("issue", "", sep="\t"), append=T)
write(file="tab2mage.txt", paste("pages", "", sep="\t"), append=T)
write(file="tab2mage.txt", paste("year", "", sep="\t"), append=T)
write(file="tab2mage.txt", paste("pubmed_id", "", sep="\t"), append=T)
write(file="tab2mage.txt", paste("\n"), append=T)

# Protocol section
write(file="tab2mage.txt", "# This section may be deleted, if all the protocols are already described in Arrayexpress", append=T)
write(file="tab2mage.txt", "Protocol section", append=T)
write(file="tab2mage.txt", paste("accession", "text", "name", "parameters", sep="\t"), append=T)
write(file="tab2mage.txt", paste("\n"), append=T)

# Hybridization section
write(file="tab2mage.txt", "Hybridization section", append=T)
write(file="tab2mage.txt", paste("File[raw]", "Array[accession]", "Array[serial]", "Protocol[treatment]", "Protocol[extraction]", "Protocol[labeling]", "Protocol[hybridization]", "Protocol[scanning]", "BioSource", "Sample", "Extract", "LabeledExtract", "Hybridization", "BioSourceMaterial", "SampleMaterial", "ExtractMaterial", "LabeledExtractMaterial", "Dye", "FactorValue[StrainOrLine]", "BioMaterialCharacteristics[Organism]", "BioMaterialCharacteristics[StrainOrLine]", "BioMaterialCharacteristics[Sex]", sep="\t"), append=T)
r<-c("\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t")
for(i in 1:length(phenodata$sample)) {
   write(file="tab2mage.txt", paste(phenodata$sample[i], r, sep=""), append=T)
}
