# TOOL ensemblfetch.R: "Retrieve data for a given organism in Ensembl" (Retrieves the genomic, cDNA or protein dataset of the given species from the Ensembl data bases.)
# OUTPUT OPTIONAL ensemblfetch.fasta
# OUTPUT OPTIONAL ensemblfetch.gtf
# OUTPUT OPTIONAL ensemblfetch_species.tsv
# OUTPUT OPTIONAL ensemblfetch.log
# PARAMETER OPTIONAL species: "Species name" TYPE STRING (Then latin name of the species for which the data is retrieved. Note that you should use under score: _ in stead of the space character in the species name. For example homo_sapiens)
# PARAMETER OPTIONAL type: "Data type to retrieve" TYPE [dna: "Genomic DNA", cdna: "cDNA trasnscipts", pep: "Protein sequences", gtf: "GTF file" ] DEFAULT dna (Sequence data type to retrieve)
# PARAMETER OPTIONAL names: "List the available species names" TYPE [ yes: "All species", nonbac: "List non-bacterial species", no: No] DEFAULT "no" (List the available species names)
# PARAMETER OPTIONAL save_log: "Collect a log file about the ensemblfetch run" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the enseblfetch run.)

# K.M 28.10.2013

# enseblfecth settings
ensemblfetch.binary <- file.path("/opt/chipster/comp/modules/admin/shell", "ensemblfetch.sh ")

if ( names == "yes" || names == "nonbac" || nchar(species) < 3 ) {
  if ( names == "nonbac" ) {	
    command.to_run <- paste(ensemblfetch.binary, " -names -bacteria no > ensemblfetch_species.tsv")
    system(command.to_run)
  }	else {
    command.to_run <- paste(ensemblfetch.binary, " -names > ensemblfetch_species.tsv" )
    system(command.to_run)
  }
} else {
  command.to_run <- paste(ensemblfetch.binary, " -type ", type, " -out ensemblfetch.fasta ", species, " > ensemblfetch.log")
  system(command.to_run)	
  if ( type == "gtf" ) {
	  system("mv ensemblfetch.fasta ensemblfetch.gtf")
  }
  
}




system ("ls -l >> ensemblfetch.log")

if ( save_log == "no") {
	system ("rm -f ensemblfetch.log")
}