# ANALYSIS Utilities/"Export GEO's SOFT format" (Writes out a text file in GEO's SOFT format. This file is suitable for 
# batch submission to GEO database. Before submission to GEO, you need to fill in the missing information of sample
# descriptions and the series description.)
# INPUT GENE_EXPRS normalized.tsv, GENERIC phenodata.tsv OUTPUT soft.txt 


# Writes out an empty SOFT formatted file for that is suitable for batch submission to NCBI's GEO database once filled in
#
# JTT 7.72006
# GEO SOFT validator:
# http://www.ncbi.nlm.nih.gov/projects/geo/submission/depslip.cgi?subm=0

# Renaming variables
meth<-c("blank")
file<-c("normalized.tsv")

# Reads data
dat<-read.table(file, header=T, sep="\t", row.names=1)
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Writes a blank SOFT file with obligatory fields only
if(meth=="blank" & ncol(calls)>0) {
   for(i in 1:length(phenodata$sample)) {
      write(file="soft.txt", paste("^SAMPLE = ", phenodata$sample[i], sep=""), append=T)
      write(file="soft.txt", paste("!Sample_title = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_source_name_ch1 = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_organism_ch1 = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_characteristics_ch1 = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_molecule_ch1 = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_extract_protocol = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_label_ch1 = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_label_protocol = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_hyb_protocol = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_scan_protocol = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_description = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_data_processing = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_platform_id = ", "", sep=""), append=T)
      write(file="soft.txt", paste("#ID_REF = ", "", sep=""), append=T)
      write(file="soft.txt", paste("#VALUE = ", "", sep=""), append=T)
      write(file="soft.txt", paste("#ABS_CALL = ", "the call in an absolute analysis that indicates if the transcript was present (P), absent (A), or marginal (M)" , sep=""), append=T)
      write(file="soft.txt", paste("!Sample_table_begin" , "", sep=""), append=T)
      dat2<-data.frame(row.names(dat2), dat2[,1], calls[,i])
      names(dat2)<-c("ID_REF", "VALUE", "ABS_CALL")
      write.table(dat2, "soft.txt", sep="\t", quote=F, row.names=F, col.names=T, append=T)
      write(file="soft.txt", paste("!Sample_table_end" , sep=""), append=T)
   }
   # Series description
   write(file="soft.txt", paste("^SERIES = ", phenodata$chiptype[1], "_",  substr(as.character(Sys.time()), 1, 10), "_", substr(date(), 12,19), sep=""), append=T)
   write(file="soft.txt", paste("!Series_title = ", "", sep=""), append=T)
   write(file="soft.txt", paste("!Series_type = ", "", sep=""), append=T)
   write(file="soft.txt", paste("!Series_summary = ", "", sep=""), append=T)
   write(file="soft.txt", paste("!Series_overall_design = ", "", sep=""), append=T)
   for(i in 1:length(phenodata$sample)) {
      write(file="soft.txt", paste("!Series_sample_id = ", phenodata$sample[i], sep=""), append=T)
   }
}

if(meth=="blank" & ncol(calls)==0) {
   for(i in 1:length(phenodata$sample)) {
      write(file="soft.txt", paste("^SAMPLE = ", phenodata$sample[i], sep=""), append=T)
      write(file="soft.txt", paste("!Sample_title = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_source_name_ch1 = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_organism_ch1 = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_characteristics_ch1 = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_molecule_ch1 = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_extract_protocol = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_label_ch1 = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_label_protocol = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_hyb_protocol = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_scan_protocol = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_description = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_data_processing = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_platform_id = ", "", sep=""), append=T)
      write(file="soft.txt", paste("#ID_REF = ", "", sep=""), append=T)
      write(file="soft.txt", paste("#VALUE = ", "", sep=""), append=T)
      write(file="soft.txt", paste("!Sample_table_begin" , "", sep=""), append=T)
      dat2<-data.frame(row.names(dat2), dat2[,1])
      names(dat2)<-c("ID_REF", "VALUE")
      write.table(dat2, "soft.txt", sep="\t", quote=F, row.names=F, col.names=T, append=T)
      write(file="soft.txt", paste("!Sample_table_end" , sep=""), append=T)
   }
   # Series description
   write(file="soft.txt", paste("^SERIES = ", phenodata$chiptype[1], "_",  substr(as.character(Sys.time()), 1, 10), "_", substr(date(), 12,19), sep=""), append=T)
   write(file="soft.txt", paste("!Series_title = ", "", sep=""), append=T)
   write(file="soft.txt", paste("!Series_type = ", "", sep=""), append=T)
   write(file="soft.txt", paste("!Series_summary = ", "", sep=""), append=T)
   write(file="soft.txt", paste("!Series_overall_design = ", "", sep=""), append=T)
   for(i in 1:length(phenodata$sample)) {
      write(file="soft.txt", paste("!Series_sample_id = ", phenodata$sample[i], sep=""), append=T)
   }
}
