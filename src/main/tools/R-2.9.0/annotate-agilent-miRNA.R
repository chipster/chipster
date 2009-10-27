# ANALYSIS Annotation/"Agilent miRNA" (Annotates the miRNA probes, and creates a web-page of links to databases with
# additional information.)
# INPUT GENE_EXPRS normalized.tsv OUTPUT annot.html


# Agilent miRNA annotation
# JTT 20.3.2009
# modified MG 21.10.2009

#Reads the data
dat<-read.table("normalized.tsv", sep="\t", header=T)

# Extracts the identifiers
id<-as.character(rownames(dat))


# Creates an HTML page
id[id=="empty"]<-""
id2<-gsub("\\*", "", id)

write(x="<HTML>", file="annot.html", append=T) 
write(x="<BODY>", file="annot.html", append=T)
write(x="<TABLE border=1>", file="annot.html", append=T)
write(x="<CAPTION> Annotations </CAPTION>", file="annot.html", append=T)
write(x="<TR> <TH>miRBase ID</TH> <TH>MIRANDA prediction</TH> <TH>TargetScan prediction</TH></TR>", file="annot.html", append=T)

for(i in 1:length(id)) {
   mirbase<-paste("http://microrna.sanger.ac.uk/cgi-bin/sequences/mirna_entry.pl?acc=", id2[i], sep="")
   miranda<-paste("http://microrna.sanger.ac.uk/cgi-bin/targets/v5/hit_list.pl?genome_id=native&mirna_id=", id2[i], sep="")
   targetscan<-paste("http://www.targetscan.org/cgi-bin/vert_50/targetscan.cgi?mirg=", id2[i], sep="")

   write(x=paste("<TR> <TD><A HREF=", '"', mirbase, '"', ">", id[i], "</A> </TD> <TD><A HREF=", '"', miranda, '"', ">", id[i], "</A> </TD> <TD><A HREF=", '"', targetscan, '"', ">", id[i], "</A> </TD> </TR>", sep=""), file="annot.html", append=T)
}

write(x="</TABLE>", file="annot.html", append=T)
write(x="</BODY>", file="annot.html", append=T)
write(x="</HTML>", file="annot.html", append=T)      
