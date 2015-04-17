# TOOL annotate-agilent-miRNA.R: "miRNA annotations" (Creates a web-page of links to miRNA databases for annotation information.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# OUTPUT annot.html: annot.html 

# JTT 20.3.2009
# MG 21.10.2009 Modified
# AMS 15.04.2015 Fixed links in output

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
write(x="<TR> <TH>miRBase</TH> <TH>MicroCosm predictions</TH> <TH>TargetScan predictions</TH></TR>", file="annot.html", append=T)

for(i in 1:length(id)) {
	# remove -3p, -5p from end of name
	shortid <- strsplit(id2[i], "-.p$")
	# remove also version number from end of name
	shorterid <- strsplit(as.character(shortid[1]), "-.$")
	
	mirbase<-paste("http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=", shortid[1], sep="")
	miranda<-paste("http://www.ebi.ac.uk/enright-srv/microcosm/cgi-bin/targets/v5/hit_list.pl?genome_id=native&mirna_id=", shorterid[1], sep="")
	targetscan<-paste("http://www.targetscan.org/cgi-bin/vert_50/targetscan.cgi?mirg=", shorterid[1], sep="")

	write(x=paste("<TR> <TD><A HREF=", '"', mirbase, '"', ">", shortid[1], "</A> </TD> <TD><A HREF=", '"', miranda, '"', ">", shorterid[1], "</A> </TD> <TD><A HREF=", '"', targetscan, '"', ">", shorterid[1], "</A> </TD> </TR>", sep=""), file="annot.html", append=T)
}

write(x="</TABLE>", file="annot.html", append=T)
write(x="</BODY>", file="annot.html", append=T)
write(x="</HTML>", file="annot.html", append=T)      
