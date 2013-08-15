# TOOL mothur-generate-count-table: mothur-generate-count-table (Generate a table of count per taxonomic group from the combined results.)
# INPUT all.grp: "groups file" TYPE GENERIC
# INPUT all.tax: "taxonomy file" TYPE GENERIC
# OUTPUT counttable.tsv: counttable.tsv
# PARAMETER cutlevel: cutlevel TYPE INTEGER FROM 0 TO 9 DEFAULT 0 (Cut levels for taxonomic names. 0 mean retain full names, e.g. Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae;Slackia;unclassified.)


# JTT 2012-11-05


# Reads the data and tabulates it
grp<-read.table("all.grp", header=F, sep="\t")
tax<-read.table("all.tax", header=F, sep="\t")
dat<-merge(tax, grp, by.x="V1", by.y="V1")
dat2<-dat
dat2$V2.y<-gsub(".[[:digit:]]{1,}.?[[:digit:]]?)", "", as.character(dat2$V2.y))

# Cutting the taxonomic names
if(cutlevel==0) {
	dat2$newnames<-dat2$V2.y
} else {
	sp<-strsplit(dat2$V2.y, ";")
	sp2<-rep(NA, nrow(dat2))
	for(i in 1:nrow(dat2)) {
		sp2[i]<-paste(sp[[i]][1:cutlevel], collapse=";")
	}
	dat2$newnames<-sp2
}

# Creating the count table
tab<-table(dat2$V2.x, dat2$newnames)
tab2<-as.data.frame.matrix(tab)

# Writing the table to disk
write.table(tab, "counttable.tsv", col.names=T, row.names=T, sep="\t", quote=FALSE)