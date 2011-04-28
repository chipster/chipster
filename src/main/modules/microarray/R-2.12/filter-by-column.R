# TOOL filter-by-column.R: "Filter using a column" (Allows the user to filter the genes on the basis of one numerical column. Altered by Mandy to include discrete values)
# INPUT normalized.tsv: normalized.tsv TYPE GENERIC 
# OUTPUT column-filter3.tsv: column-filter3.tsv 
# PARAMETER column: column TYPE COLUMN_SEL (Data column to filter by)
# PARAMETER Type: Type TYPE [Discrete: Discrete, Continuous: Continuous] DEFAULT Continuous (Choose Discrete or Continuous)
# PARAMETER cutoff: cutoff TYPE DECIMAL FROM -100000 TO 100000 DEFAULT 1 (Cut-off for filtering)
# PARAMETER discrete.match: discrete.match TYPE STRING DEFAULT empty (String to search for)
# PARAMETER smaller.or.larger: smaller.or.larger TYPE [equal-to: equal-to, smaller-than: smaller-than, larger-than: larger-than, within: within, outside: outside] DEFAULT smaller-than (Smaller or larger than the cutoff is filtered. Use the within or outside options to filter symmmetrically around two cut-offs, useful for example when searching for 2-fold up- and down-regulated genes.)

# Filtering using a column
# JTT 31.1.2008
# Amanda Miotto 16.1.2009 added the discrete option
#
# MG, 2.3.2010 added the option to filter "within" and "outside" a symmetrical
# range of values

# Parameter settings (default) for testing purposes
#column<-c("chip.microarray001.cel")
#Type<-c("Continuous")
#cutoff<-c(7)
#discrete.match<-c("empty")
#smaller.or.larger<-c("larger-than")


# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

if(Type=="Continuous") {
	# Extract the data to a vector
	f<-dat[,grep(column, colnames(dat))]
	# Filters the data
	if(smaller.or.larger=="equal-to") {
		dat2<-dat[which(f==cutoff),]
	}
	if(smaller.or.larger=="smaller-than") {
		dat2<-dat[which(f<=cutoff),]
	}
	if(smaller.or.larger=="larger-than") {
		dat2<-dat[which(f>=cutoff),]
	}
	if(smaller.or.larger=="outside") {
		cutoff_2 <- -cutoff
		dat2<-dat[which(f>=cutoff | f<=cutoff_2),]
	}
	if(smaller.or.larger=="within") {
		cutoff_2 <- -cutoff
		dat2<-dat[-which(f>=cutoff | f <=cutoff_2),]
	}
	
	# Writing the data to disk
	write.table(dat2, "column-filter3.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}

if (Type=="Discrete"){
	dat2<-dat[grep(discrete.match, dat[,(as.vector(grep(column,names(dat))))]),]
	write.table(dat2, "column-filter3.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}
