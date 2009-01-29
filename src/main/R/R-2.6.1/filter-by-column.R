# ANALYSIS Utilities/"Filter using a column" (Allows the user to filter the genes on the basis of one numerical column. Altered by Mandy to include discrete values)
# INPUT GENERIC normalized.tsv OUTPUT column-filter3.tsv
# PARAMETER column COLUMN_SEL (Data column to filter by)
# PARAMETER Type [Discrete, Continuous] DEFAULT Continuous (Choose Discrete or Continuous)
# PARAMETER cutoff DECIMAL FROM -100000 TO 100000 DEFAULT 1 (Cut-off for filtering)
# PARAMETER discrete.match STRING DEFAULT empty (String to search for)
# PARAMETER smaller.or.larger [equal-to, smaller-than, larger-than] DEFAULT smaller-than (Smaller or larger than the cutoff is filtered)


# Two-group parametric and non-parametric tests
# JTT 31.1.2008
# Amanda Miotto 16.1.2009 added the discrete option

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
   # Writing the data to disk
   write.table(dat2, "column-filter3.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}

if (Type=="Discrete"){
   dat2<-dat[grep(discrete.match, dat[,(as.vector(grep(column,names(dat))))]),]
   write.table(dat2, "column-filter3.tsv", sep="\t", row.names=T, col.names=T, quote=F)
}
