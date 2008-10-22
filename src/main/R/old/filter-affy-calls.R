# ANALYSIS "Expression"/"Flag filtering" (P/M/A flag filtering)
# INPUT GENE_EXPRS calls.txt OUTPUT call-filtered.txt

# JTT, 13.7.2005
# heavily modified by JTT, 30.11.2005
# a few fixes by JTT, 1.12.2005

# Loads the normalized data
# This needs a parameter "file" that tells which data file should be loaded.
# Datafile (calls.txt) is produced by MAS5 normalization script (norm-affy.mas5.R)
# At the moment datafile can be only be "calls.txt"
file<-c("calls.txt")
data<-read.table(file, header=T, sep="\t")

# Older filter, removed 30.11.2005
# Filter
# Finds genes that are present in all samples
# len<-length(data)
# By which call are we filtering? - This is needed as a parameter f
# f<-c("P")
# for(i in 1:len) {
# data[,i]<-ifelse(data[,i]==f, as.character(data[,i]), NA);
# }
# data.P<-na.omit(data)

# This terrible and slow two-loop thingy filters the data by affy-calls
# Needs optimization, is too slow (for 22000+ genes takes several minutes)
# OK, the latter sloop was replace with row-wise apply function, this is now very much faster (takes about a second)

# Finds genes that are present/absent/marginal in at least p samples
len<-length(data)
# By which call are we filtering? - This is needed as a parameter f
f<-c("P")
# In how many sample should the gene have a call of f? - This is needed as a parameter p
# p<-c(0.50)
p<-c(3) # This is 50% of the demodata
# An assessory variable, the number of chips is rounded down
# pp<-floor(len*p)
# Find and replace all the samples with flag f with 1, otherwise recode them with 0
for(i in 1:len) {
data[,i]<-ifelse(data[,i]==f, 1, 0);
}
# len2<-dim(data)[[1]]
# s<-rep(0,len2)
# Calculate a row-wise sum of recoded values
s<-apply(data, MARGIN=1, FUN="sum")
# for(i in 1:len2) {
# s[i]<-sum(data[i,]);
# }
# Replace all the rows with a sum lower than p with NAs
s2<-ifelse(s<=p, NA, s)
# Bind data and variable s2
data2<-data.frame(data, s2)
# Remove rown with missing data, and save row names (gene names) into a variable data.calls
data.calls<-row.names(na.omit(data2))
# data.calls2<-data.calls[-(len+1)]
# Saving the results into a text file
write.table(data.calls, "call-filtered.txt", sep="\t", row.names=T, col.names=T)

