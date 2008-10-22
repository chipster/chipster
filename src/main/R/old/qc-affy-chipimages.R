# Affymetrix chip images reproduced from CEL-files
# JTT 12.7.

# loads the needed libraries
library(affy)

# Reads data
data<-ReadAffy()

#produces chip images and saves then as jpg files in the data folder
len<-length(data)
cwd=getwd()
for(i in 1:len) {
file<-paste(i,"jpg", sep=".");
jpeg(file.path(cwd, file), quality=100);
image(data[,i]);
dev.off();
}