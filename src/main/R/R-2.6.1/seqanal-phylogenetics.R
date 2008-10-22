# ANALYSIS "Miscellaneous"/Phylogenetics (Takes a sequence alignment in Clustal-format, and performs a phylogenetic
# analysis. Currently only distance and likelihood analyses are working.)
# INPUT GENERIC aligned-seqs.txt OUTPUT trees.txt, trees.png, log.txt
# PARAMETER method [distance, likelihood] DEFAULT distance (Which analysis method to use)
# PARAMETER ras INTEGER FROM 1 TO 10000 DEFAULT 1 (Number of random addition sequences)
# PARAMETER bootstraps INTEGER FROM 0 TO 10000 DEFAULT 10 (Number of bootstrap/jackknifing replicates)
# PARAMETER datatype [DNA, protein] DEFAULT DNA (Are sequences DNA or protein sequences)
# PARAMETER image.width INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# JTT 25.2.2008

# Sets the path to ClustalW executable
path.clustal<-c("/v/linux26_x86_64/appl/molbio/clustal/clustalw1.83.linux/clustalw")

# Sets the path to RAxML executable
path.raxml<-c("/v/linux26_x86_64/appl/molbio/raxml/RAxML-7.0.0/raxmlHPC")

# Sets the path to seqret executable
path.seqret<-c("/v/linux26_x86_64/appl/molbio/emboss/EMBOSS-5.0.0/emboss/seqret")

# Renaming variables
w<-image.width
h<-image.height

# Loads the libraries
library(ape)

# Convert the alignment to a correct format
if(method=="likelihood") {
   write("#!/bin/tcsh", file="b")
   write(paste(path.clustal, " <<EOF", sep=""), file="b", append=T)
   write("1", file="b", append=T)
   write("aligned-seqs.txt", file="b", append=T)
   write("2", file="b", append=T)
   write("9", file="b", append=T)
   write("4", file="b", append=T)
   write("1", file="b", append=T)
   write("0", file="b", append=T)
   write("aligned-seqs.txt", file="b", append=T)
   write("", file="b", append=T)
   write("", file="b", append=T)
   write("x", file="b", append=T)
   write("EOF", file="b", append=T)
   system("sh b")
}


# Sanity checks
if(method=="parsimony" | method=="bayesian") {
   # system(paste(path.seqret, " clustal::aligned-seqs.txt nexus::aligned-seqs.txt", sep=""))
   stop("These methods are not yet implemented in Chipster!")
}


# Initiates the result files
write("", file="trees.txt")
write("", file="log.txt")


# Runs the distance analysis and plots the tree
if(method=="distance") {
   # Using ClustalW for inferring the distance tree
   write("#!/bin/tcsh", file="d")
   write(paste(path.clustal, " <<EOF", sep=""), file="d", append=T)
   write("4", file="d", append=T)
   write("1", file="d", append=T)
   write("aligned-seqs.txt", file="d", append=T)
   write("3", file="d", append=T)
   write("6", file="d", append=T)
   write("5", file="d", append=T)
   write("", file="d", append=T)
   write("5", file="d", append=T)
   write("", file="d", append=T)
   write(floor(runif(1)*1000), file="d", append=T)
   write(bootstraps, file="d", append=T)
   write("", file="d", append=T)
   write("x", file="d", append=T)
   write("EOF", file="d", append=T)
   system("sh d")
   system("cp aligned-seqs.phb trees.txt")
   write("Clustal run successful!", file="log.txt")
   # Drawing the tree
   tr<-read.tree("trees.txt")
   bitmap(file="trees.png", width=w/72, height=h/72)
   plot(tr, font=1)
   add.scale.bar(x=0, y=0.5, length=round(max(tr$edge.length)/2, digits=2))
   boot.values<-tr$node.label
   boot.values[1]<-NA
   boot.values<-round(as.numeric(boot.values)/bootstraps*100)
   boot.values<-as.character(boot.values) 
   nodelabels(boot.values, frame="n", cex=1, adj=c(-0.2))
   dev.off()
}


# Runs the likelihood analysis
#if(method=="likelihood" & datatype=="DNA") {
#   for(i in 1:ras) {
#      system(paste("raxmlHPC -m GTRGAMMA -s aligned-seqs.txt -n chipster", i, sep=""))
#      system(paste("cat RAxML_result.chipster", i, " >>trees.txt", sep=""))
#      system(paste("grep Inference RAxML_info.chipster", i, " | head -n 1 >>log.txt", sep=""))
#   }
#   if(bootstraps>0) {
#      system(paste("raxmlHPC -m GTRGAMMA -s aligned-seqs.txt -n chipsterboot -# ", bootstraps, " -b ", floor(runif(1)*1000), sep=""))
#   }
#}
#
#if(method=="likelihood" & datatype=="protein") {
#   for(i in 1:ras) {
#      system(paste("raxmlHPC -m PROTGAMMAJTT -s aligned-seqs.txt -n chipster", i, sep=""))
#      system(paste("cat RAxML_result.chipster", i, " >>trees.txt", sep=""))
#      system(paste("grep Inference RAxML_info.chipster", i, " | head -n 1 >>log.txt", sep=""))
#   }
#   if(bootstraps>0) {
#      system(paste("raxmlHPC -m PROTGAMMAJTT -s aligned-seqs.txt -n chipsterboot -# ", bootstraps, " -b ", floor(runif(1)*1000), sep=""))
#   }
#}


# Runs the likelihood analysis and plots the tree
if(method=="likelihood") {
   # Running the analysis
   if(datatype=="DNA") {
      system(paste(path.raxml, " -m GTRGAMMA -s aligned-seqs.txt -n chipster -# ", ras, sep=""))
      system(paste(path.raxml, " -m GTRGAMMA -s aligned-seqs.txt -n chipsterboot -# ", bootstraps, " -b ", floor(runif(1)*1000), sep=""))
   }
   if(datatype=="protein") {
      system(paste(path.raxml, " -m GTRGAMMA -s aligned-seqs.txt -n chipster -# ", ras, sep=""))
      system(paste(path.raxml, " -m GTRGAMMA -s aligned-seqs.txt -n chipsterboot -# ", bootstraps, " -b ", floor(runif(1)*1000), sep=""))
   }
   if(ras>1) {
      system("grep Best RAxML_info.chipster > best")
      best.tree<-gsub(pattern=":", replacement="", x=scan(file="best", what="asis")[6])
      system(paste("cp ", "RAxML_result.chipster.RUN.", best.tree, " trees.txt", sep=""))
      system(paste("cp ", "RAxML_log.chipster.RUN.", best.tree, " log.txt", sep=""))
   } else {
      system(paste("cp ", "RAxML_result.chipster", " trees.txt", sep=""))
      system(paste("cp ", "RAxML_info.chipster", " log.txt", sep=""))
   }
   # system("raxmlHPC -f b -z RAxML_bootstrap.chipsterboot -t trees.txt -m GTRGAMMA -n chipsterboottree -s aligned-seqs.txt")
   # system("cat RAxML_bipartitions.chipsterboottree >>trees.txt")
   # Plotting the tree
   tr<-read.tree("trees.txt")
   tr.boot<-read.tree("RAxML_bootstrap.chipsterboot")
   boot.values<-c(prop.clades(tr, tr.boot))
   boot.values<-boot.values/bootstraps*100
   boot.values[1]<-NA
   bitmap(file="trees.png", width=w/72, height=h/72)
   plot(tr, font=1)
   add.scale.bar(x=0, y=0.5, length=round(max(tr$edge.length)/2, digits=2))
   nodelabels(boot.values, frame="n", cex=1, adj=c(-0.2))
   dev.off()
}




# if(method=="likelihood") {
#   # Running the analysis
#   if(datatype=="DNA") {
#      system(paste("phyml aligned-seqs.txt 0 i 1 ", bootstraps, " GTR e e 4 e BIONJ y y", sep=""))
#   }
#   if(datatype=="protein") {
#      system(paste("phyml aligned-seqs.txt 1 i 1 ", bootstraps, " JTT e 4 e BIONJ y y", sep=""))
#   }



# Runs the parsimony analysis
##!/bin/tcsh
#tnt <<EOF
#p aligned-seqs.txt
#taxname =
#out 0
#hold 100
#mult=tbr replic 100 hold 1;
#tsave *trees.txt; save; tsave /;
#resample=jak probability 34 frequency;
#quit;
#EOF




