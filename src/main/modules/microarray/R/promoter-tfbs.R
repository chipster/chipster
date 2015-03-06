# TOOL promoter-tfbs.R: Weeder (Finds common sequence motifs in the promoters of input genes. Promoter sequences are automatically retrieved from a central database. Currently works only for human, mouse, rat, drosophila, and yeast data.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT seqs.html: seqs.html 
# OUTPUT seqs.txt.wee: seqs.txt.wee 
# PARAMETER species: Species TYPE [human: human, mouse: mouse, rat: rat, drosophila: drosophila, yeast: yeast] DEFAULT human ()
# PARAMETER promoter.size: "Promoter size" TYPE [1000: small, 2000: medium, 5000: large] DEFAULT 1000 (Length of upstream sequences. Small=1000 bp, medium=2000 bp and large=5000 bp)
# PARAMETER multiple.promoters: "Retrieve multiple promoters per gene" TYPE [yes: yes, no: no] DEFAULT no (In the case where the gene has more than one transcription start site, print them all)
# PARAMETER strands: Strands TYPE [single: single, both: both] DEFAULT single (Analyze both strands of DNA)
# PARAMETER appears.more.than.once: "Appears more than once" TYPE [yes: yes, no: no] DEFAULT no (Could the motif appear more than once in every sequence)
# PARAMETER no.of.motifs: "No of motifs" TYPE INTEGER FROM 1 TO 100 DEFAULT 10 (Number of motifs to return)
# PARAMETER percentage: Percentage TYPE INTEGER FROM 1 TO 100 DEFAULT 50 (Percentage of sequences the motif should appear)
# PARAMETER tfsize: "TFBS size" TYPE [small: small, medium: medium] DEFAULT small (Transcription factor binding site size)

# 17.11.2006 JTT
# 27.09.2013 MK

# Sets up the path to the promoter sequences
path.seq<-c(file.path(chipster.tools.path, "weeder", "seqs"))

# Sets up the paths to weeder executables
path.weeder.tfbs<-c(file.path(chipster.tools.path, "weeder", "Weeder1.4.2", "weederTFBS.out"))
path.weeder.advicer<-c(file.path(chipster.tools.path, "weeder", "Weeder1.4.2", "adviser.out"))

# Sets up the paths to weeder frequency files
path.weeder.freq<-c(file.path(chipster.tools.path, "weeder", "Weeder1.4.2", "FreqFiles"))

# Renaming variable
size<-promoter.size
once<-appears.more.than.once
no<-no.of.motifs

# Retrieving the sequences. The function is available in common/R-2.12/promoter.utils.
source(file.path(chipster.common.path, "promoter-utils.R"))
seqs <- retreive_promoters(species, promoter.size, multiple.promoters, "normalized.tsv", "phenodata.tsv")

# Write sequences on disk
for(i in 1:length(seqs$seq.names)) {
      write(file="seqs.txt", seqs$seq.names[i], append=T)
      write(file="seqs.txt", seqs$seq[i], append=T)
}

# Weeder needs frequency files - copying those in place
system("mkdir FreqFiles")
system(paste("ln -s ", file.path(path.weeder.freq, "*"), " ", file.path("FreqFiles", "."), sep=""))

# Assembling Weeder arguments
arg<-c()
if(percentage!=50) {
   arg<-paste(arg, " -R ", percentage, sep="")
} else {
   arg<-paste(arg, " -R ", 50, sep="")
}
if(strands=="both") {
   arg<-paste(arg, " -S ", sep="")
}
if(once=="yes") {
   arg<-paste(arg, "-M ", sep="")
}
if(no!=10) {
   arg<-paste(arg, "-T ", no, sep="")
}

# Running Weeder to find tentative TFBSs

if(species=="human") { weeder.org <- "HS" } 
if(species=="mouse") { weeder.org <- "MM" }
if(species=="rat") { weeder.org <- "RN" }
if(species=="yeast") { weeder.org <- "SC" }

if(tfsize=="small") {
   system(paste(path.weeder.tfbs, " -f seqs.txt -O ", weeder.org, " -W 6 -e 1", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O ", weeder.org, " -W 8 -e 2", arg, sep=""))
}
if(tfsize=="medium") {
   system(paste(path.weeder.tfbs, " -f seqs.txt -O ", weeder.org, " -W 6 -e 1", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O ", weeder.org, " -W 8 -e 2", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O ", weeder.org, " -W 10 -e 3", arg, sep=""))
}
if(tfsize=="large") {
   system(paste(path.weeder.tfbs, " -f seqs.txt -O ", weeder.org, " -W 12 -e 4", arg, sep=""))
}
if(tfsize=="extra") {
   system(paste(path.weeder.tfbs, " -f seqs.txt -O ", weeder.org, " -W 6 -e 1", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O ", weeder.org, " -W 8 -e 3", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O ", weeder.org, " -W 10 -e 4", arg, sep=""))
   system(paste(path.weeder.tfbs, " -f seqs.txt -O ", weeder.org, " -W 12 -e 4", arg, sep=""))
}
if(strands=="both") {
   system(paste(path.weeder.advicer, "  seqs.txt S",sep=""))
} else {
   system(paste(path.weeder.advicer, "  seqs.txt",sep=""))
}

# fix a broken HTML markup
system("sed 's/<\\/body><\\/html>//' seqs.txt.html > seqs.html")
