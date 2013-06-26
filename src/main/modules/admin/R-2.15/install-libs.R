## Install packages, and dependencies, from CRAN
install.packages(repos="http://ftp.sunet.se/pub/lang/CRAN", dependencies=TRUE, pkgs="locfit")

install.packages(repos="http://ftp.sunet.se/pub/lang/CRAN", dependencies=TRUE, pkgs="MKmisc")
install.packages(repos="http://ftp.sunet.se/pub/lang/CRAN", dependencies=TRUE, pkgs="e1071")
install.packages(repos="http://ftp.sunet.se/pub/lang/CRAN", dependencies=TRUE, pkgs="GeneCycle")
install.packages(repos="http://ftp.sunet.se/pub/lang/CRAN", dependencies=TRUE, pkgs="fastICA")

install.packages(repos="http://ftp.sunet.se/pub/lang/CRAN", c('flexmix', 'R2HTML', 'snowfall'))

install.packages(repos="http://ftp.sunet.se/pub/lang/CRAN", dependencies=TRUE, pkgs="png")


# zinba and dependencies
install.packages(repos="http://ftp.sunet.se/pub/lang/CRAN", dependencies=TRUE, c('multicore','doMC','foreach','quantreg','R.utils'))
system("wget http://zinba.googlecode.com/files/zinba_2.02.03.tar.gz")
install.packages("zinba_2.02.03.tar.gz", repos=NULL)


## Install packages from Bioconductor

source("http://bioconductor.org/biocLite.R")

# Bionductor core packages
biocLite()

# Bioconductor specific packages
biocLite("DESeq")
biocLite("RPA")

# hdrcde
biocLite("hdrcde")

# methylumi packages
biocLite(c("lumi", "methylumi")) # install hdrcde dependencies manually
biocLite("IlluminaHumanMethylation450k.db") # annotation package, not needed if all annotation packages from the repository are installed
biocLite("IlluminaHumanMethylation27k.db") # annotation package, not needed if all annotation packages from the repository are installed

biocLite("VariantAnnotation")
biocLite("edgeR")
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")

biocLite("DEXSeq")

biocLite("BSgenome.Hsapiens.UCSC.hg19")

biocLite("vegan")
biocLite("rich")
biocLite("BiodiversityR")
biocLite("pegas")
biocLite("labdsv")

biocLite(c('CGHregions', 'CGHcall', 'CGHbase', 'GOstats', 'impute'))


# SNP 5.0 / 6.0
biocLite("crlmm")
biocLite("oligo")
biocLite("pd.mapping50k.hind240")
biocLite("pd.mapping50k.xba240")
biocLite("pd.mapping250k.nsp")
biocLite("pd.mapping250k.sty")
biocLite("pd.genomewidesnp.5")
biocLite("pd.genomewidesnp.6")

biocLite("genomewidesnp6Crlmm")
biocLite("genomewidesnp5Crlmm")

biocLite("human1mduov3bCrlmm")          #Illumina
biocLite("human1mv1cCrlmm")             #Illumina
biocLite("human370quadv3cCrlmm")        #Illumina
biocLite("human370v1cCrlmm")            #Illumina
biocLite("human550v3bCrlmm")            #Illumina
biocLite("human610quadv1bCrlmm")        #Illumina
biocLite("human650v3aCrlmm")            #Illumina
biocLite("human660quadv1aCrlmm")        #Illumina
biocLite("humanomni1quadv1bCrlmm")      #Illumina
biocLite("humanomniexpress12v1bCrlmm")  #Illumina

biocLite("sva")


