## Install packages, and dependencies, from CRAN
install.packages(repos="http://ftp.sunet.se/pub/lang/CRAN", dependencies=TRUE, pkgs="locfit")

install.packages(repos="http://ftp.sunet.se/pub/lang/CRAN", dependencies=TRUE, pkgs="MKmisc")
install.packages(repos="http://ftp.sunet.se/pub/lang/CRAN", dependencies=TRUE, pkgs="e1071")
install.packages(repos="http://ftp.sunet.se/pub/lang/CRAN", dependencies=TRUE, pkgs="GeneCycle")
install.packages(repos="http://ftp.sunet.se/pub/lang/CRAN", dependencies=TRUE, pkgs="fastICA")

install.packages(repos="http://ftp.sunet.se/pub/lang/CRAN", c('flexmix', 'R2HTML', 'snowfall'))


# zinba and dependencies
install.packages(repos="http://ftp.sunet.se/pub/lang/CRAN", dependencies=TRUE, c('multicore','doMC','foreach','quantreg','R.utils'))
system("wget http://zinba.googlecode.com/files/zinba_2.01.tar.gz")
install.packages("zinba_2.01.tar.gz", repos=NULL)


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

