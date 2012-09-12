## Install packages, and dependencies, from CRAN
install.packages(repos="http://ftp.sunet.se/pub/lang/CRAN", dependencies=TRUE, pkgs="locfit")

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
