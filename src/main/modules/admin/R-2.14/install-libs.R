## Use install utility functions
source("/opt/chipster/comp/modules/admin/R-2.14/smip.R")

## Install packages, and dependencies, from CRAN
install.packages(repos="http://ftp.sunet.se/pub/lang/CRAN", dependencies=TRUE, pkgs="locfit")

## Install packages from Bioconductor

source("http://bioconductor.org/biocLite.R")

# Bionductor core packages
biocLite()

# Bioconductor specific packages
biocLite("DESeq")
biocLite("RPA")
biocLite(c("lumi", "methylumi", "annotate")) # install hdrcde dependencies manually
biocLite("IlluminaHumanMethylation450k.db") # annotation package, not needed if all annotation packages from the repository are installed
biocLite("IlluminaHumanMethylation27k.db") # annotation package, not needed if all annotation packages from the repository are installed

# hdrcde
biocLite(c("locfit", "ash", "ks")) # install hdrcde dependencies manually
smart.install.packages(url.package="http://cran.r-project.org/src/contrib/00Archive/hdrcde/hdrcde_2.15.tar.gz") # install hdrcde manually because correct version is needed for R-2.14

# EOF

