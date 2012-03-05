## Install packages, and dependencies, from CRAN
install.packages(repos="http://ftp.sunet.se/pub/lang/CRAN", dependencies=TRUE, pkgs="locfit")

## Install packages from Bioconductor
source("http://bioconductor.org/biocLite.R")
# Bionductor core packages
biocLite()
# Bioconductor specific packages
biocLite("DESeq")
biocLite("RPA")

# EOF

