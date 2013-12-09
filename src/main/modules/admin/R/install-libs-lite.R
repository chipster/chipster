# Lighweight version, to be used before R-3.0 is rolled completely to production.

# Configure paths and repos (change if you need)
chipster.path <- "/opt/chipster/"
repo.cran <- "http://ftp.sunet.se/pub/lang/CRAN"
repo.bioc <- "http://www.bioconductor.org"

# Use smart.* install utility functions
# They skip all packages that already have been installed
source(file(chipster.path, "comp/modules/admin/R-3.0/smip.R"))

smart.install.packages(bioconductor.package="edgeR", mirror=repo.bioc)
smart.install.packages(bioconductor.package="Heatplus", mirror=repo.bioc)
