# This install script is run in bare R installation and 
# it installs all packages required to run Chipster.
# The script uses install functions that check each package
# before installation, meaning that it can be rerun if needed
# and only missing packages are installed.

# Determine the path of the executing script
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
# Make path name absolute
script.basename <- normalizePath(script.basename)

# Use smart.* install utility functions
# They skip all packages that already have been installed
source(paste(script.basename, "/smip.R", sep=""));

# Configure paths and repos (change if you need)
repo.cran <- "http://ftp.sunet.se/pub/lang/CRAN"
repo.bioc <- "http://www.bioconductor.org"

#check where this script resides
#relative.script.dir <- dirname(parent.frame(2)$ofile)
#absolute.script.dir <- normalizePath(relative.script.dir)
#source(paste(absolute.script.dir, "/smip.R", sep=""))

# Use smart.* install utility functions
# They skip all packages that already have been installed
#source("smip.R")

#Unlike R, RScript does not seem to load the method-package, why some try-catches can crash
library(methods)

# Install packages, and their dependencies, from CRAN
smart.install.packages(package="gplots", mirror=repo.cran)

# Install packages, and their dependencies, from Bioconductor
smart.install.packages(bioconductor.package="DESeq2", mirror=repo.bioc)

smart.install.packages(bioconductor.package="illuminaHumanv1.db", mirror=repo.bioc)
smart.install.packages(bioconductor.package="illuminaHumanv2.db", mirror=repo.bioc)
smart.install.packages(bioconductor.package="illuminaHumanv3.db", mirror=repo.bioc)
smart.install.packages(bioconductor.package="illuminaHumanv4.db", mirror=repo.bioc)
smart.install.packages(bioconductor.package="illuminaMousev1.db", mirror=repo.bioc)
smart.install.packages(bioconductor.package="illuminaMousev1p1.db", mirror=repo.bioc)
smart.install.packages(bioconductor.package="illuminaMousev2.db", mirror=repo.bioc)
smart.install.packages(bioconductor.package="illuminaRatv1.db", mirror=repo.bioc)


# required by illumina annotations below
smart.install.packages(bioconductor.package="org.Hs.eg.db", mirror=repo.bioc)
smart.install.packages(bioconductor.package="org.Mm.eg.db", mirror=repo.bioc)
smart.install.packages(bioconductor.package="org.Rn.eg.db", mirror=repo.bioc)

smart.install.packages(url.package="ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaMousev1BeadID.db_1.8.0.tar.gz");
smart.install.packages(url.package="ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaMousev2BeadID.db_1.8.0.tar.gz");
smart.install.packages(url.package="ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaMousev1p1BeadID.db_1.8.0.tar.gz");
smart.install.packages(url.package="ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaRatv1BeadID.db_1.8.0.tar.gz");
smart.install.packages(url.package="ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaHumanv1BeadID.db_1.8.0.tar.gz");
smart.install.packages(url.package="ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaHumanv2BeadID.db_1.8.0.tar.gz");
smart.install.packages(url.package="ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaHumanv3BeadID.db_1.8.0.tar.gz");
smart.install.packages(url.package="ftp://ctan.uib.no/pub/bioconductor/2.7/data/annotation/src/contrib/illuminaHumanv4BeadID.db_1.8.0.tar.gz");
