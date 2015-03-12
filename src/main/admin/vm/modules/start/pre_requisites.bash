##depends:start/sys_update.bash

## Pre-requisites:

aptitude -y --without-recommends install openjdk-7-jre

# For now OpenJDK 6 is required, as is Ubuntu's default-jre,
# so we will set OpenJDK 7 as our selected choice
update-java-alternatives -s java-1.7.0-openjdk-amd64


# GNU Plot (w/o X)
# aptitude -y install gnuplot-nox # (55 packages)
aptitude -y --without-recommends install gnuplot-nox # (11 packages)

# Ghostscript
# Needed by (at least) qc-illumina.R
aptitude -y install ghostscript

## FASTX
aptitude -y install fastx-toolkit

## OpenMPI
#aptitude -y --without-recommends install openmpi-bin
