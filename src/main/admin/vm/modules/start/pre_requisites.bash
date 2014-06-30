##depends:start/sys_update.bash

## Pre-requisites:

# Java
# sun-java6-jre (prefered would be openjdk-7-jre)
# Bypasses manual interaction to accept license
#echo 'sun-java6-jre shared/accepted-sun-dlj-v1-1 select true' | /usr/bin/debconf-set-selections
# aptitude -y install sun-java6-jre # (40 packages)
#aptitude -y --without-recommends install sun-java6-jre # (7 packages)
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
